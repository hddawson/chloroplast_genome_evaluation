library(arrow)
library(data.table)
library(Biostrings)

# ---------------------------------------------------------------------
# LOAD COMMON DATA
# ---------------------------------------------------------------------

data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")

# Prepare PC table
pcs_IDS <- aln_index$index
scores <- as.data.table(ev_pcs$x)
setnames(scores, paste0("PC", seq_len(ncol(scores))))
scores[, ID := pcs_IDS]

n_pcs <- 1000L
pc_names <- paste0("PC", seq_len(n_pcs))
scores <- scores[, c("ID", pc_names), with = FALSE]

# Pre-join phenotype + PCs once
pheno_pcs <- data[, .(ID, pheno = get(pheno_col))][scores, on = "ID", nomatch = 0]
stopifnot(nrow(pheno_pcs) > 0)

# ---------------------------------------------------------------------
# PREPARE FILE LISTS
# ---------------------------------------------------------------------

emb_files <- list.files("data/embeddings/", pattern = "_residue_embeddings\\.parquet$", full.names = TRUE)
aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

get_gene <- function(path) sub("_residue_embeddings\\.parquet|_AA_aligned\\.fasta", "", basename(path))
emb_genes <- get_gene(emb_files)
aln_genes <- get_gene(aln_files)
genes_to_process <- intersect(emb_genes, aln_genes)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

all_results_list <- vector("list", length(genes_to_process))

# ---------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------

set.seed(123)

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  
  out_path <- paste0("results/embStandalone/tmp_results_nested_", gene, ".rds")
  
  if (file.exists(out_path)) {
    message("Skipping ", gene, " (output already exists)")
    next
  }
  
  emb_file <- file.path("data/embeddings/", paste0(gene, "_residue_embeddings.parquet"))
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # --- Read embeddings ---
  df_emb <- tryCatch(as.data.table(read_parquet(emb_file)), error = function(e) NULL)
  if (is.null(df_emb)) {
    message("Skipping ", gene, ": could not read embeddings file")
    next
  }
  
  # --- Read alignment ---
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # --- Strip IDs ---
  names(aln) <- sub("\\|.*", "", names(aln))
  df_emb[, ID := sub("\\|.*", "", ID)]
  
  # -------------------------------------------------------------------
  # BUILD RESIDUE LOOKUP (data.table-friendly)
  # -------------------------------------------------------------------
  
  aln_mat <- as.matrix(aln)
  seq_ids <- names(aln)
  
  # Flatten alignment into long form (fast)
  residue_lookup <- rbindlist(
    lapply(seq_along(seq_ids), function(i) {
      residues <- aln_mat[i, ]
      non_gap <- which(residues != "-")
      data.table(
        ID = seq_ids[i],
        Residue_Index = seq_len(length(non_gap)),
        Aligned_Position = non_gap
      )
    }),
    use.names = TRUE
  )
  
  # -------------------------------------------------------------------
  # JOIN EMBEDDINGS + ALIGNMENT + PHENOTYPE+PCs
  # -------------------------------------------------------------------
  
  setkey(df_emb, ID, Residue_Index)
  setkey(residue_lookup, ID, Residue_Index)
  df_joined <- residue_lookup[df_emb, nomatch = 0]
  
  # Add phenotypes & PCs
  setkey(pheno_pcs, ID)
  df_joined <- df_joined[pheno_pcs, nomatch = 0]
  
  if (nrow(df_joined) == 0) {
    message("Skipping ", gene, ": no joined data")
    next
  }
  
  
  emb_cols <- grep("^embedding_", names(df_joined), value = TRUE)
  stopifnot(length(emb_cols) > 0)
  
  # -------------------------------------------------------------------
  # SCALE PCs ONCE PER GENE (using base scale() for numeric speed)
  # -------------------------------------------------------------------
  
  common_ids <- unique(df_joined$ID)
  pheno_pcs_sub <- pheno_pcs[ID %in% common_ids]
  X_pcs_scaled <- scale(as.matrix(pheno_pcs_sub[, ..pc_names]))
  rownames(X_pcs_scaled) <- pheno_pcs_sub$ID
  
  # -------------------------------------------------------------------
  # LOOP OVER POSITIONS
  # -------------------------------------------------------------------
  
  positions <- sort(unique(df_joined$Aligned_Position))
  results_list <- vector("list", length(positions))
  
  for (j in seq_along(positions)) {
    
    
    pos <- positions[j]
    sub <- df_joined[Aligned_Position == pos]
    
    # Check sample size
    if (nrow(sub) < 8500L) next
    
    #this particular residue breaks my linalg solver
    if (pos == 19 & gene=="psbJ") next
    
    # Check residue variation
    residues <- aln_mat[match(sub$ID, names(aln)), pos]
    residue_table <- table(residues)
    if (length(residue_table) < 2L) next
    
    y <- sub$pheno
    X_pcs <- X_pcs_scaled[sub$ID, , drop = FALSE]
    X_emb <- scale(as.matrix(sub[, ..emb_cols]))
    
    # Skip if bad data before PCA
    if (any(!is.finite(X_emb))) {
      message("Skipping ", gene, " pos=", pos, " — non-finite values in X_emb")
      next
    }
    
    pca <- tryCatch({
      prcomp(X_emb, rank. = 10, scale. = TRUE)
    }, error = function(e) {
      message("Skipping ", gene, " pos=", pos, " — PCA failed: ", conditionMessage(e))
      return(NULL)
    })
    
    pv <- pca$sdev^2 / sum(pca$sdev^2)
    X_emb <- scale(as.matrix(pca$x))
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    var0 <- which(apply(X_aa, 2, var) == 0)
    if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # -----------------------------------------------------------------
    # FIT THREE MODELS (tryCatch for safety)
    # -----------------------------------------------------------------
    
    fit_reduced <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    fit_partial <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
    if (is.null(fit_partial)) next
    
    fit_full <- tryCatch(lm(y ~ X_emb + X_pcs), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_partial <- summary(fit_partial)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    p_res <- tryCatch(anova(fit_reduced, fit_partial)[2, "Pr(>F)"], error = function(e) NA_real_)
    p_emb <- tryCatch(anova(fit_partial, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    results_list[[j]] <- data.table(
      Gene = gene,
      Aligned_Position = pos,
      N = nrow(sub),
      R2_reduced = r2_reduced,
      R2_res = r2_partial,
      R2_emb = r2_full,
      P_res = p_res,
      P_emb = p_emb
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (pos %% 50 ==0) {cat(pos)}
  
  if (length(results_list) > 0) {
    all_results_list[[gene]] <- rbindlist(results_list)
    saveRDS(results_list, paste0("results/embStandalone/tmp_results_nested_", gene, ".rds"))
    message("Completed ", gene, ": ", length(results_list), " positions analyzed")
  } else {
    message("No positions retained for gene ", gene)
  }
}

saveRDS(all_results_list, "results/embStandalone/results_list_nested_all_genes.rds")
message("\nSaved nested model results for ", length(all_results_list), " genes")


gene <- "atpA"
res_list <- readRDS(paste0("results/tmp_results_nested_", gene, ".rds"))
res_df <- do.call(rbind, res_list)
summary(res_df)

pairs(res_df[,-1])

plot(-log10(res_df$P_res))
plot(-log10(res_df$P_emb))
plot(-log10(res_df$P_res),-log10(res_df$P_emb))
text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)

gene <- "atpA"
res_list <- readRDS(paste0("results/tmp_results_nested_", gene, ".rds"))
res_df <- do.call(rbind, res_list)
summary(res_df)

pairs(res_df[,-1])

hist(-log10(res_df$P_res))
hist(-log10(res_df$P_emb))
text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
hist(res_df$R2_full/res_df$R2_partial)
hist(res_df$R2_partial/res_df$R2_reduced)
hist(res_df$R2_full/res_df$R2_reduced)
plot(res_df$R2_full/res_df$R2_reduced, res_df$R2_full/res_df$R2_partial)
plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene)
lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")

gene <- "atpB"
res_list <- readRDS(paste0("results/tmp_results_nested_", gene, ".rds"))
res_df <- do.call(rbind, res_list)
summary(res_df)

hist(-log10(res_df$P_res))
hist(-log10(res_df$P_emb))
plot(-log10(res_df$P_res),-log10(res_df$P_emb))

text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
hist(res_df$R2_full/res_df$R2_partial)
hist(res_df$R2_partial/res_df$R2_reduced)
hist(res_df$R2_full/res_df$R2_reduced)
plot(res_df$R2_full/res_df$R2_reduced, res_df$R2_full/res_df$R2_partial)
plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene)
lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")

gene <- "atpB"
res_list <- readRDS(paste0("results/tmp_results_nested_", gene, ".rds"))
res_df <- do.call(rbind, res_list)
summary(res_df)

hist(-log10(res_df$P_res))
hist(-log10(res_df$P_emb))
plot(-log10(res_df$P_res),-log10(res_df$P_emb))

text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
hist(res_df$R2_full/res_df$R2_partial)
hist(res_df$R2_partial/res_df$R2_reduced)
hist(res_df$R2_full/res_df$R2_reduced)
plot(res_df$R2_full/res_df$R2_reduced, res_df$R2_full/res_df$R2_partial)
plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene)
lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")

gene <- "atpB"
res_list <- readRDS(paste0("results/tmp_results_nested_", gene, ".rds"))
res_df <- do.call(rbind, res_list)
summary(res_df)
plot(-log10(res_df$P_res),-log10(res_df$P_emb))
text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene)
lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")

res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

for (f in res_files) {
  res_list <- readRDS(f)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df <- do.call(rbind, res_list)
  res_df <- na.omit(res_df)
  #plot(-log10(res_df$P_res),-log10(res_df$P_emb))
  #text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
  plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene_name,
       ylim=c(0,2+max(c(-log10(res_df$P_res), -log10(res_df$P_emb)))),
       ylab="-log10(p)", xlab="Alignment Index")
  points(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
  lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")
  #bf <- 0.05 / (2 * max(res_df$Aligned_Position))
  abline(h=-log10(bf), col="black", lty=2)
  
  points(res_df[-log10(res_df$P_res) > -log10(bf),]$Aligned_Position,
         -log10(res_df$P_res[-log10(res_df$P_res) > -log10(bf)]), col="orange",pch=16)

  
  legend(
    "topleft",
    legend = c("Residue p (residues | PCs)",
                   "Embedding p (embeddings | residues + PCs)",
                   "Bonferroni threshold"),
    col = c("orange", "tomato", "black"),
    lty = c(1, 1, 2),
    lwd = c(2, 2, 1),
    bty = "y",  # remove box
    cex = 0.8
  )
}

plot(res_df$R2_partial)


library(ggplot2)

all_res <- lapply(res_files, function(f) {
  res_list <- readRDS(f)
  res_df <- do.call(rbind, res_list)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df$gene <- gene_name
  res_df
}) |> do.call(rbind, args = _)

# Melt into long format for ggplot
all_long <- tidyr::pivot_longer(all_res,
                                cols = c(P_res, P_emb),
                                names_to = "type",
                                values_to = "pval")

library(tidyverse)

res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

# Load and combine all results
all_res <- lapply(res_files, function(f) {
  res_list <- readRDS(f)
  res_df <- do.call(rbind, res_list)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df$gene <- gene_name
  res_df
}) |> bind_rows()

# Sort genes by name (optional)
all_res <- all_res %>% arrange(gene, Aligned_Position)

# Assign cumulative positions along x-axis
gene_lengths <- all_res %>%
  group_by(gene) %>%
  summarise(len = max(Aligned_Position, na.rm = TRUE))

gene_offsets <- c(0, cumsum(head(gene_lengths$len, -1)))
names(gene_offsets) <- gene_lengths$gene

all_res <- all_res %>%
  mutate(pos_global = Aligned_Position + gene_offsets[gene])

# Bonferroni threshold (using total sites)
#bf <- 0.05 / (2 * sum(gene_lengths$len))
bf <- 0.05 / (1 * sum(gene_lengths$len))
ylim_max <- 1 + max(-log10(c(all_res$P_res, all_res$P_emb)), na.rm = TRUE)
ylim_max <- 1 + max(-log10(c(all_res$P_res)), na.rm = TRUE)

# Plot setup
plot(all_res$pos_global, -log10(all_res$P_res),
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "-log10(p)",
     ylim = c(0, ylim_max),
     main = "Chloroplast Protein Association Study",
     xaxt = "n")

# Alternate colors for genes
gene_colors <- setNames(rep(c("steelblue3", "grey60"),
                            length.out = length(unique(all_res$gene))),
                        unique(all_res$gene))


# Draw lines per gene
for (g in unique(all_res$gene)) {
  sub <- subset(all_res, gene == g)
  points(sub$pos_global, -log10(sub$P_res), col=gene_colors[g], lwd=1.5)
  #lines(sub$pos_global, -log10(sub$P_emb), col=adjustcolor(gene_colors[g], alpha.f=0.5), lty=2)
  points(sub[-log10(sub$P_res) > -log10(bf),]$pos_global,
         -log10(sub$P_res[-log10(sub$P_res) > -log10(bf)]), col=gene_colors[g],pch=16)
}

# Add Bonferroni threshold line
abline(h = -log10(bf), col = "black", lty = 2)

# Add gene labels at midpoints
midpoints <- gene_lengths %>%
  mutate(offset = gene_offsets,
         mid = offset + len / 2)
axis(1, at = midpoints$mid, labels = midpoints$gene, las = 2, cex.axis = 0.7)




library(jsonlite)

# --- Load UniProt JSON ---
j <- fromJSON("data/psbA_TM_uniprot.json")

# Extract start, end for Transmembrane regions
tm_uniprot <- j$features %>%
  filter(type == "Transmembrane") %>%
  transmute(
    start = location$start$value,
    end   = location$end$value
  )
# --- Map to your alignment ---
at_ID <- "AP000423.1"
at_df <- df_joined %>%
  filter(ID == at_ID) %>%
  select(Residue_Index, Aligned_Position)

at_start <- dplyr::rename(at_df, start_aln = Aligned_Position)
at_end   <- dplyr::rename(at_df, end_aln = Aligned_Position)

tm_df <- tm_uniprot %>%
  left_join(at_start, by = c("start" = "Residue_Index")) %>%
  left_join(at_end,   by = c("end"   = "Residue_Index")) %>%
  filter(!is.na(start_aln), !is.na(end_aln)) %>%
  mutate(Label = "TM helix")

res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

res_list <- readRDS("results/tmp_results_nested_psbA.rds")
gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
res_df <- do.call(rbind, res_list)
  

plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main="psbA residue associations",
     ylim=c(0,1+max(c(-log10(res_df$P_res), -log10(res_df$P_emb)))),
     ylab="-log10(p)", xlab="Alignment Index")
lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")
bf <- 0.05 / (2 * max(res_df$Aligned_Position))
abline(h=-log10(bf), col="black", lty=2)

points(res_df[-log10(res_df$P_res) > -log10(bf),]$Aligned_Position,
       -log10(res_df$P_res[-log10(res_df$P_res) > -log10(bf)]), col="orange",pch=16)
#the big one is at position 155 in the protein, 207 in aln T-->A
#the second one is position 243, 303 in aln



legend(
  "topleft",
  legend = c("Residue effect (residues | PCs)",
             "Gene Bonferroni threshold"),
  col = c("orange", "black"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 1),
  bty = "n",  # remove box
  cex = 0.8
)

e_prior <- 0
for (i in 1:nrow(tm_df)) {
  s <- tm_df[i,]$start_aln
  e <- tm_df[i,]$end_aln
  abline(v=s, col="coral")
  abline(v=e, col="coral")
  #sset <- a[a$Aligned_Position <= s &a$Aligned_Position >=e_prior,]
  #m <- mean(sset$LogLik_Ratio)
  #segments(e_prior, m, s, m, col="tomato")
  e_prior <-e
}
s <- e_prior
sset <- a[a$Aligned_Position <= max(a$Aligned_Position) & a$Aligned_Position >=s,]
m <- mean(sset$CVMRatio)
segments(e_prior, m, max(a$Aligned_Position), m, col="tomato")
text(305, 6, "DE-loop", col="tomato")

feat_json <- 'data/psbA_features.json'
feats <- as.data.frame(fromJSON(feat_json) ) 
feat_df <- feats %>%
  transmute(Residue_Index = features.location$start$value,
            Label = features.ligand$name) %>%
  left_join(at_df, by="Residue_Index") %>%
  filter(!is.na(Aligned_Position))

for (i in 1:nrow(feat_df)) {
  abline(v = feat_df[i,]$Aligned_Position, col="green", lty=2, lwd=0.6)
  text(feat_df[i,]$Aligned_Position, 1.01, feat_df[i,]$Label, 
       srt=90, adj=c(0,0), cex=0.8, col="black")
}

legend(
  "topleft",
  legend = c("Residue effect (residues | PCs)",
             "Embedding effect (embeddings | residues + PCs)",
             "Bonferroni threshold"),
  col = c("orange", "tomato", "black"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 1),
  bty = "n",  # remove box
  cex = 0.8
)


### align results to the residue index of the chloroplast proteins 
#data$Organism[grep("Zea mays", data$Organism)]


# Reference ID for mapping
at_ID <- "ON012183.1"

# Get result files
res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)
stopifnot(length(res_files) > 0)

# Prepare output list
mapped_results_list <- vector("list", length(res_files))

for (i in seq_along(res_files)) {
  f <- res_files[i]
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  
  message("Processing gene: ", gene_name)
  
  # Load results
  res_list <- readRDS(f)
  res_df <- do.call(rbind, res_list)
  
  # Load alignment for this gene
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene_name, "_AA_aligned.fasta"))
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  
  if (is.null(aln)) {
    message("  Skipping ", gene_name, ": could not read alignment")
    next
  }
  
  # Strip IDs
  names(aln) <- sub("\\|.*", "", names(aln))
  
  # Check if reference ID exists in alignment
  if (!at_ID %in% names(aln)) {
    message("  Skipping ", gene_name, ": reference ID not found in alignment")
    next
  }
  
  # Build residue index mapping for reference sequence
  aln_mat <- as.matrix(aln)
  ref_seq <- aln_mat[at_ID, ]
  non_gap <- which(ref_seq != "-")
  
  at_df <- data.table(
    Aligned_Position = non_gap,
    Residue_Index = seq_along(non_gap),
    Residue = ref_seq[non_gap]
  )
  
  # Join results with reference mapping
  setkey(at_df, Aligned_Position)
  setkey(res_df, Aligned_Position)
  
  mapped <- at_df[res_df, nomatch = 0]
  
  if (nrow(mapped) == 0) {
    message("  Skipping ", gene_name, ": no mapped positions")
    next
  }
  
  # Select and order columns
  mapped_results_list[[i]] <- mapped[, .(
    Gene = gene_name,
    Residue_Index,
    Aligned_Position,
    Residue,
    N,
    R2_reduced,
    R2_partial,
    R2_full,
    P_res,
    P_emb
  )]
  
  message("  Mapped ", nrow(mapped), " positions")
}

# Combine all genes
mapped_results_list <- Filter(Negate(is.null), mapped_results_list)
stopifnot(length(mapped_results_list) > 0)

all_mapped <- rbindlist(mapped_results_list)

# Save as CSV
out_file <- "/workdir/hdd29/chloroplast_genome_evaluation/data/speciesWork/Salix/Salix_reference_mapped_results.csv"
fwrite(all_mapped, out_file)
message("\nSaved reference-mapped results to: ", out_file)
message("Total rows: ", nrow(all_mapped))

quantile(all_mapped$P_res,0.1)

qqplot(-log10(ppoints(length(all_mapped$P_res))),
       -log10(sort(all_mapped$P_res)),
       xlab="Expected -log10(p)",
       ylab="Observed -log10(p)",
       main="Chloroplast residue QQ plot")
abline(0,1,col="red")
cutoff = -log10(quantile(all_mapped$P_res,0.1))
abline(v=cutoff,col="blue")
text(cutoff, 10, "Proposed cuttoff: top decile", 
     srt=90, adj=c(0,0), cex=1, col="blue")



res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

for (f in res_files) {
  res_list <- readRDS(f)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df <- do.call(rbind, res_list)
  summary(res_df)
  #plot(-log10(res_df$P_res),-log10(res_df$P_emb))
  #text(-log10(res_df$P_res),-log10(res_df$P_emb), res_df$Aligned_Position)
  plot(res_df$Aligned_Position, -log10(res_df$P_res), col="white", main=gene_name,
       ylim=c(0,1+max(c(-log10(res_df$P_res), -log10(res_df$P_emb)))),
       ylab="-log10(p)", xlab="Alignment Index")
  lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
  lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")
  bf <- 0.05 / (2 * max(res_df$Aligned_Position))
  abline(h=-log10(bf), col="black", lty=2)
  
  points(res_df[-log10(res_df$P_res) > -log10(bf),]$Aligned_Position,
         -log10(res_df$P_res[-log10(res_df$P_res) > -log10(bf)]), col="orange",pch=16)

  
  legend(
    "topleft",
    legend = c("Residue effect (P_res: residues | PCs)",
                   "Embedding effect (P_emb: embeddings | residues + PCs)",
                   "Bonferroni threshold (0.05 / (2 × seqLen))"),
    col = c("orange", "tomato", "black"),
    lty = c(1, 1, 2),
    lwd = c(2, 2, 1),
    bty = "n",  # remove box
    cex = 0.8
  )
}

plot(res_df$R2_partial)


library(ggplot2)

all_res <- lapply(res_files, function(f) {
  res_list <- readRDS(f)
  res_df <- do.call(rbind, res_list)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df$gene <- gene_name
  res_df
}) |> do.call(rbind, args = _)

# Melt into long format for ggplot
all_long <- tidyr::pivot_longer(all_res,
                                cols = c(P_res, P_emb),
                                names_to = "type",
                                values_to = "pval")

library(tidyverse)

res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

# Load and combine all results
all_res <- lapply(res_files, function(f) {
  res_list <- readRDS(f)
  res_df <- do.call(rbind, res_list)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df$gene <- gene_name
  res_df
}) |> bind_rows()

# Sort genes by name (optional)
all_res <- all_res %>% arrange(gene, Aligned_Position)

# Assign cumulative positions along x-axis
gene_lengths <- all_res %>%
  group_by(gene) %>%
  summarise(len = max(Aligned_Position, na.rm = TRUE))

gene_offsets <- c(0, cumsum(head(gene_lengths$len, -1)))
names(gene_offsets) <- gene_lengths$gene

all_res <- all_res %>%
  mutate(pos_global = Aligned_Position + gene_offsets[gene])

# Bonferroni threshold (using total sites)
bf <- 0.05 / (1 * sum(gene_lengths$len))
ylim_max <- 1 + max(-log10(all_res$P_res), na.rm = TRUE)

# Plot setup
plot(all_res$pos_global, -log10(all_res$P_res),
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "-log10(p)",
     ylim = c(0, ylim_max),
     main = "Chloroplast Protein Association Study",
     xaxt = "n")

# Alternate colors for genes
gene_colors <- setNames(rep(c("steelblue3", "grey60"),
                            length.out = length(unique(all_res$gene))),
                        unique(all_res$gene))

# Draw lines per gene
for (g in unique(all_res$gene)) {
  sub <- subset(all_res, gene == g)
  points(sub$pos_global, -log10(sub$P_res), col=gene_colors[g], lwd=1.5)
  #lines(sub$pos_global, -log10(sub$P_emb), col=adjustcolor(gene_colors[g], alpha.f=0.5), lty=2)
  points(sub[-log10(sub$P_res) > -log10(bf),]$pos_global,
         -log10(sub$P_res[-log10(sub$P_res) > -log10(bf)]), col=gene_colors[g],pch=16)
}

# Add Bonferroni threshold line
abline(h = -log10(bf), col = "black", lty = 2)

# Add gene labels at midpoints
midpoints <- gene_lengths %>%
  mutate(offset = gene_offsets,
         mid = offset + len / 2)
axis(1, at = midpoints$mid, labels = midpoints$gene, las = 2, cex.axis = 0.7)




library(jsonlite)

# --- Load UniProt JSON ---
j <- fromJSON("data/psbA_TM_uniprot.json")

# Extract start, end for Transmembrane regions
tm_uniprot <- j$features %>%
  filter(type == "Transmembrane") %>%
  transmute(
    start = location$start$value,
    end   = location$end$value
  )
# --- Map to your alignment ---
at_ID <- "AP000423.1"
at_df <- df_joined %>%
  filter(ID == at_ID) %>%
  select(Residue_Index, Aligned_Position)

at_start <- dplyr::rename(at_df, start_aln = Aligned_Position)
at_end   <- dplyr::rename(at_df, end_aln = Aligned_Position)

tm_df <- tm_uniprot %>%
  left_join(at_start, by = c("start" = "Residue_Index")) %>%
  left_join(at_end,   by = c("end"   = "Residue_Index")) %>%
  filter(!is.na(start_aln), !is.na(end_aln)) %>%
  mutate(Label = "TM helix")

res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

res_list <- readRDS("results/tmp_results_nested_psbA.rds")
gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
res_df <- do.call(rbind, res_list)







# --- Plot with TM regions ---
ggplot(a, aes(Aligned_Position, CVMRatio_smooth)) +
  geom_line(color="#2C7BB6", linewidth=0.8) +
  geom_rect(
    data=tm_df,
    aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
    inherit.aes=FALSE,
    fill="orange", alpha=0.15
  ) +
  geom_text(
    data=tm_df,
    aes(x=(start+end)/2, y=max(a$CVMRatio_smooth, na.rm=TRUE)*0.95, label=Label),
    color="orange", angle=90, size=3, vjust=0
  ) +
  theme_minimal(base_size=12) +
  labs(x="Aligned Position", y=expression(R^2)) +
  scale_x_continuous(limits=c(0,400))


# ----  rbcL 

# --- Map to your alignment ---
os_ID <- data$ID[grep("Oryza sativa", data$Organism)]
os_df <- df_joined %>%
  filter(ID == os_ID) %>%
  select(Residue_Index, Aligned_Position)

os_start <- dplyr::rename(os_df, start_aln = Aligned_Position)
os_end   <- dplyr::rename(os_df, end_aln = Aligned_Position)


res_list <- readRDS("results/tmp_results_nested_rbcL.rds")
res_df <- do.call(rbind, res_list)

feat_json <- 'data/rbcL_features.json'
feats <- as.data.frame(fromJSON(feat_json) ) 
feat_df <- feats %>%
  transmute(Residue_Index = features.location$start$value,
            Label = features.ligand$name,
            Desc = features.description) %>%
  left_join(os_df, by="Residue_Index") %>%
  filter(!is.na(Aligned_Position))

plot(res_df$Aligned_Position, -log10(res_df$P_res), col="orange",
     main="Rubisco large subunit (rbcL) temp associations",
     ylim=c(0,1+max(c(-log10(res_df$P_res)))),
     ylab="-log10(p)", xlab="Alignment Index")
##lines(res_df$Aligned_Position, -log10(res_df$P_res), col="orange")
#lines(res_df$Aligned_Position, -log10(res_df$P_emb), col="tomato")
bf <- 0.05 / (1 * max(res_df$Aligned_Position))
abline(h=-log10(bf), col="black", lty=2)

points(res_df[-log10(res_df$P_res) > -log10(bf),]$Aligned_Position,
       -log10(res_df$P_res[-log10(res_df$P_res) > -log10(bf)]), col="orange",pch=16)
#text(305, 6, "DE-loop", col="tomato")

for (i in 1:nrow(feat_df)) {
  abline(v = feat_df[i,]$Aligned_Position, col="blue", lty=2, lwd=0.6)
  #text(feat_df[i,]$Aligned_Position, 1.01, feat_df[i,]$Label, 
  #     srt=90, adj=c(0,0), cex=0.8, col="black")
  #text(feat_df[i,]$Aligned_Position-1, 15, feat_df[i,]$Desc, 
  #     srt=90, adj=c(0,0), cex=2, col="black")
}

legend(
  "topleft",
  legend = c("Residue p value (residues | PCs)",
             "Active, binding, or functional site",
             "Bonferroni threshold"),
  col = c("orange", "blue", "black"),
  lty = c(1, 2, 2),
  lwd = c(2, 2, 1),
  cex = 0.8
)

#big guys are at aln 439 and 420








