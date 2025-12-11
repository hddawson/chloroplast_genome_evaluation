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

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes_to_process <- get_gene(aln_files)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

# ---------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------

set.seed(123)

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  
  out_path <- paste0("results/residue_models/", gene, "_effects.rds")
  
  if (file.exists(out_path)) {
    message("Skipping ", gene, " (output already exists)")
    next
  }
  
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # --- Read alignment ---
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # --- Strip IDs ---
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  seq_ids <- names(aln)
  
  # Get common IDs
  common_ids <- intersect(seq_ids, pheno_pcs$ID)
  if (length(common_ids) < 8500L) {
    message("Skipping ", gene, ": insufficient overlap with phenotype data")
    next
  }
  
  # Subset to common IDs
  aln_mat <- aln_mat[match(common_ids, seq_ids), , drop = FALSE]
  pheno_pcs_sub <- pheno_pcs[ID %in% common_ids]
  setkey(pheno_pcs_sub, ID)
  pheno_pcs_sub <- pheno_pcs_sub[match(common_ids, ID)]
  
  y <- pheno_pcs_sub$pheno
  X_pcs <- scale(as.matrix(pheno_pcs_sub[, ..pc_names]))
  
  # -------------------------------------------------------------------
  # LOOP OVER POSITIONS
  # -------------------------------------------------------------------
  
  positions <- seq_len(ncol(aln_mat))
  results_list <- vector("list", length(positions))
  
  for (pos in positions) {
    
    residues <- aln_mat[, pos]
    residue_table <- table(residues)
    
    # Skip gaps-only or monomorphic positions
    if (length(residue_table) < 2L || all(names(residue_table) == "-")) next
    
    # Skip positions with only gaps
    non_gap <- residues != "-"
    if (sum(non_gap) < 8500L) next
    
    # Check special case
    if (pos == 19 & gene == "psbJ") next
    
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    var0 <- which(apply(X_aa, 2, var) == 0)
    if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # -----------------------------------------------------------------
    # FIT MODELS
    # -----------------------------------------------------------------
    
    fit_reduced <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    p_res <- tryCatch(anova(fit_reduced, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    # Extract only residue effect sizes and SEs
    coef_summary <- summary(fit_full)$coefficients
    aa_coefs <- grep("^X_aa", rownames(coef_summary), value = TRUE)
    
    if (length(aa_coefs) > 0) {
      effects_dt <- data.table(
        Residue = sub("res_factor", "", aa_coefs),
        Effect = coef_summary[aa_coefs, "Estimate"],
        SE = coef_summary[aa_coefs, "Std. Error"],
        T_value = coef_summary[aa_coefs, "t value"],
        P_value = coef_summary[aa_coefs, "Pr(>|t|)"]
      )
    } else {
      effects_dt <- NULL
    }
    
    results_list[[pos]] <- list(
      Gene = gene,
      Aligned_Position = pos,
      N = sum(non_gap),
      R2_reduced = r2_reduced,
      R2_full = r2_full,
      P_res = p_res,
      residue_counts = as.list(residue_table),
      effects = effects_dt,
      IDs = common_ids
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) > 0) {
    dir.create("results/residue_models", showWarnings = FALSE, recursive = TRUE)
    saveRDS(results_list, out_path)
    message("\nCompleted ", gene, ": ", length(results_list), " positions with effects")
  } else {
    message("No positions retained for gene ", gene)
  }
}

message("\nDone processing all genes")

library(data.table)
library(ggplot2)

# Load all model files
model_files <- list.files("results/residue_models/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

# Extract summary stats from each gene
all_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      N = m$N,
      R2_reduced = m$R2_reduced,
      R2_full = m$R2_full,
      Delta_R2 = m$R2_full - m$R2_reduced,
      P = m$P_res
    )
  }))
}))

stopifnot(nrow(all_results) > 0)

# Calculate -log10(P)
all_results[, neg_log10_P := -log10(P)]
all_results[is.infinite(neg_log10_P), neg_log10_P := NA]

# Create output directory
dir.create("results/manhattan_plots", showWarnings = FALSE, recursive = TRUE)

# Plot for each gene
genes <- unique(all_results$Gene)

for (gene in genes) {
  gene_data <- all_results[Gene == gene]
  
  p <- ggplot(gene_data, aes(x = Position, y = neg_log10_P)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
      title = paste0("Manhattan Plot: ", gene),
      x = "Aligned Position",
      y = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  ggsave(
    filename = paste0("results/manhattan_plots/", gene, "_manhattan.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("Saved Manhattan plot for ", gene)
}

message("\nCompleted all Manhattan plots")

all_results[, norm_pos := Position / max(Position), by = Gene]
all_results[, neg_log10_P := -log10(P)]

# optional cap for extreme values
all_results[is.infinite(neg_log10_P), neg_log10_P := NA]

ggplot(all_results[!is.na(neg_log10_P)],
       aes(x = norm_pos, weight = neg_log10_P)) +
  geom_density(adjust = 1.2, linewidth = 1) +
  labs(
    x = "Normalized position (0–1)",
    y = "Weighted density",
    title = "Smoothed significance density along normalized gene length"
  ) +
  theme_bw()

x <- all_results$norm_pos
w <- all_results$neg_log10_P
d <- density(x, weights = w / sum(w, na.rm=TRUE), na.rm=TRUE)

plot(d, type="l", lwd=2,
     xlab="Normalized position (0–1)",
     ylab="Weighted density")
library(stats)

ord <- order(all_results$norm_pos)
x <- all_results$norm_pos[ord]
y <- all_results$neg_log10_P[ord]

k <- 200
sm <- filter(y, rep(1/k, k), sides=2)

plot(x, sm, type="l", lwd=2,
     xlab="Normalized position",
     ylab="Running-average significance")

