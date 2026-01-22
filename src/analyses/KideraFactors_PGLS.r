library(ape)
library(data.table)
library(caper)
library(Peptides)
library(pheatmap)

# ==== CONFIGURATION ====
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
gwas_dir <- "results/residue_models_triple/"
partition_map_file <- "raxml_input/partitionMap.rds"

stopifnot(file.exists(aln_file), file.exists(data_file), file.exists(tree_file))

# ==== LOAD GWAS RESULTS (from previous script) ====
cat("Loading GWAS results...\n")
partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_effects <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    data.table(
      Gene = m$Gene,
      GenePos = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs,
      P_aa_only = m$P_aa_only
    )
  }), fill = TRUE)
}), fill = TRUE)

gwas_effects <- unique(gwas_effects)
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

# Define thresholds
thresh_control <- quantile(gwas_effects$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(gwas_effects$P_aa_only, 0.20, na.rm = TRUE)

gwas_effects[, sig_class := fcase(
  P_aa_with_pcs < thresh_control & P_aa_only < thresh_nocontrol, "sig_both",
  P_aa_with_pcs < thresh_control & P_aa_only >= thresh_nocontrol, "sig_control",
  P_aa_with_pcs >= thresh_control & P_aa_only < thresh_nocontrol, "sig_nocontrol",
  default = "not_sig"
)]

cat("\nSite classification:\n")
print(gwas_effects[, .N, by = sig_class][order(-N)])

# ==== LOAD DATA ====
cat("Loading alignment...\n")
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)

cat("Loading phenotype data...\n")
data <- as.data.table(arrow::read_parquet(data_file))
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID

cat("Loading tree...\n")
tree <- read.tree(tree_file)
tree$node.label <- NULL

# ==== FUNCTION: CALCULATE KIDERA FACTORS FOR SITES ====
calc_kidera <- function(aln_mat, positions) {
  stopifnot(all(positions <= ncol(aln_mat)))
  
  sub_mat <- aln_mat[, positions, drop = FALSE]
  
  # Collapse each row to a sequence string, replacing gaps/ambiguous with ""
  seqs <- apply(sub_mat, 1, function(row) {
    row[!row %in% c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- ""
    paste(row, collapse = "")
  })
  
  # Calculate Kidera factors for each sequence
  kf_list <- lapply(seqs, function(s) {
    if (nchar(s) < 1) return(rep(NA, 10))
    unlist(kideraFactors(s))
  })
  
  kf_mat <- do.call(rbind, kf_list)
  rownames(kf_mat) <- rownames(aln_mat)
  colnames(kf_mat) <- paste0("KF", 1:10)
  
  as.data.table(kf_mat, keep.rownames = "taxon")
}

# ==== CALCULATE KIDERA BY SIG CLASS ====
cat("\nCalculating Kidera factors by sig_class...\n")

sig_classes <- c("sig_both", "sig_control", "sig_nocontrol", "not_sig")
kidera_list <- list()

for (sc in sig_classes) {
  pos <- gwas_effects[sig_class == sc, GlobalPos]
  pos <- pos[pos <= ncol(aln_mat)]
  
  if (length(pos) < 10) {
    cat("Skipping", sc, "- too few positions\n")
    next
  }
  
  cat(sc, ":", length(pos), "positions\n")
  kf <- calc_kidera(aln_mat, pos)
  kf[, sig_class := sc]
  kidera_list[[sc]] <- kf
}

# ==== PREPARE PGLS DATA ====
cat("\nPreparing PGLS data...\n")

common_taxa <- Reduce(intersect, list(
  rownames(aln_mat),
  names(pheno)[!is.na(pheno)],
  tree$tip.label
))
cat("Common taxa:", length(common_taxa), "\n")
stopifnot(length(common_taxa) > 100)

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))
pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree_pruned$tip.label)
tree_pruned <- root(tree_pruned, outgroup = pinales_in_tree, resolve.root = TRUE)

# ==== RUN PGLS BY SIG CLASS ====
cat("\n=== PGLS Results (Kidera Factors) ===\n")

pgls_results_kidera <- list()

for (sc in names(kidera_list)) {
  cat("\n---", sc, "---\n")
  
  kf <- kidera_list[[sc]]
  kf <- kf[taxon %in% common_taxa]
  kf[, temp := pheno[taxon]]
  kf <- kf[!is.na(temp)]
  
  # Check for valid Kidera values
  kf <- kf[complete.cases(kf[, paste0("KF", 1:10), with = FALSE])]
  
  if (nrow(kf) < 100) {
    cat("Too few complete cases, skipping\n")
    next
  }
  
  form <- as.formula(paste("temp ~", paste(paste0("KF", 1:10), collapse = " + ")))
  
  df <- as.data.frame(kf)
  rownames(df) <- df$taxon
  
  comp_data <- comparative.data(
    phy = tree_pruned,
    data = df,
    names.col = taxon,
    vcv = TRUE,
    warn.dropped = FALSE
  )
  
  pgls_fit <- tryCatch({
    pgls(form, data = comp_data, lambda = "ML")
  }, error = function(e) {
    cat("PGLS failed:", e$message, "\n")
    NULL
  })
  
  if (is.null(pgls_fit)) next
  
  summ <- summary(pgls_fit)
  cat("Lambda:", pgls_fit$param["lambda"], "\n")
  cat("R-squared:", summ$r.squared, "\n")
  
  coefs <- as.data.table(summ$coefficients, keep.rownames = "term")
  setnames(coefs, c("term", "estimate", "se", "t", "p"))
  
  sig_coefs <- coefs[p < 0.05 & term != "(Intercept)"]
  if (nrow(sig_coefs) > 0) {
    cat("\nSignificant Kidera factors (p < 0.05):\n")
    print(sig_coefs[order(p)])
  }
  
  pgls_results_kidera[[sc]] <- list(
    model = pgls_fit,
    summary = summ,
    coefs = coefs,
    lambda = pgls_fit$param["lambda"],
    r2 = summ$r.squared,
    n_sites = gwas_effects[sig_class == sc, .N]
  )
}

saveRDS(pgls_results_kidera, "data/tmp/pgls_kidera_results.rds")

pgls_results_kidera <- readRDS("data/tmp/pgls_kidera_results.rds")
# ==== HEATMAP ====
coef_dt <- rbindlist(lapply(names(pgls_results_kidera), function(sc) {
  res <- pgls_results_kidera[[sc]]
  coefs <- res$coefs
  coefs[, sig_class := sc]
  coefs
}))
coef_dt <- coef_dt[term != "(Intercept)"]

summary(pgls_results_kidera[[1]])
summary(pgls_results_kidera$not_sig$model)
summary(pgls_results_kidera$sig_nocontrol$model)
summary(pgls_results_kidera$sig_control$model)
summary(pgls_results_kidera$sig_both$model)

pgls_results_kidera$sig_both

coef_mat <- dcast(coef_dt, term ~ sig_class, value.var = "estimate", fill = 0)
rn <- coef_mat$term
coef_mat[, term := NULL]
coef_mat <- as.matrix(coef_mat)
rownames(coef_mat) <- rn
coef_mat
p_mat <- dcast(coef_dt, term ~ sig_class, value.var = "p", fill = 1)
p_mat[, term := NULL]
p_mat <- as.matrix(p_mat)
rownames(p_mat) <- rn

sig_stars <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_stars[p_mat < 0.05] <- "*"
sig_stars[p_mat < 0.01] <- "**"
sig_stars[p_mat < 0.001] <- "***"

max_abs <- max(abs(coef_mat))
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(coef_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = sig_stars,
         fontsize_number = 12,
         breaks = breaks,
         color = colors,
         main = "PGLS Coefficients for Kidera Factors\n(* p<0.05, ** p<0.01, *** p<0.001)")


# ---- plot regressions ---- 
library(ggplot2)
library(patchwork)

# Reload necessary data
pgls_results_kidera <- readRDS("data/tmp/pgls_kidera_results.rds")
data <- as.data.table(arrow::read_parquet(data_file))
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID

# Reconstruct kidera_list from saved results or recalculate
# (We need the raw Kidera values for plotting)
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)

partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
gwas_effects <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    data.table(Gene = m$Gene, GenePos = m$Aligned_Position,
               P_aa_with_pcs = m$P_aa_with_pcs, P_aa_only = m$P_aa_only)
  }), fill = TRUE)
}), fill = TRUE)
gwas_effects <- unique(gwas_effects)
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

thresh_control <- quantile(gwas_effects$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(gwas_effects$P_aa_only, 0.20, na.rm = TRUE)
gwas_effects[, sig_class := fcase(
  P_aa_with_pcs < thresh_control & P_aa_only < thresh_nocontrol, "sig_both",
  P_aa_with_pcs < thresh_control & P_aa_only >= thresh_nocontrol, "sig_control",
  P_aa_with_pcs >= thresh_control & P_aa_only < thresh_nocontrol, "sig_nocontrol",
  default = "not_sig"
)]

# Recalculate Kidera factors
calc_kidera <- function(aln_mat, positions) {
  stopifnot(all(positions <= ncol(aln_mat)))
  sub_mat <- aln_mat[, positions, drop = FALSE]
  seqs <- apply(sub_mat, 1, function(row) {
    row[!row %in% c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- ""
    paste(row, collapse = "")
  })
  kf_list <- lapply(seqs, function(s) {
    if (nchar(s) < 1) return(rep(NA, 10))
    unlist(kideraFactors(s))
  })
  kf_mat <- do.call(rbind, kf_list)
  rownames(kf_mat) <- rownames(aln_mat)
  colnames(kf_mat) <- paste0("KF", 1:10)
  as.data.table(kf_mat, keep.rownames = "taxon")
}

sig_classes <- c("sig_both", "sig_control", "sig_nocontrol", "not_sig")
kidera_list <- list()
for (sc in sig_classes) {
  pos <- gwas_effects[sig_class == sc, GlobalPos]
  pos <- pos[pos <= ncol(aln_mat)]
  if (length(pos) >= 10) {
    kidera_list[[sc]] <- calc_kidera(aln_mat, pos)
  }
}

# Build combined data for plotting
plot_data <- rbindlist(lapply(names(kidera_list), function(sc) {
  kf <- copy(kidera_list[[sc]])
  kf[, phenotype := pheno[taxon]]
  kf[, sig_class := sc]
  kf[complete.cases(kf)]
}))

# Melt to long format for faceting
kf_cols <- paste0("KF", 1:10)
plot_long <- melt(plot_data, id.vars = c("taxon", "phenotype", "sig_class"),
                  measure.vars = kf_cols, variable.name = "KF", value.name = "value")

# Extract p-values for annotation
coef_dt <- rbindlist(lapply(names(pgls_results_kidera), function(sc) {
  res <- pgls_results_kidera[[sc]]
  coefs <- res$coefs
  coefs[, sig_class := sc]
  coefs
}))
coef_dt <- coef_dt[term != "(Intercept)"]
setnames(coef_dt, "term", "KF")

plot_long <- merge(plot_long, coef_dt[, .(KF, sig_class, p, estimate)],
                   by = c("KF", "sig_class"), all.x = TRUE)

# Create faceted plot
p <- ggplot(plot_long, aes(x = value, y = phenotype)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
  facet_grid(sig_class ~ KF, scales = "free_x") +
  geom_text(data = unique(plot_long[, .(KF, sig_class, p, estimate)]),
            aes(x = -Inf, y = Inf, 
                label = sprintf("β=%.2f\np=%.3f", estimate, p)),
            hjust = -0.1, vjust = 1.2, size = 2.5) +
  labs(x = "Kidera Factor Value", y = "Phenotype (bio_8_p50)",
       title = "Kidera Factors vs Phenotype by Significance Class") +
  theme_bw() +
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 6))
p
ggsave("results/kidera_regression_plots.pdf", p, width = 16, height = 10)
cat("Saved: results/kidera_regression_plots.pdf\n")

kf_labels <- c(
  KF1 = "KF1: Helix/bend preference",
  KF2 = "KF2: Side-chain size",
  KF3 = "KF3: Extended structure preference",
  KF4 = "KF4: Hydrophobicity",
  KF5 = "KF5: Double-bend preference",
  KF6 = "KF6: Partial specific volume",
  KF7 = "KF7: Flat extended preference",
  KF8 = "KF8: Occurrence in alpha region",
  KF9 = "KF9: pK-C",
  KF10 = "KF10: Surrounding hydrophobicity"
)

plot_long[, KF := factor(KF, levels = names(kf_labels), labels = kf_labels)]

# Update annotation data too
annot_data <- unique(plot_long[, .(KF, sig_class, p, estimate)])

# Create faceted plot
p <- ggplot(plot_long, aes(x = value, y = phenotype)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 0.8) +
  facet_grid(sig_class ~ KF, scales = "free_x") +
  geom_text(data = annot_data,
            aes(x = -Inf, y = Inf, 
                label = sprintf("β=%.2f\np=%.3f", estimate, p)),
            hjust = -0.1, vjust = 1.2, size = 2.5) +
  labs(x = "Kidera Factor Value", y = "Phenotype (bio_8_p50)",
       title = "Kidera Factors vs Phenotype by Significance Class") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 7, hjust = 0),
        strip.text.y = element_text(size = 8),
        axis.text = element_text(size = 6))
p
ggsave("results/kidera_regression_plots.pdf", p, width = 18, height = 10)
# Extract model R² and p-values for sig_class labels
model_stats <- rbindlist(lapply(names(pgls_results_kidera), function(sc) {
  res <- pgls_results_kidera[[sc]]
  # Get overall model p-value from F-statistic
  f_stat <- res$summary$fstatistic
  model_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  data.table(
    sig_class = sc,
    r2 = res$r2,
    model_p = model_p
  )
}))

# Create new labels
model_stats[, sig_label := sprintf("%s (R²=%.2f, p=%.2f)", sig_class, r2, model_p)]

# Map to plot_long
plot_long <- merge(plot_long, model_stats[, .(sig_class, sig_label)], by = "sig_class", all.x = TRUE)
plot_long[, sig_class := factor(sig_label, levels = model_stats[order(sig_class), sig_label])]

# Update annotation data
annot_data <- unique(plot_long[, .(KF, sig_class, p, estimate)])

# Create faceted plot
p <- ggplot(plot_long, aes(x = value, y = phenotype)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 0.8) +
  facet_grid(sig_class ~ KF, scales = "free_x") +
  geom_text(data = annot_data,
            aes(x = -Inf, y = Inf, 
                label = sprintf("β=%.2f\np=%.3f", estimate, p)),
            hjust = -0.1, vjust = 1.2, size = 2.5) +
  labs(x = "Kidera Factor Value", y = "Phenotype (bio_8_p50)",
       title = "Kidera Factors vs Phenotype by Significance Class") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 7, hjust = 0),
        strip.text.y = element_text(size = 8),
        axis.text = element_text(size = 6))


p


# ---- kf5 by gene ----

# Get positions for each sig_class
calc_kf5_by_gene <- function(aln_mat, positions, partition_map) {
  stopifnot(all(positions <= ncol(aln_mat)))
  
  # Map positions back to genes
  pos_genes <- partition_map[GlobalPos %in% positions, .(Gene, GlobalPos)]
  
  # Calculate KF5 per gene
  gene_kf5 <- rbindlist(lapply(unique(pos_genes$Gene), function(g) {
    gene_pos <- pos_genes[Gene == g, GlobalPos]
    if (length(gene_pos) < 1) return(NULL)
    
    sub_mat <- aln_mat[, gene_pos, drop = FALSE]
    seqs <- apply(sub_mat, 1, function(row) {
      row[!row %in% c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- ""
      paste(row, collapse = "")
    })
    
    kf5_vals <- sapply(seqs, function(s) {
      if (nchar(s) < 1) return(NA)
      kideraFactors(s)[[1]][5]  # KF5 only
    })
    
    data.table(taxon = names(kf5_vals), Gene = g, KF5 = kf5_vals, n_sites = length(gene_pos))
  }))
  
  gene_kf5
}

# Analyze KF5 by gene for each sig_class
kf5_gene_stats <- rbindlist(lapply(names(kidera_list), function(sc) {
  pos <- gwas_effects[sig_class == sc, GlobalPos]
  pos <- pos[pos <= ncol(aln_mat)]
  if (length(pos) < 10) return(NULL)
  
  gene_kf5 <- calc_kf5_by_gene(aln_mat, pos, partition_map)
  gene_kf5[, sig_class := sc]
  
  # Summary stats per gene
  gene_kf5[!is.na(KF5), .(
    mean_KF5 = mean(KF5),
    sd_KF5 = sd(KF5),
    range_KF5 = max(KF5) - min(KF5),
    n_sites = unique(n_sites),
    n_taxa = .N
  ), by = .(Gene, sig_class)]
}))

# Which genes contribute most to KF5 variance?
cat("\n=== Genes contributing to KF5 variance ===\n")
hist(kf5_gene_stats$mean_KF5)
hist(kf5_gene_stats$sd_KF5)
plot(kf5_gene_stats$mean_KF5,kf5_gene_stats$sd_KF5)
boxplot(mean_KF5 ~ Gene, kf5_gene_stats[sig_class=="not_sig"])
boxplot(mean_KF5 ~ Gene, kf5_gene_stats[sig_class=="sig_both"])
boxplot(mean_KF5 ~ Gene, kf5_gene_stats[sig_class=="sig_control"])
boxplot(mean_KF5 ~ Gene, kf5_gene_stats[sig_class=="sig_nocontrol"])
print(kf5_gene_stats[order(-sd_KF5)][, .SD[1:5], by = sig_class])

# Plot KF5 distribution by gene for sig_nocontrol (your interesting class)
kf5_by_gene <- rbindlist(lapply(c("not_sig", "sig_both"), function(sc) {
  pos <- gwas_effects[sig_class == sc, GlobalPos]
  pos <- pos[pos <= ncol(aln_mat)]
  gene_kf5 <- calc_kf5_by_gene(aln_mat, pos, partition_map)
  gene_kf5[, sig_class := sc]
  gene_kf5
}))

head(kf5_by_gene)
summary(kf5_by_gene)
kf5_by_gene_clean <- kf5_by_gene

kf5_by_gene_clean$taxon <- gsub("\\.KF5$", "", kf5_by_gene_clean$taxon)

kf5_by_gene_clean[, phenotype := pheno[taxon]]
hist(kf5_by_gene_clean$KF5)
hist(kf5_by_gene_clean$phenotype)
sum(is.na(kf5_by_gene_clean$phenotype))
#kf5_by_gene <- kf5_by_gene[!is.na(KF5) & !is.na(phenotype)]

# Top variable genes
top_genes <- kf5_gene_stats[sig_class == "sig_nocontrol"][order(-sd_KF5), head(Gene, 6)]

p_genes <- ggplot(kf5_by_gene_clean[Gene %in% top_genes], 
                  aes(x = KF5, y = phenotype)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  facet_grid(sig_class ~ Gene, scales = "free") +
  labs(x = "KF5 (Double-bend preference)", y = "Phenotype",
       title = "Gene-specific KF5 vs Phenotype") +
  theme_bw()

table(kf5_by_gene$Gene)

p_genes <- ggplot(kf5_by_gene_clean, 
                  aes(x = KF5, y = phenotype)) +
  geom_point(alpha = 0.3, size = 0.8) +
  #geom_smooth(method = "lm", se = TRUE, color = "red") +
  facet_grid(sig_class ~ Gene, scales = "free") +
  labs(x = "KF5 (Double-bend preference)", y = "Phenotype",
       title = "Gene-specific KF5 vs Phenotype") +
  theme_bw()

p_genes

ggsave("results/kf5_by_gene_allw.pdf", p_genes, width = 50, height = 10, limitsize=F)

p_genes <- ggplot(kf5_by_gene_clean, 
                  aes(x = KF5, y = phenotype, color = sig_class)) +
  geom_point(alpha = 0.3, size = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Gene, scales = "free_x") +
  labs(x = "KF5 (Double-bend preference)", y = "Phenotype",
       title = "Gene-specific KF5 vs Phenotype", color = "Sig Class") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("results/kf5_by_gene_all_x.pdf", p_genes, width = 12, height = 10)


p_genes

ggsave("results/kf5_by_gene.pdf", p_genes, width = 14, height = 6)
cat("Saved: results/kf5_by_gene.pdf\n")

# ----kf5 stats----
# ==== STATISTICAL TESTS FOR KF5 BY GENE ====

# 1. Per-gene correlations between KF5 and phenotype
cat("\n=== KF5-Phenotype correlations by Gene and Sig Class ===\n")
kf5_cors <- kf5_by_gene_clean[, .(
  r = cor(KF5, phenotype, use = "complete.obs"),
  p = cor.test(KF5, phenotype)$p.value,
  n = .N
), by = .(Gene, sig_class)]

# Significant correlations
cat("\nSignificant correlations (p < 0.05):\n")
print(kf5_cors[p < 0.05][order(p)])

# 2. Which genes show opposite effects in different sig classes?
cat("\n=== Genes with divergent KF5 effects across sig classes ===\n")
kf5_wide <- dcast(kf5_cors, Gene ~ sig_class, value.var = "r")
kf5_wide[, range_r := do.call(pmax, c(.SD, na.rm = TRUE)) - do.call(pmin, c(.SD, na.rm = TRUE)), 
         .SDcols = patterns("sig_")]
print(kf5_wide[order(-range_r)][1:10])

# 3. Mean KF5 differences across sig classes (is selection shifting loop propensity?)
cat("\n=== Mean KF5 by gene and sig class ===\n")
kf5_means <- kf5_by_gene_clean[, .(
  mean_KF5 = mean(KF5, na.rm = TRUE),
  sd_KF5 = sd(KF5, na.rm = TRUE)
), by = .(Gene, sig_class)]

kf5_means_wide <- dcast(kf5_means, Gene ~ sig_class, value.var = "mean_KF5")
print(kf5_means_wide)

# 4. ANOVA: does KF5 differ by sig_class within genes?
cat("\n=== ANOVA: KF5 variation by sig_class within genes ===\n")
anova_results <- kf5_by_gene_clean[, {
  if (length(unique(sig_class)) > 1 && .N > 30) {
    fit <- aov(KF5 ~ sig_class, data = .SD)
    s <- summary(fit)[[1]]
    .(F_val = s$`F value`[1], p = s$`Pr(>F)`[1])
  } else {
    .(F_val = NA_real_, p = NA_real_)
  }
}, by = Gene]

cat("\nGenes where KF5 differs significantly by sig_class:\n")
print(anova_results[p < 0.05][order(p)])

# 5. Biological interpretation: which genes are these?
cat("\n=== Gene functions (plastid context) ===\n")
gene_functions <- data.table(
  Gene = c("rbcL", "matK", "psbA", "psbB", "psbC", "psbD", "atpA", "atpB", 
           "atpE", "atpF", "atpH", "atpI", "petA", "petB", "petD", "ndhA",
           "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI",
           "ndhJ", "ndhK", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps2", "rps3",
           "rps4", "rps7", "rps8", "rps11", "rps12", "rps14", "rps15", "rps16",
           "rps18", "rps19", "rpl2", "rpl14", "rpl16", "rpl20", "rpl22", "rpl23",
           "rpl32", "rpl33", "rpl36", "ycf1", "ycf2", "ycf3", "ycf4", "clpP",
           "accD", "ccsA", "cemA", "infA"),
  Function = c("RuBisCO large subunit", "Maturase", "PSII reaction center D1", 
               "PSII CP47", "PSII CP43", "PSII D2", "ATP synthase CF1 alpha",
               "ATP synthase CF1 beta", "ATP synthase CF1 epsilon", "ATP synthase CF0 B",
               "ATP synthase CF0 C", "ATP synthase CF0 A", "Cytochrome b6f complex",
               "Cytochrome b6", "Cytochrome b6f subunit IV", "NADH dehydrogenase",
               "NADH dehydrogenase", "NADH dehydrogenase", "NADH dehydrogenase",
               "NADH dehydrogenase", "NADH dehydrogenase", "NADH dehydrogenase",
               "NADH dehydrogenase", "NADH dehydrogenase", "NADH dehydrogenase",
               "NADH dehydrogenase", "RNA polymerase alpha", "RNA polymerase beta",
               "RNA polymerase beta'", "RNA polymerase beta''", rep("Ribosomal protein S", 12),
               rep("Ribosomal protein L", 11), "Ycf1 (translocon)", "Ycf2 (ATPase)",
               "Ycf3 (PSI assembly)", "Ycf4 (PSI assembly)", "Clp protease",
               "Acetyl-CoA carboxylase", "Cytochrome c biogenesis", "Envelope membrane protein",
               "Translation initiation")
)

# Match to our significant genes
sig_genes <- kf5_cors[p < 0.05, unique(Gene)]
cat("\nFunctions of genes with significant KF5-phenotype correlations:\n")
print(gene_functions[Gene %in% sig_genes])

# 6. Summary: does KF5 (loop/turn propensity) relate to temperature adaptation?
cat("\n=== Biological interpretation summary ===\n")
cat("Phenotype: Mean temperature of wettest quarter (bio_8)\n")
cat("KF5: Double-bend preference (β-turns, loops)\n\n")

# Overall correlation direction
overall_cor <- kf5_by_gene_clean[, cor(KF5, phenotype, use = "complete.obs"), by = sig_class]
print(overall_cor)
cat("\nPositive correlation = higher KF5 (more loops/turns) in warmer climates\n")
cat("Negative correlation = lower KF5 (fewer loops/turns) in warmer climates\n")
