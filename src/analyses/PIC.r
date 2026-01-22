#!/usr/bin/env Rscript
# ==============================================================================
# PIC-Based Analysis of GWAS Site-Phenotype Associations
# Tests whether GWAS-identified residue effects predict phenotype using
# phylogenetically independent contrasts
# ==============================================================================

library(ape)
library(data.table)
library(phangorn)
library(arrow)

# ==== CONFIGURATION ====
# File paths - UPDATE THESE
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

# ==== LOAD AND PREPARE DATA ====

cat("Loading tree...\n")
tree <- read.tree(tree_file)
stopifnot(file.exists(tree_file))

cat("Loading phenotype data...\n")
data <- read_parquet(data_file)
setDT(data)

# Root tree on Pinales
pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree$tip.label)
tree <- root(tree, outgroup = pinales_in_tree, resolve.root = TRUE)
stopifnot(is.rooted(tree))

cat("Tree:", Ntip(tree), "tips\n")

# ==== LOAD TIP SEQUENCES ====

cat("Loading alignment...\n")
tip_aln <- read.FASTA(aln_file, type = "AA")
tip_seqs_list <- sapply(tip_aln, function(x) paste(rawToChar(x, multiple = TRUE), collapse = ""))
names(tip_seqs_list) <- names(tip_aln)

stopifnot(all(tree$tip.label %in% names(tip_seqs_list)))

# Pad sequences to same length if needed
max_len <- max(nchar(tip_seqs_list))
tip_seqs_list <- sapply(tip_seqs_list, function(s) {
  if (nchar(s) < max_len) {
    paste0(paste(rep("-", max_len - nchar(s)), collapse = ""), s)
  } else s
}, USE.NAMES = TRUE)


stopifnot(length(unique(nchar(tip_seqs_list))) == 1)
aln_length <- nchar(tip_seqs_list[1])
cat("Alignment length:", aln_length, "\n")

# ==== LOAD GWAS RESULTS ====

cat("Loading GWAS results...\n")
partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_effects <- rbindlist(lapply(model_files, function(f) {
  cat(f, "\n")
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    dt <- copy(m$effects)
    dt[, Gene := m$Gene]
    dt[, GenePos := m$Aligned_Position]
    dt[, P_aa_with_pcs := m$P_aa_with_pcs]
    dt[, P_aa_only := m$P_aa_only]
    dt[, residue := gsub("X_aares_factor|X_aa", "", Residue)]
    dt[, .(Gene, GenePos, residue, Effect, SE, P_residue = P_value, P_aa_with_pcs, P_aa_only)]
  }))
}))

# Add GlobalPos
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

cat("GWAS effects loaded:", nrow(gwas_effects), "residue-site combinations\n")
cat("Unique sites:", uniqueN(gwas_effects$GlobalPos), "\n")

hist(gwas_effects$Effect)
fit <- lm(Effect ~ residue + P_aa_with_pcs, gwas_effects)
summary(fit)

site_pvals <- unique(gwas_effects[, .(Gene, GenePos, P_aa_with_pcs, P_aa_only)])
thresh_control <- quantile(site_pvals$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(site_pvals$P_aa_only, 0.20, na.rm = TRUE)

gwas_effects[, sig_class := {
  sig_ctrl <- P_aa_with_pcs < thresh_control
  sig_noctrl <- P_aa_only < thresh_nocontrol
  
  fcase(
    sig_ctrl & sig_noctrl, "sig_both",
    sig_ctrl & !sig_noctrl, "sig_control",
    !sig_ctrl & sig_noctrl, "sig_nocontrol",
    default = "not_sig"
  )
}]

gwas_effects$variant_sig_class <- gwas_effects$sig_class
gwas_effects$variant_sig_class[gwas_effects$P_residue >= 0.05] <- "not_sig"

table(gwas_effects$variant_sig_class)
table(gwas_effects$sig_class)

fit <- lm(Effect ~ residue*variant_sig_class, gwas_effects)
summary(fit)
# ==== PREPARE PHENOTYPE ====

pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID
tip_pheno <- pheno[tree$tip.label]

# Check for missing phenotypes
n_missing <- sum(is.na(tip_pheno))
cat("Tips with missing phenotype:", n_missing, "of", length(tip_pheno), "\n")

# ==== FUNCTION: PIC ANALYSIS FOR ONE SITE ====

#' Run PIC-based analysis for a single alignment position
#' 
#' For a given site, this:
#' 1. Extracts the AA state at each tip
#' 2. Assigns numeric "effect score" based on GWAS residue effects
#' 3. Computes PICs for both effect scores and phenotype
#' 4. Tests correlation between PICs (should be 0 under null)
#'
#' @param global_pos Integer, position in supermatrix
#' @param tree Phylo object (rooted, ultrametric not required)
#' @param tip_seqs Named character vector of sequences
#' @param tip_pheno Named numeric vector of phenotypes
#' @param effects_dt data.table with columns: residue, Effect
#' @return List with correlation, p-value, n_tips, etc.

pic_site_analysis <- function(global_pos, tree, tip_seqs, tip_pheno, effects_dt) {
  
  # Extract AA at this position for all tips
  tip_aa <- sapply(tip_seqs[tree$tip.label], function(s) substr(s, global_pos, global_pos))
  
  # Map AA to effect scores
  setkey(effects_dt, residue)
  effect_scores <- effects_dt[tip_aa, Effect]
  names(effect_scores) <- tree$tip.label
  
  # Identify tips with valid data (non-gap AA, known effect, non-NA phenotype)
  valid_aa <- tip_aa != "-" & tip_aa != "X" & tip_aa != "B" & tip_aa != "J" & tip_aa != "?"
  valid_effect <- !is.na(effect_scores)
  valid_pheno <- !is.na(tip_pheno[tree$tip.label])
  valid <- valid_aa & valid_effect & valid_pheno
  
  n_valid <- sum(valid)
  if (n_valid < 20) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      n_residues = NA_integer_,
      cor = NA_real_,
      pval = NA_real_,
      status = "insufficient_tips"
    ))
  }
  
  # Count unique residues with effects
  n_residues <- uniqueN(tip_aa[valid])
  if (n_residues < 2) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      n_residues = n_residues,
      cor = NA_real_,
      pval = NA_real_,
      status = "monomorphic"
    ))
  }
  
  # Prune tree to valid tips
  tips_to_drop <- tree$tip.label[!valid]
  if (length(tips_to_drop) > 0) {
    tree_sub <- drop.tip(tree, tips_to_drop)
  } else {
    tree_sub <- tree
  }
  
  # Subset data to match pruned tree
  effect_sub <- effect_scores[tree_sub$tip.label]
  pheno_sub <- tip_pheno[tree_sub$tip.label]
  
  # Check for variance
  if (var(effect_sub, na.rm = TRUE) < 1e-10 || var(pheno_sub, na.rm = TRUE) < 1e-10) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      n_residues = n_residues,
      cor = NA_real_,
      pval = NA_real_,
      status = "no_variance"
    ))
  }
  
  # Compute PICs
  pic_effect <- tryCatch(
    pic(effect_sub, tree_sub, scaled = TRUE, var.contrasts = FALSE),
    error = function(e) NULL
  )
  pic_pheno <- tryCatch(
    pic(pheno_sub, tree_sub, scaled = TRUE, var.contrasts = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(pic_effect) || is.null(pic_pheno)) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      n_residues = n_residues,
      cor = NA_real_,
      pval = NA_real_,
      status = "pic_error"
    ))
  }
  
  # Correlation through origin (standard PIC regression)
  # Under Brownian motion, PICs should have mean 0 and be uncorrelated under null
  ct <- cor.test(pic_effect, pic_pheno)
  
  list(
    global_pos = global_pos,
    n_tips = n_valid,
    n_residues = n_residues,
    n_contrasts = length(pic_effect),
    cor = ct$estimate,
    pval = ct$p.value,
    status = "success"
  )
}

# ==== RUN PIC ANALYSIS FOR ALL GWAS SITES ====

# Get unique sites with their effects
sites <- unique(gwas_effects$GlobalPos)
cat("\nRunning PIC analysis for", length(sites), "sites...\n")

# Prepare effects lookup by site
effects_by_site <- split(gwas_effects[, .(residue, Effect)], gwas_effects$GlobalPos)

# Progress tracking
pb <- txtProgressBar(min = 0, max = length(sites), style = 3)

pic_results <- rbindlist(lapply(seq_along(sites), function(i) {
  pos <- sites[i]
  effects_dt <- effects_by_site[[as.character(pos)]]
  
  result <- pic_site_analysis(
    global_pos = pos,
    tree = tree,
    tip_seqs = tip_seqs_list,
    tip_pheno = tip_pheno,
    effects_dt = effects_dt
  )
  
  if (i %% 100 == 0) setTxtProgressBar(pb, i)
  
  as.data.table(result)
}), fill = TRUE)
close(pb)


summary(pic_results)

plot(hist(pic_results$cor, breaks = 50, main = "Distribution of PIC correlations", xlab = "PIC correlation"))
plot(hist(-log10(pic_results$pval), breaks = 50, main = "Distribution of PIC p-values", xlab = "-log10(p-value)"))

# 1. Is mean cor significantly > 0 across all sites?
t.test(pic_results$cor, mu = 0, na.rm = TRUE)

# ==== ADD GWAS SIGNIFICANCE CLASS ====

# Get site-level p-values
site_pvals <- unique(gwas_effects[, .(GlobalPos, sig_class)])
pic_results <- merge(pic_results, site_pvals, by.x = "global_pos", by.y = "GlobalPos", all.x = TRUE)


# 2. Stratify by GWAS significance - this is the key test
pic_results[status == "success", .(
  n = .N,
  mean_cor = mean(cor),
  prop_positive = mean(cor > 0),
  prop_sig = mean(pval < 0.05)
), by = sig_class][order(-mean_cor)]

# 3. Are GWAS-sig sites better than background?
wilcox.test(
  pic_results[sig_class == "sig_both", cor],
  pic_results[sig_class == "not_sig", cor],
  alternative = "greater"
)

# ==== SUMMARIZE RESULTS ====

cat("\n=== PIC Analysis Summary ===\n")
cat("Sites analyzed:", nrow(pic_results), "\n")
cat("Successful:", sum(pic_results$status == "success"), "\n")

cat("\nStatus breakdown:\n")
print(pic_results[, .N, by = status][order(-N)])

cat("\nSignificance class distribution:\n")
print(pic_results[status == "success", .N, by = sig_class][order(-N)])

# ==== KEY TEST: Do GWAS effects predict phenotype? ====

cat("\n=== Core Question: Do GWAS residue effects correlate with phenotype (PIC)? ===\n")

# For each significance class, test if mean correlation differs from 0
pic_summary <- pic_results[status == "success", .(
  n_sites = .N,
  mean_cor = mean(cor, na.rm = TRUE),
  se_cor = sd(cor, na.rm = TRUE) / sqrt(.N),
  median_cor = median(cor, na.rm = TRUE),
  prop_positive = mean(cor > 0, na.rm = TRUE),
  prop_sig = mean(pval < 0.05, na.rm = TRUE)
), by = sig_class]

pic_summary[, `:=`(
  ci_lo = mean_cor - 1.96 * se_cor,
  ci_hi = mean_cor + 1.96 * se_cor,
  tstat = mean_cor / se_cor,
  pval_vs_zero = 2 * pt(-abs(mean_cor / se_cor), df = n_sites - 1)
)]

print(pic_summary[order(-mean_cor)])

# Binomial test on proportion positive for sig_both
binom.test(
  sum(pic_results[sig_class == "sig_both" & status == "success", cor > 0]),
  pic_results[sig_class == "sig_both" & status == "success", .N],
  p = 0.5
)

# Compare proportions
prop.test(
  c(pic_results[sig_class == "sig_both" & status == "success", sum(cor > 0)],
    pic_results[sig_class == "not_sig" & status == "success", sum(cor > 0)]),
  c(pic_results[sig_class == "sig_both" & status == "success", .N],
    pic_results[sig_class == "not_sig" & status == "success", .N])
)

# ==== META-ANALYSIS: OVERALL EFFECT ====

cat("\n=== Meta-analysis: Overall correlation across sites ===\n")

# Simple meta-analysis: Fisher's method for combining p-values
# But more informative: weighted mean correlation

successful <- pic_results[status == "success"]

# Weight by sqrt(n_contrasts)
successful[, weight := sqrt(n_contrasts)]
weighted_cor <- successful[, sum(cor * weight) / sum(weight)]
cat("Weighted mean correlation:", round(weighted_cor, 4), "\n")

# One-sample t-test: is mean correlation different from 0?
t_test <- t.test(successful$cor, mu = 0)
cat("t-test vs 0: t =", round(t_test$statistic, 3), ", p =", format.pval(t_test$p.value), "\n")

# ==== STRATIFIED: GWAS-significant sites should show stronger effect ====

cat("\n=== Comparison: GWAS-significant vs non-significant sites ===\n")

if (any(pic_results$sig_class == "sig_both")) {
  sig_cors <- successful[sig_class == "sig_both", cor]
  nonsig_cors <- successful[sig_class == "not_sig", cor]
  
  wt <- wilcox.test(sig_cors, nonsig_cors, alternative = "greater")
  cat("Wilcoxon test (sig > nonsig): W =", wt$statistic, ", p =", format.pval(wt$p.value), "\n")
  
  cat("Mean cor (both sig):", round(mean(sig_cors, na.rm = TRUE), 4), "\n")
  cat("Mean cor (not sig):", round(mean(nonsig_cors, na.rm = TRUE), 4), "\n")
}

# ==== SAVE RESULTS ====

saveRDS(pic_results, "pic_results.rds")
cat("\nResults saved to pic_results.rds\n")

# ==== DIAGNOSTIC PLOTS ====

library(ggplot2)

# Histogram of correlations by sig class
p1 <- ggplot(successful, aes(x = cor, fill = sig_class)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~sig_class, scales = "free_y") +
  labs(x = "PIC correlation (effect score ~ phenotype)",
       y = "Number of sites",
       title = "Distribution of PIC correlations by GWAS significance") +
  theme_minimal()

# Volcano-style plot
p2 <- ggplot(successful, aes(x = cor, y = -log10(pval), color = sig_class)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "PIC correlation",
       y = "-log10(p-value)",
       title = "Site-level PIC correlations") +
  theme_minimal()

p2
# Mean correlation by sig class with CIs
p3 <- ggplot(pic_summary, aes(x = reorder(sig_class, -mean_cor), y = mean_cor)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "GWAS significance class",
       y = "Mean PIC correlation",
       title = "Mean correlation by GWAS significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p3
# Save plots

ggsave("pic_histogram.png", p1, width = 10, height = 6, dpi = 150)
ggsave("pic_volcano.png", p2, width = 8, height = 6, dpi = 150)
ggsave("pic_summary.png", p3, width = 6, height = 5, dpi = 150)

cat("\nPlots saved.\n")