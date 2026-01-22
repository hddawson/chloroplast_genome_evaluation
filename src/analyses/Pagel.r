#!/usr/bin/env Rscript
# ==============================================================================
# Classical Convergence Test: Pagel's Discrete Correlation
# 
# Tests whether amino acid state and phenotype evolve in a correlated fashion
# using Pagel (1994) method. This is the classical test for correlated evolution
# of two discrete traits on a phylogeny.
#
# For each GWAS site, we ask: do AA transitions depend on phenotypic state?
# ==============================================================================

library(ape)
library(data.table)
library(phytools)
library(arrow)

# ==== CONFIGURATION ====
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

stopifnot(file.exists(tree_file), file.exists(aln_file), file.exists(data_file))

# ==== LOAD TREE AND DATA ====

cat("Loading tree and data...\n")
tree <- read.tree(tree_file)
data <- read_parquet(data_file)
setDT(data)

# Root on Pinales
pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree$tip.label)
tree <- root(tree, outgroup = pinales_in_tree, resolve.root = TRUE)
stopifnot(is.rooted(tree))

n_tips <- Ntip(tree)
cat("Tree:", n_tips, "tips\n")

# ==== LOAD PHENOTYPE ====

pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID
tip_pheno <- pheno[tree$tip.label]

# --- EDA: Phenotype distribution ---
par(mfrow = c(1, 2))
hist(tip_pheno, breaks = 50, main = "Phenotype distribution", 
     xlab = "Mean temp wettest quarter", col = "steelblue")
abline(v = median(tip_pheno, na.rm = TRUE), col = "red", lwd = 2)

# Discretize phenotype for Pagel's test
pheno_median <- median(tip_pheno, na.rm = TRUE)
tip_pheno_binary <- ifelse(tip_pheno > pheno_median, "warm", "cold")
tip_pheno_binary[is.na(tip_pheno_binary)] <- NA

barplot(table(tip_pheno_binary), main = "Discretized phenotype", 
        col = c("dodgerblue", "tomato"), ylab = "N tips")
par(mfrow = c(1, 1))

cat("Phenotype median:", round(pheno_median, 2), "\n")
cat("Warm tips:", sum(tip_pheno_binary == "warm", na.rm = TRUE), "\n")
cat("Cold tips:", sum(tip_pheno_binary == "cold", na.rm = TRUE), "\n")

# ==== LOAD TIP SEQUENCES ====

cat("\nLoading alignment...\n")
tip_aln <- read.FASTA(aln_file, type = "AA")
tip_seqs_list <- sapply(tip_aln, function(x) paste(rawToChar(x, multiple = TRUE), collapse = ""))
names(tip_seqs_list) <- names(tip_aln)

stopifnot(all(tree$tip.label %in% names(tip_seqs_list)))

# Pad sequences
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

cat("\nLoading GWAS results...\n")
partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_effects <- rbindlist(lapply(model_files, function(f) {
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

gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

cat("GWAS effects loaded:", nrow(gwas_effects), "residue-site combinations\n")
cat("Unique sites:", uniqueN(gwas_effects$GlobalPos), "\n")

# --- EDA: GWAS effect distribution ---
par(mfrow = c(1, 2))
hist(gwas_effects$Effect, breaks = 50, main = "GWAS residue effects", 
     xlab = "Effect size", col = "darkgreen")
abline(v = 0, col = "red", lwd = 2)

hist(-log10(gwas_effects$P_aa_with_pcs), breaks = 50, 
     main = "Site-level significance", xlab = "-log10(P)", col = "purple")
par(mfrow = c(1, 1))

# ==== DEFINE SIGNIFICANCE CLASSES ====

site_pvals <- unique(gwas_effects[, .(Gene, GenePos, GlobalPos, P_aa_with_pcs, P_aa_only)])
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

cat("\nSig class distribution:\n")
print(table(gwas_effects$sig_class))

# ==== HELPER: GET AA STATE FOR A SITE ====

get_site_aa <- function(pos, tip_seqs, tree) {
  aa <- sapply(tip_seqs[tree$tip.label], function(s) substr(s, pos, pos))
  names(aa) <- tree$tip.label
  aa
}

# ==== HELPER: BINARIZE AA STATE ====
# For Pagel's test, we need binary traits
# Strategy: "warm-associated" AA vs "cold-associated" AA based on GWAS effects

binarize_aa <- function(aa_vector, effects_dt) {
  # Get effect for each AA
  setkey(effects_dt, residue)
  
  # Identify warm vs cold AAs based on effect sign
  aa_effects <- effects_dt[, .(mean_effect = mean(Effect, na.rm = TRUE)), by = residue]
  warm_aas <- aa_effects[mean_effect > 0, residue]
  cold_aas <- aa_effects[mean_effect <= 0, residue]
  
  # Binarize
  aa_binary <- fcase(
    aa_vector %in% warm_aas, "warm_aa",
    aa_vector %in% cold_aas, "cold_aa",
    default = NA_character_
  )
  names(aa_binary) <- names(aa_vector)
  
  list(binary = aa_binary, warm_aas = warm_aas, cold_aas = cold_aas)
}

# ==== PAGEL'S TEST FOR ONE SITE ====

#' Run Pagel's test for correlated evolution at a single site
#' 
#' Tests if AA state and phenotype evolve independently (model "independent")
#' or in a correlated fashion (model "dependent")
#'
#' @param global_pos Position in alignment
#' @param tree Phylo object
#' @param tip_seqs Named character vector
#' @param tip_pheno_binary Named character vector ("warm"/"cold")
#' @param effects_dt data.table with residue effects
#' @return List with test results

pagel_site_test <- function(global_pos, tree, tip_seqs, tip_pheno_binary, effects_dt) {
  
  # Get AA at this site
  aa_raw <- get_site_aa(global_pos, tip_seqs, tree)
  
  # Binarize based on GWAS effects
  aa_bin_result <- binarize_aa(aa_raw, effects_dt)
  aa_binary <- aa_bin_result$binary
  
  # Identify valid tips (non-gap AA, non-NA phenotype)
  valid_aa <- !is.na(aa_binary)
  valid_pheno <- !is.na(tip_pheno_binary[tree$tip.label])
  valid <- valid_aa & valid_pheno
  
  n_valid <- sum(valid)
  
  # Need sufficient tips and variation in both traits
  if (n_valid < 50) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      status = "insufficient_tips",
      p_pagel = NA_real_,
      lr = NA_real_,
      aic_ind = NA_real_,
      aic_dep = NA_real_
    ))
  }
  
  # Check for variation
  aa_table <- table(aa_binary[valid])
  pheno_table <- table(tip_pheno_binary[tree$tip.label][valid])
  
  if (length(aa_table) < 2 || min(aa_table) < 10) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      status = "insufficient_aa_variation",
      p_pagel = NA_real_,
      lr = NA_real_,
      aic_ind = NA_real_,
      aic_dep = NA_real_
    ))
  }
  
  if (length(pheno_table) < 2 || min(pheno_table) < 10) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      status = "insufficient_pheno_variation",
      p_pagel = NA_real_,
      lr = NA_real_,
      aic_ind = NA_real_,
      aic_dep = NA_real_
    ))
  }
  
  # Prune tree to valid tips
  tips_to_drop <- tree$tip.label[!valid]
  if (length(tips_to_drop) > 0) {
    tree_sub <- drop.tip(tree, tips_to_drop)
  } else {
    tree_sub <- tree
  }
  
  # Prepare trait vectors (must be named, match tree)
  x <- aa_binary[tree_sub$tip.label]
  y <- tip_pheno_binary[tree_sub$tip.label]
  
  # Convert to factor
  x <- as.factor(x)
  y <- as.factor(y)
  
  # Run Pagel's test
  fit <- tryCatch({
    fitPagel(tree_sub, x, y, method = "fitMk")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(list(
      global_pos = global_pos,
      n_tips = n_valid,
      status = "fit_error",
      p_pagel = NA_real_,
      lr = NA_real_,
      aic_ind = NA_real_,
      aic_dep = NA_real_
    ))
  }
  
  # Extract results
  list(
    global_pos = global_pos,
    n_tips = n_valid,
    n_warm_aa = aa_table["warm_aa"],
    n_cold_aa = aa_table["cold_aa"],
    status = "success",
    p_pagel = fit$P,
    lr = fit$LR,
    aic_ind = fit$independent.AIC,
    aic_dep = fit$dependent.AIC,
    delta_aic = fit$independent.AIC - fit$dependent.AIC  # positive = dependent better
  )
}

# ==== RUN PAGEL'S TEST FOR ALL GWAS SITES ====

# For speed, focus on sites with GWAS signal first
# Then sample from non-significant for comparison

sig_sites <- unique(gwas_effects[sig_class == "sig_both", GlobalPos])
control_sites <- unique(gwas_effects[sig_class == "sig_control", GlobalPos])
nocontrol_sites <- unique(gwas_effects[sig_class == "sig_nocontrol", GlobalPos])
nonsig_sites <- unique(gwas_effects[sig_class == "not_sig", GlobalPos])

cat("\nSites by class:\n")
cat("sig_both:", length(sig_sites), "\n")
cat("sig_control:", length(control_sites), "\n")
cat("sig_nocontrol:", length(nocontrol_sites), "\n")
cat("not_sig:", length(nonsig_sites), "\n")

# Sample nonsig for tractability (Pagel's test is slow)
set.seed(42)
n_nonsig_sample <- min(100, length(nonsig_sites))
nonsig_sample <- sample(nonsig_sites, n_nonsig_sample)

all_test_sites <- c(sig_sites, control_sites, nocontrol_sites, nonsig_sample)
cat("\nTotal sites to test:", length(all_test_sites), "\n")

# Prepare effects lookup
effects_by_site <- split(gwas_effects[, .(residue, Effect)], gwas_effects$GlobalPos)

# Progress tracking
cat("\nRunning Pagel's tests (this may take a while)...\n")
pb <- txtProgressBar(min = 0, max = length(all_test_sites), style = 3)

pagel_results <- rbindlist(lapply(seq_along(all_test_sites), function(i) {
  pos <- all_test_sites[i]
  effects_dt <- effects_by_site[[as.character(pos)]]
  
  result <- pagel_site_test(
    global_pos = pos,
    tree = tree,
    tip_seqs = tip_seqs_list,
    tip_pheno_binary = tip_pheno_binary,
    effects_dt = effects_dt
  )
  
  if (i %% 10 == 0) setTxtProgressBar(pb, i)
  
  as.data.table(result)
}), fill = TRUE)
close(pb)

saveRDS(pagel_results, "data/tmp/pagel_results.rds")
cat("SAVEDEM!")

# ==== ADD SIGNIFICANCE CLASS ====

site_class_lookup <- unique(gwas_effects[, .(GlobalPos, sig_class)])
pagel_results <- merge(pagel_results, site_class_lookup, 
                       by.x = "global_pos", by.y = "GlobalPos", all.x = TRUE)

# Mark sampled nonsig
pagel_results[global_pos %in% nonsig_sample & sig_class == "not_sig", 
              sig_class := "not_sig_sampled"]

cat("\n=== Pagel Results Summary ===\n")
cat("Sites tested:", nrow(pagel_results), "\n")
cat("Successful:", sum(pagel_results$status == "success"), "\n")

# --- EDA: Status breakdown ---
par(mfrow = c(1, 1))
status_tab <- table(pagel_results$status)
barplot(status_tab, main = "Pagel test status", las = 2, col = "lightblue",
        cex.names = 0.7)

# ==== ANALYZE RESULTS ====

successful <- pagel_results[status == "success"]
cat("\nSuccessful tests:", nrow(successful), "\n")

# --- EDA: P-value distribution ---
par(mfrow = c(2, 2))

hist(successful$p_pagel, breaks = 30, main = "Pagel P-values (all sites)",
     xlab = "P-value", col = "lightgreen")
abline(v = 0.05, col = "red", lwd = 2)

hist(-log10(successful$p_pagel), breaks = 30, main = "-log10(P) distribution",
     xlab = "-log10(P)", col = "lightgreen")
abline(v = -log10(0.05), col = "red", lwd = 2)

# Delta AIC distribution
hist(successful$delta_aic, breaks = 30, main = "Delta AIC (Ind - Dep)",
     xlab = "Delta AIC (positive = dependent better)", col = "salmon")
abline(v = 0, col = "blue", lwd = 2)
abline(v = 2, col = "red", lwd = 2, lty = 2)  # AIC threshold

# By sig class
boxplot(delta_aic ~ sig_class, data = successful, 
        main = "Delta AIC by GWAS class",
        col = c("tomato", "orange", "gold", "gray"),
        las = 2, cex.axis = 0.8)
abline(h = 0, col = "blue", lwd = 2)
abline(h = 2, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# ==== STRATIFIED SUMMARY ====

cat("\n=== Results by GWAS significance class ===\n")

pagel_summary <- successful[, .(
  n_sites = .N,
  mean_p = mean(p_pagel, na.rm = TRUE),
  median_p = median(p_pagel, na.rm = TRUE),
  prop_sig = mean(p_pagel < 0.05, na.rm = TRUE),
  mean_delta_aic = mean(delta_aic, na.rm = TRUE),
  prop_aic_support = mean(delta_aic > 2, na.rm = TRUE)  # AIC > 2 favors dependent
), by = sig_class]

print(pagel_summary[order(-prop_sig)])

# --- EDA: Proportion significant by class ---
par(mfrow = c(1, 2))

barplot(pagel_summary[order(-prop_sig), prop_sig], 
        names.arg = pagel_summary[order(-prop_sig), sig_class],
        main = "Prop. significant (P < 0.05)",
        col = "forestgreen", ylim = c(0, max(pagel_summary$prop_sig) * 1.2),
        las = 2, cex.names = 0.8)
abline(h = 0.05, col = "red", lwd = 2, lty = 2)  # expected under null

barplot(pagel_summary[order(-prop_aic_support), prop_aic_support], 
        names.arg = pagel_summary[order(-prop_aic_support), sig_class],
        main = "Prop. AIC supports correlation",
        col = "purple", ylim = c(0, max(pagel_summary$prop_aic_support) * 1.2),
        las = 2, cex.names = 0.8)

par(mfrow = c(1, 1))

# ==== STATISTICAL TESTS ====

cat("\n=== Statistical comparisons ===\n")

# Test: are sig_both sites more likely to show correlated evolution?
if (sum(successful$sig_class == "sig_both") > 10 & 
    sum(successful$sig_class == "not_sig_sampled") > 10) {
  
  # Compare proportion significant
  sig_both_prop <- successful[sig_class == "sig_both", mean(p_pagel < 0.05)]
  nonsig_prop <- successful[sig_class == "not_sig_sampled", mean(p_pagel < 0.05)]
  
  sig_both_n <- successful[sig_class == "sig_both", .N]
  nonsig_n <- successful[sig_class == "not_sig_sampled", .N]
  
  prop_test <- prop.test(
    c(round(sig_both_prop * sig_both_n), round(nonsig_prop * nonsig_n)),
    c(sig_both_n, nonsig_n),
    alternative = "greater"
  )
  
  cat("\nProportion with significant Pagel P:\n")
  cat("sig_both:", round(sig_both_prop, 4), "(n=", sig_both_n, ")\n")
  cat("not_sig:", round(nonsig_prop, 4), "(n=", nonsig_n, ")\n")
  cat("Prop test (sig_both > not_sig): p =", format.pval(prop_test$p.value), "\n")
  
  # Wilcoxon on delta AIC
  wt <- wilcox.test(
    successful[sig_class == "sig_both", delta_aic],
    successful[sig_class == "not_sig_sampled", delta_aic],
    alternative = "greater"
  )
  cat("\nWilcoxon test on delta AIC (sig_both > not_sig): p =", 
      format.pval(wt$p.value), "\n")
  
  # Mean delta AIC comparison
  cat("\nMean delta AIC:\n")
  cat("sig_both:", round(successful[sig_class == "sig_both", mean(delta_aic)], 3), "\n")
  cat("not_sig:", round(successful[sig_class == "not_sig_sampled", mean(delta_aic)], 3), "\n")
}

# ==== CONTINGENCY TABLE ====

cat("\n=== Contingency: GWAS sig vs Pagel sig ===\n")

successful[, pagel_sig := p_pagel < 0.05]
successful[, gwas_sig := sig_class %in% c("sig_both", "sig_control")]

cont_table <- table(
  GWAS = successful$gwas_sig,
  Pagel = successful$pagel_sig
)
print(cont_table)

fisher_test <- fisher.test(cont_table)
cat("\nFisher's exact test p:", format.pval(fisher_test$p.value), "\n")
cat("Odds ratio:", round(fisher_test$estimate, 3), "\n")

# --- EDA: Mosaic plot ---
mosaicplot(cont_table, main = "GWAS sig vs Pagel sig",
           color = c("lightblue", "salmon"), shade = FALSE)

# ==== TOP HITS ====

cat("\n=== Top 20 sites by Pagel significance ===\n")
top_hits <- successful[order(p_pagel)][1:20, .(global_pos, sig_class, n_tips, 
                                               p_pagel, delta_aic)]
print(top_hits)

# --- EDA: Top hits by class ---
top_by_class <- successful[order(p_pagel), .SD[1:5], by = sig_class]
cat("\n=== Top 5 per class ===\n")
print(top_by_class[, .(sig_class, global_pos, p_pagel, delta_aic)])

# ==== SAVE RESULTS ====

saveRDS(pagel_results, "pagel_results.rds")
cat("\nResults saved to pagel_results.rds\n")

# ==== FINAL SUMMARY PLOT ====

par(mfrow = c(2, 2))

# 1. P-value by class
boxplot(p_pagel ~ sig_class, data = successful,
        main = "Pagel P-value by GWAS class",
        col = c("tomato", "orange", "gold", "gray"),
        las = 2, cex.axis = 0.8, ylab = "P-value")
abline(h = 0.05, col = "red", lwd = 2, lty = 2)

# 2. -log10(P) by class  
boxplot(-log10(p_pagel) ~ sig_class, data = successful,
        main = "-log10(P) by GWAS class",
        col = c("tomato", "orange", "gold", "gray"),
        las = 2, cex.axis = 0.8, ylab = "-log10(P)")
abline(h = -log10(0.05), col = "red", lwd = 2, lty = 2)

# 3. Correlation: GWAS effect range vs Pagel significance
site_effect_range <- gwas_effects[, .(effect_range = max(Effect) - min(Effect)), by = GlobalPos]
successful_merged <- merge(successful, site_effect_range, 
                           by.x = "global_pos", by.y = "GlobalPos")

plot(successful_merged$effect_range, -log10(successful_merged$p_pagel),
     pch = 16, col = adjustcolor("steelblue", 0.5),
     xlab = "GWAS effect range at site",
     ylab = "-log10(Pagel P)",
     main = "Effect range vs convergence signal")
abline(h = -log10(0.05), col = "red", lwd = 2, lty = 2)

# Add correlation
cor_test <- cor.test(successful_merged$effect_range, 
                     -log10(successful_merged$p_pagel))
legend("topright", 
       legend = paste0("r = ", round(cor_test$estimate, 3), 
                       "\np = ", format.pval(cor_test$p.value)),
       bty = "n")

# 4. Summary barplot
summary_mat <- as.matrix(pagel_summary[, .(prop_sig, prop_aic_support)])
rownames(summary_mat) <- pagel_summary$sig_class
barplot(t(summary_mat), beside = TRUE,
        col = c("forestgreen", "purple"),
        legend.text = c("P < 0.05", "AIC support"),
        args.legend = list(x = "topright", bty = "n"),
        main = "Convergence evidence by GWAS class",
        las = 2, cex.names = 0.8)

par(mfrow = c(1, 1))

cat("\n=== Analysis complete ===\n")