#!/usr/bin/env Rscript
# ==============================================================================
# Sister Clade Comparison: Most Rigorous Phylogenetic Test
# 
# Key insight: For each substitution on the tree, we can compare the mean 
# phenotype of the descendant clade vs its sister clade. These are phylogenetically
# independent comparisons that control for shared ancestry.
#
# This is essentially a tree-based analog of PICs but focused on discrete changes.
# ==============================================================================

library(ape)
library(data.table)
library(phangorn)
library(arrow)

# ==== CONFIGURATION ====
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
states_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralStates"
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

stopifnot(file.exists(tree_file), file.exists(states_file), file.exists(aln_file))

# ==== LOAD TREE AND PHENOTYPE DATA ====

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
n_internal <- Nnode(tree)
cat("Tree:", n_tips, "tips,", n_internal, "internal nodes\n")

# Phenotype
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID
tip_pheno <- pheno[tree$tip.label]

n_missing <- sum(is.na(tip_pheno))
cat("Tips with missing phenotype:", n_missing, "of", length(tip_pheno), "\n")

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

# ==== LOAD ASR STATES ====

cat("Loading ASR states...\n")
asr_lines <- readLines(states_file)
asr_data <- strsplit(asr_lines, "\t")
asr_states <- data.table(
  node = sapply(asr_data, `[`, 1),
  sequence = sapply(asr_data, `[`, 2)
)

stopifnot(all(grepl("^Node", asr_states$node)))

# Build seq_lookup: tips = 1:n_tips, internals = (n_tips+1):(n_tips+n_internal)
seq_lookup <- character(n_tips + n_internal)
seq_lookup[1:n_tips] <- tip_seqs_list[tree$tip.label]
internal_nums <- as.integer(gsub("Node", "", asr_states$node))
seq_lookup[n_tips + internal_nums] <- asr_states$sequence

# Check alignment lengths match
stopifnot(nchar(seq_lookup[1]) == nchar(seq_lookup[n_tips + 1]))
cat("Sequences loaded. Alignment length:", nchar(seq_lookup[1]), "\n")

# ==== LOAD GWAS RESULTS ====

cat("Loading GWAS results...\n")
partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

# Load and flatten GWAS effects - matching PIC script exactly
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

# ==== DEFINE SIGNIFICANCE CLASSES ====

# Thresholds based on unique site-level p-values (not residue-level)
site_pvals <- unique(gwas_effects[, .(Gene, GenePos, GlobalPos, P_aa_with_pcs, P_aa_only)])
thresh_control <- quantile(site_pvals$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(site_pvals$P_aa_only, 0.20, na.rm = TRUE)

cat("Threshold (control):", thresh_control, "\n")
cat("Threshold (no control):", thresh_nocontrol, "\n")

# Site-level sig_class
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

# Variant-level sig_class: demote if P_residue >= 0.05
gwas_effects[, variant_sig_class := sig_class]
gwas_effects[P_residue >= 0.05, variant_sig_class := "not_sig"]

cat("\nSite-level sig_class distribution:\n")
print(table(gwas_effects$sig_class))

cat("\nVariant-level sig_class distribution:\n")
print(table(gwas_effects$variant_sig_class))

# Create lookup keyed by GlobalPos and residue
setkey(gwas_effects, GlobalPos, residue)

# ==== IDENTIFY SISTER CLADES ====

#' For each internal node, identify its two child clades
get_sister_info <- function(tree) {
  edge <- tree$edge
  
  # For each internal node, find its children
  internal_nodes <- unique(edge[, 1])
  
  sister_info <- rbindlist(lapply(internal_nodes, function(parent) {
    children <- edge[edge[, 1] == parent, 2]
    
    if (length(children) != 2) return(NULL)  # Skip non-bifurcating
    
    data.table(
      parent_node = parent,
      child1 = children[1],
      child2 = children[2]
    )
  }))
  
  sister_info
}

sister_pairs <- get_sister_info(tree)
cat("\nSister pairs identified:", nrow(sister_pairs), "\n")

# ==== GET CLADE TIP PHENOTYPES ====

#' For each node, get mean phenotype of all descendant tips
get_clade_pheno <- function(node, tree, tip_pheno) {
  if (node <= Ntip(tree)) {
    # It's a tip
    return(tip_pheno[tree$tip.label[node]])
  }
  
  desc_tips <- Descendants(tree, node, type = "tips")[[1]]
  mean(tip_pheno[tree$tip.label[desc_tips]], na.rm = TRUE)
}

cat("Computing clade phenotypes...\n")
all_nodes <- unique(c(sister_pairs$child1, sister_pairs$child2))
clade_phenos <- sapply(all_nodes, function(n) get_clade_pheno(n, tree, tip_pheno))
names(clade_phenos) <- all_nodes

sister_pairs[, pheno1 := clade_phenos[as.character(child1)]]
sister_pairs[, pheno2 := clade_phenos[as.character(child2)]]

# ==== EXTRACT SUBSTITUTIONS WITH SISTER CONTEXT ====

#' For each position, find edges where a substitution occurred
#' Then identify the sister clade and compute phenotype difference

extract_subs_with_sisters <- function(pos, seq_lookup, tree, sister_pairs, gwas_effects) {
  
  edge <- tree$edge
  n_edges <- nrow(edge)
  
  results <- list()
  
  for (i in seq_len(n_edges)) {
    parent <- edge[i, 1]
    child <- edge[i, 2]
    
    if (seq_lookup[parent] == "" || seq_lookup[child] == "") next
    
    parent_aa <- substr(seq_lookup[parent], pos, pos)
    child_aa <- substr(seq_lookup[child], pos, pos)
    
    # Skip if same or gap/ambiguous
    if (parent_aa == child_aa) next
    if (parent_aa %in% c("-", "X", "B", "J", "?")) next
    if (child_aa %in% c("-", "X", "B", "J", "?")) next
    
    # Get GWAS effects for these residues
    from_row <- gwas_effects[.(pos, parent_aa)]
    to_row <- gwas_effects[.(pos, child_aa)]
    
    if (nrow(from_row) == 0 || nrow(to_row) == 0) next
    if (is.na(from_row$Effect[1]) || is.na(to_row$Effect[1])) next
    
    # Find sister clade
    parent_row <- sister_pairs[parent_node == parent]
    
    if (nrow(parent_row) == 0) next
    
    if (parent_row$child1 == child) {
      sister_node <- parent_row$child2
      sister_pheno <- parent_row$pheno2
      focal_pheno <- parent_row$pheno1
    } else {
      sister_node <- parent_row$child1
      sister_pheno <- parent_row$pheno1
      focal_pheno <- parent_row$pheno2
    }
    
    results[[length(results) + 1]] <- data.table(
      position = pos,
      edge_idx = i,
      parent_node = parent,
      child_node = child,
      from_aa = parent_aa,
      to_aa = child_aa,
      from_effect = from_row$Effect[1],
      to_effect = to_row$Effect[1],
      from_P_residue = from_row$P_residue[1],
      to_P_residue = to_row$P_residue[1],
      effect_change = to_row$Effect[1] - from_row$Effect[1],
      focal_pheno = focal_pheno,
      sister_pheno = sister_pheno,
      pheno_contrast = focal_pheno - sister_pheno,
      # Site-level classification
      sig_class = from_row$sig_class[1],
      # Variant-level: use the 'to' residue's classification
      variant_sig_class = to_row$variant_sig_class[1]
    )
  }
  
  if (length(results) == 0) return(NULL)
  rbindlist(results)
}

# ==== RUN ANALYSIS ====

gwas_positions <- unique(gwas_effects$GlobalPos)
cat("\nExtracting substitutions at", length(gwas_positions), "GWAS positions...\n")

pb <- txtProgressBar(min = 0, max = length(gwas_positions), style = 3)

sister_subs <- rbindlist(lapply(seq_along(gwas_positions), function(i) {
  pos <- gwas_positions[i]
  result <- extract_subs_with_sisters(pos, seq_lookup, tree, sister_pairs, gwas_effects)
  if (i %% 500 == 0) setTxtProgressBar(pb, i)
  result
}), fill = TRUE)
close(pb)

saveRDS(sister_subs, "data/tmp/sister_subs_intermediate.rds")

cat("\nSubstitutions with sister context:", nrow(sister_subs), "\n")

# ==== SUMMARIZE BY SIGNIFICANCE CLASS ====

cat("\n=== Substitutions by site-level sig_class ===\n")
print(sister_subs[, .N, by = sig_class][order(-N)])

cat("\n=== Substitutions by variant-level sig_class ===\n")
print(sister_subs[, .N, by = variant_sig_class][order(-N)])

# ==== CORE TEST: Effect change vs phenotype contrast ====

cat("\n=== Sister Clade Analysis ===\n")

# Key prediction: if a substitution increases effect (to_effect > from_effect),
# the focal clade should be warmer than its sister (pheno_contrast > 0)

cor_test <- cor.test(sister_subs$effect_change, sister_subs$pheno_contrast)
cat("\nOverall correlation (effect_change vs pheno_contrast):\n")
cat("r =", round(cor_test$estimate, 4), ", p =", format.pval(cor_test$p.value), "\n")

# Stratified by site-level GWAS class
cat("\n=== Correlation by site-level sig_class ===\n")
strat_cor_site <- sister_subs[, {
  if (.N < 50) list(n = .N, cor = NA_real_, pval = NA_real_)
  else {
    ct <- cor.test(effect_change, pheno_contrast)
    list(n = .N, cor = ct$estimate, pval = ct$p.value)
  }
}, by = sig_class]

strat_cor_site[, `:=`(
  ci_lo = cor - 1.96 * sqrt((1 - cor^2) / (n - 2)),
  ci_hi = cor + 1.96 * sqrt((1 - cor^2) / (n - 2))
)]

print(strat_cor_site[order(-cor)])

plot(sister_subs$effect_change[sister_subs$variant_sig_class=="sig_both" &
  abs(sister_subs$effect_change) > 10 & abs(sister_subs$pheno_contrast) > 2],
     sister_subs$pheno_contrast[sister_subs$variant_sig_class=="sig_both" &
  abs(sister_subs$effect_change) > 10 & abs(sister_subs$pheno_contrast) > 2])

specific_dt <- sister_subs[sister_subs$variant_sig_class=="sig_both" &
  abs(sister_subs$effect_change) > 10 & abs(sister_subs$pheno_contrast) > 2,]

summary(sister_subs$effect_change[sister_subs$variant_sig_class=="sig_both" &
  abs(sister_subs$effect_change) > 10])

unique(sister_subs$variant_sig_class)
#[1] "not_sig"       "sig_nocontrol" "sig_control"   "sig_both"  
par(mfrow=c(2,2), mar=c(4,4,2,1), oma=c(0,0,2,0))

plot(sister_subs$effect_change[sister_subs$variant_sig_class=="sig_both"],
     sister_subs$pheno_contrast[sister_subs$variant_sig_class=="sig_both"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_both")

plot(sister_subs$effect_change[sister_subs$variant_sig_class=="sig_control"],
     sister_subs$pheno_contrast[sister_subs$variant_sig_class=="sig_control"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_control")

plot(sister_subs$effect_change[sister_subs$variant_sig_class=="sig_nocontrol"],
     sister_subs$pheno_contrast[sister_subs$variant_sig_class=="sig_nocontrol"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_nocontrol")

plot(sister_subs$effect_change[sister_subs$variant_sig_class=="not_sig"],
     sister_subs$pheno_contrast[sister_subs$variant_sig_class=="not_sig"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="not_sig")

mtext("Effect change vs phenotype contrast", outer=TRUE, cex=1.1)
par(mfrow=c(1,1))


good_sigs <- sister_subs[from_P_residue < 0.05 & to_P_residue < 0.05,]

par(mfrow=c(2,2), mar=c(4,4,2,1), oma=c(0,0,2,0))

plot(good_sigs$to_effect[good_sigs$variant_sig_class=="sig_both"],
     good_sigs$pheno_contrast[good_sigs$variant_sig_class=="sig_both"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_both")

plot(good_sigs$to_effect[good_sigs$variant_sig_class=="sig_control"],
     good_sigs$pheno_contrast[good_sigs$variant_sig_class=="sig_control"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_control")

plot(good_sigs$to_effect[good_sigs$variant_sig_class=="sig_nocontrol"],
     good_sigs$pheno_contrast[good_sigs$variant_sig_class=="sig_nocontrol"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="sig_nocontrol")

plot(good_sigs$to_effect[good_sigs$variant_sig_class=="not_sig"],
     good_sigs$pheno_contrast[good_sigs$variant_sig_class=="not_sig"],
     pch=16, cex=0.3,
     xlab="Effect change", ylab="Phenotype contrast",
     main="not_sig")

mtext("Effect change vs phenotype contrast (only consider significant gwas estimates)", outer=TRUE, cex=1.1)
par(mfrow=c(1,1))




# Stratified by variant-level sig_class
cat("\n=== Correlation by variant-level sig_class ===\n")
strat_cor_var <- good_sigs[, {
  if (.N < 50) list(n = .N, cor = NA_real_, pval = NA_real_)
  else {
    ct <- cor.test(effect_change, pheno_contrast)
    list(n = .N, cor = ct$estimate, pval = ct$p.value)
  }
}, by = variant_sig_class]

strat_cor_var[, `:=`(
  ci_lo = cor - 1.96 * sqrt((1 - cor^2) / (n - 2)),
  ci_hi = cor + 1.96 * sqrt((1 - cor^2) / (n - 2))
)]

print(strat_cor_var[order(-cor)])

# ==== SIGN CONCORDANCE ====

# Simpler test: does the sign match?
effect_thresh <- sd(sister_subs$effect_change, na.rm = TRUE) * 0.5
pheno_thresh <- sd(sister_subs$pheno_contrast, na.rm = TRUE) * 0.5

cat("\nEffect threshold (0.5 SD):", round(effect_thresh, 4), "\n")
cat("Pheno threshold (0.5 SD):", round(pheno_thresh, 4), "\n")

sister_subs[, effect_dir := fcase(
  effect_change > effect_thresh, "warm",
  effect_change < -effect_thresh, "cold",
  default = "neutral"
)]

sister_subs[, pheno_dir := fcase(
  pheno_contrast > pheno_thresh, "warm",
  pheno_contrast < -pheno_thresh, "cold",
  default = "neutral"
)]

cat("\n=== Direction distributions ===\n")
cat("Effect direction:\n")
print(table(sister_subs$effect_dir))
cat("\nPheno direction:\n")
print(table(sister_subs$pheno_dir))

# Filter to directional changes
sister_dir <- sister_subs[effect_dir != "neutral" & pheno_dir != "neutral"]
cat("\nDirectional substitutions:", nrow(sister_dir), "\n")

cat("\n=== Sign concordance (directional only) ===\n")
conc_table <- table(effect = sister_dir$effect_dir, pheno = sister_dir$pheno_dir)
print(conc_table)

concordance <- sum(diag(conc_table)) / sum(conc_table)
cat("\nOverall concordance:", round(concordance, 4), "\n")

# Binomial test vs 0.5
n_concordant <- sum(diag(conc_table))
n_total <- sum(conc_table)
bt <- binom.test(n_concordant, n_total, p = 0.5)
cat("Binomial test vs 0.5: p =", format.pval(bt$p.value), "\n")

# ==== STRATIFIED CONCORDANCE ====

# By site-level sig_class
cat("\n=== Concordance by site-level sig_class ===\n")
conc_by_site_class <- sister_dir[, {
  ct <- table(effect_dir, pheno_dir)
  n_conc <- sum(diag(ct))
  n_tot <- sum(ct)
  prop <- n_conc / n_tot
  bt <- if (n_tot >= 10) binom.test(n_conc, n_tot, p = 0.5) else list(p.value = NA, conf.int = c(NA, NA))
  list(
    n = n_tot,
    concordance = prop,
    pval = bt$p.value,
    ci_lo = bt$conf.int[1],
    ci_hi = bt$conf.int[2]
  )
}, by = sig_class]

print(conc_by_site_class[order(-concordance)])

# By variant-level sig_class
cat("\n=== Concordance by variant-level sig_class ===\n")
conc_by_var_class <- sister_dir[, {
  ct <- table(effect_dir, pheno_dir)
  n_conc <- sum(diag(ct))
  n_tot <- sum(ct)
  prop <- n_conc / n_tot
  bt <- if (n_tot >= 10) binom.test(n_conc, n_tot, p = 0.5) else list(p.value = NA, conf.int = c(NA, NA))
  list(
    n = n_tot,
    concordance = prop,
    pval = bt$p.value,
    ci_lo = bt$conf.int[1],
    ci_hi = bt$conf.int[2]
  )
}, by = variant_sig_class]

print(conc_by_var_class[order(-concordance)])

# ==== COMPARISON: GWAS-sig vs background ====

cat("\n=== Is GWAS-significant better than background? ===\n")

# Site-level comparison
if (any(sister_dir$sig_class == "sig_both")) {
  sig_conc <- sister_dir[sig_class == "sig_both", mean(effect_dir == pheno_dir)]
  nonsig_conc <- sister_dir[sig_class == "not_sig", mean(effect_dir == pheno_dir)]
  
  sig_n <- sister_dir[sig_class == "sig_both", .N]
  nonsig_n <- sister_dir[sig_class == "not_sig", .N]
  
  sig_succ <- round(sig_conc * sig_n)
  nonsig_succ <- round(nonsig_conc * nonsig_n)
  
  prop_test <- prop.test(c(sig_succ, nonsig_succ), c(sig_n, nonsig_n), alternative = "greater")
  
  cat("\nSite-level comparison:\n")
  cat("sig_both concordance:", round(sig_conc, 4), "(n=", sig_n, ")\n")
  cat("not_sig concordance:", round(nonsig_conc, 4), "(n=", nonsig_n, ")\n")
  cat("Prop test (sig > nonsig): p =", format.pval(prop_test$p.value), "\n")
}

# Variant-level comparison
if (any(sister_dir$variant_sig_class == "sig_both")) {
  sig_conc_v <- sister_dir[variant_sig_class == "sig_both", mean(effect_dir == pheno_dir)]
  nonsig_conc_v <- sister_dir[variant_sig_class == "not_sig", mean(effect_dir == pheno_dir)]
  
  sig_n_v <- sister_dir[variant_sig_class == "sig_both", .N]
  nonsig_n_v <- sister_dir[variant_sig_class == "not_sig", .N]
  
  sig_succ_v <- round(sig_conc_v * sig_n_v)
  nonsig_succ_v <- round(nonsig_conc_v * nonsig_n_v)
  
  prop_test_v <- prop.test(c(sig_succ_v, nonsig_succ_v), c(sig_n_v, nonsig_n_v), alternative = "greater")
  
  cat("\nVariant-level comparison:\n")
  cat("sig_both concordance:", round(sig_conc_v, 4), "(n=", sig_n_v, ")\n")
  cat("not_sig concordance:", round(nonsig_conc_v, 4), "(n=", nonsig_n_v, ")\n")
  cat("Prop test (sig > nonsig): p =", format.pval(prop_test_v$p.value), "\n")
}

# ==== WILCOXON TESTS ON CORRELATIONS ====

cat("\n=== Wilcoxon tests: sig_both vs not_sig ===\n")

# Site-level
wt_site <- wilcox.test(
  sister_subs[sig_class == "sig_both", effect_change * sign(pheno_contrast)],
  sister_subs[sig_class == "not_sig", effect_change * sign(pheno_contrast)],
  alternative = "greater"
)
cat("Site-level Wilcoxon p:", format.pval(wt_site$p.value), "\n")

# Variant-level
wt_var <- wilcox.test(
  sister_subs[variant_sig_class == "sig_both", effect_change * sign(pheno_contrast)],
  sister_subs[variant_sig_class == "not_sig", effect_change * sign(pheno_contrast)],
  alternative = "greater"
)
cat("Variant-level Wilcoxon p:", format.pval(wt_var$p.value), "\n")

# ==== SAVE RESULTS ====

cat("\nResults saved to data/tmp/sister_subs_results.rds\n")

# ==== VISUALIZATION ====

library(ggplot2)

# Scatter plot of effect_change vs pheno_contrast by site-level sig_class
p1 <- ggplot(sister_subs, aes(x = effect_change, y = pheno_contrast)) +
  geom_hex(bins = 50) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10") +
  facet_wrap(~sig_class) +
  labs(x = "GWAS effect change (to - from)",
       y = "Phenotype contrast (focal - sister clade)",
       title = "Sister clade comparison by site-level significance") +
  theme_minimal()
p1
# By variant-level sig_class
p2 <- ggplot(sister_subs, aes(x = effect_change, y = pheno_contrast)) +
  geom_hex(bins = 50) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10") +
  facet_wrap(~variant_sig_class) +
  labs(x = "GWAS effect change (to - from)",
       y = "Phenotype contrast (focal - sister clade)",
       title = "Sister clade comparison by variant-level significance") +
  theme_minimal()
p2
# Concordance bar plot - site level
p3 <- ggplot(conc_by_site_class, aes(x = reorder(sig_class, -concordance), y = concordance)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.45, 0.65), labels = scales::percent_format()) +
  labs(x = "Site-level GWAS significance",
       y = "Sign concordance",
       title = "Sister clade concordance by site significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p3
# Concordance bar plot - variant level
p4 <- ggplot(conc_by_var_class, aes(x = reorder(variant_sig_class, -concordance), y = concordance)) +
  geom_col(fill = "darkgreen") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.45, 0.65), labels = scales::percent_format()) +
  labs(x = "Variant-level GWAS significance",
       y = "Sign concordance",
       title = "Sister clade concordance by variant significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p4
# Save plots
ggsave("sister_scatter_site.png", p1, width = 10, height = 8, dpi = 150)
ggsave("sister_scatter_variant.png", p2, width = 10, height = 8, dpi = 150)
ggsave("sister_concordance_site.png", p3, width = 7, height = 5, dpi = 150)
ggsave("sister_concordance_variant.png", p4, width = 7, height = 5, dpi = 150)

cat("\nPlots saved.\n")