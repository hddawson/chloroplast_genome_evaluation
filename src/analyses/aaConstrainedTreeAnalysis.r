#!/usr/bin/env Rscript
# ==============================================================================
# GWAS-ASR Substitution Analysis
# Connect GWAS significant sites to phylogenetic substitution patterns from ASR
# ==============================================================================

library(ape)
library(data.table)
library(phangorn)
library(arrow)

# ==== FILE PATHS ====
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
states_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralStates"
aln_file <- "raxml_input/superaa_collapsed.fasta"  # UPDATE: path to input AA alignment
data_file <- "data/processed_data.parquet"
# model_files <- list.files("path/to/gwas", pattern = "*.rds", full.names = TRUE)  # UPDATE

stopifnot(file.exists(tree_file), file.exists(states_file), file.exists(aln_file))

# ==== LOAD TREE ====
tree <- read.tree(tree_file)
data <- read_parquet(data_file)
setDT(data)

pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree$tip.label)
tree <- root(tree, outgroup = pinales_in_tree, resolve.root = TRUE)
stopifnot(is.rooted(tree))

n_tips <- Ntip(tree)
n_internal <- Nnode(tree)
cat("Tree: ", n_tips, "tips,", n_internal, "internal nodes\n")

# ==== LOAD SEQUENCES ====

# Tip sequences from original alignment (AA FASTA)
tip_aln <- read.FASTA(aln_file, type = "AA")
tip_seqs_list <- sapply(tip_aln, function(x) paste(rawToChar(x, multiple = TRUE), collapse = ""))
names(tip_seqs_list) <- names(tip_aln)

stopifnot(all(tree$tip.label %in% names(tip_seqs_list)))
stopifnot(length(unique(nchar(tip_seqs_list))) == 1)  # all same length


# Check length distribution
lens <- nchar(tip_seqs_list)
cat("Length distribution:\n")
print(table(lens))

# Get examples of each length
for (l in unique(lens)) {
  examples <- names(tip_seqs_list)[lens == l][1:2]
  cat("\n--- Length", l, "---\n")
  for (nm in examples) {
    seq <- tip_seqs_list[nm]
    cat(nm, "\n")
    cat("  First 20:", substr(seq, 1, 20), "\n")
    cat("  Last 20: ", substr(seq, nchar(seq) - 19, nchar(seq)), "\n")
  }
}

# Check if it's trailing gaps/characters
max_len <- max(lens)
cat("\n--- Checking trailing characters ---\n")
for (l in sort(unique(lens))) {
  if (l < max_len) {
    # Get a max-length seq and a short seq
    short_seq <- tip_seqs_list[lens == l][1]
    long_seq <- tip_seqs_list[lens == max_len][1]
    
    cat("\nComparing length", l, "vs", max_len, "\n")
    cat("Short ends with:", substr(short_seq, nchar(short_seq) - 29, nchar(short_seq)), "\n")
    cat("Long ends with: ", substr(long_seq, nchar(long_seq) - 29, nchar(long_seq)), "\n")
    cat("Long extra tail:", substr(long_seq, l + 1, max_len), "\n")
  }
}

cat("ASR sequence length:", nchar(asr_states$sequence[1]), "\n")

max_len <- max(nchar(tip_seqs_list))
tip_seqs_list <- sapply(tip_seqs_list, function(s) {
  if (nchar(s) < max_len) {
    paste0(paste(rep("-", max_len - nchar(s)), collapse = ""), s)  # prepend gaps
  } else s
}, USE.NAMES = TRUE)

stopifnot(length(unique(nchar(tip_seqs_list))) == 1)


# Internal node sequences from RAxML ASR (Node1, Node2, ...)
asr_lines <- readLines(states_file)
asr_data <- strsplit(asr_lines, "\t")
asr_states <- data.table(
  node = sapply(asr_data, `[`, 1),
  sequence = sapply(asr_data, `[`, 2)
)
stopifnot(all(grepl("^Node", asr_states$node)))
stopifnot(nrow(asr_states) == n_internal)

asr_node_nums <- as.integer(gsub("Node", "", asr_states$node))
expected_nodes <- 1:n_internal
missing <- setdiff(expected_nodes, asr_node_nums)
cat("Missing node(s):", missing, "\n")
cat("This corresponds to ape node:", n_tips + missing, "\n")
cat("Root node is:", n_tips + 1, "\n")

problem_node <- n_tips + 6093
cat("Node", problem_node, "appears in edges:\n")
print(tree$edge[tree$edge[,1] == problem_node | tree$edge[,2] == problem_node, ])

# Is it a parent or child?
is_parent <- any(tree$edge[,1] == problem_node)
is_child <- any(tree$edge[,2] == problem_node)
cat("Is parent:", is_parent, "Is child:", is_child, "\n")

# What are its descendants?
desc <- Descendants(tree, problem_node, type = "tips")[[1]]
cat("N descendant tips:", length(desc), "\n")
cat("Descendant tips:", head(tree$tip.label[desc]), "\n")
found_ids <- tree$tip.label[desc]
data$Organism[data$ID%in%found_ids]

# Identify outgroup tips and all edges within/leading to outgroup
outgroup_tips <- which(tree$tip.label %in% pinales_in_tree)
outgroup_mrca <- getMRCA(tree, outgroup_tips)
outgroup_nodes <- c(outgroup_tips, Descendants(tree, outgroup_mrca, type = "all")[[1]])
outgroup_nodes <- c(outgroup_nodes, outgroup_mrca)

# Edges to exclude: any edge where child is in outgroup
edges_to_exclude <- which(tree$edge[,2] %in% outgroup_nodes)
cat("Excluding", length(edges_to_exclude), "outgroup edges\n")

# Relax ASR check
stopifnot(nrow(asr_states) >= n_internal - 1)

# Build seq_lookup (missing node stays "")
seq_lookup <- character(n_tips + n_internal)
seq_lookup[1:n_tips] <- tip_seqs_list[tree$tip.label]
internal_nums <- as.integer(gsub("Node", "", asr_states$node))
seq_lookup[n_tips + internal_nums] <- asr_states$sequence

# Pass edges_to_exclude to extraction function

# Build unified seq_lookup indexed by ape node numbers
# ape: tips = 1:n_tips, internals = (n_tips+1):(n_tips+n_internal)
seq_lookup <- character(n_tips + n_internal)
seq_lookup[1:n_tips] <- tip_seqs_list[tree$tip.label]

internal_nums <- as.integer(gsub("Node", "", asr_states$node))
seq_lookup[n_tips + internal_nums] <- asr_states$sequence

stopifnot(sum(seq_lookup == "") == 1)


stopifnot(nchar(seq_lookup[1]) == nchar(seq_lookup[n_tips + 1]))  # tips and internals same length
cat("Sequences loaded. Alignment length:", nchar(seq_lookup[1]), "\n")

# ==== EXTRACT SUBSTITUTIONS ====

extract_substitutions <- function(tree, seq_lookup, exclude_edges = integer(0)) {
  edges <- tree$edge
  n_edges <- nrow(edges)
  subs_list <- vector("list", n_edges)
  
  for (i in seq_len(n_edges)) {
    if (i %in% exclude_edges) next
    if (seq_lookup[edges[i,1]] == "" || seq_lookup[edges[i,2]] == "") next

    parent <- edges[i, 1]
    child <- edges[i, 2]
    
    parent_seq <- strsplit(seq_lookup[parent], "")[[1]]
    child_seq <- strsplit(seq_lookup[child], "")[[1]]
    stopifnot(length(parent_seq) == length(child_seq))
    
    diff_pos <- which(parent_seq != child_seq)
    
    if (length(diff_pos) > 0) {
      subs_list[[i]] <- data.table(
        edge_idx = i,
        parent_node = parent,
        child_node = child,
        position = diff_pos,
        from_aa = parent_seq[diff_pos],
        to_aa = child_seq[diff_pos],
        branch_length = tree$edge.length[i]
      )
    }
  }
  rbindlist(subs_list)
}

cat("Extracting substitutions...\n")
subs <- extract_substitutions(tree, seq_lookup, edges_to_exclude)
cat("Found", nrow(subs), "substitutions across", length(unique(subs$edge_idx)), "edges\n")

# ==== PHENOTYPE FOR DIRECTIONALITY ====


library(phytools)

# tip phenotype
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID
tip_pheno <- pheno[tree$tip.label]

# ML ancestral reconstruction (BM)
fit <- fastAnc(tree, tip_pheno, vars = TRUE, CI = FALSE)

node_pheno <- c(tip_pheno, fit$ace)
names(node_pheno) <- c(tree$tip.label,
                       as.character((length(tree$tip.label)+1):
                                      (length(tree$tip.label)+tree$Nnode)))

## 2) Map edge indices → parent/child phenotypes
edge <- tree$edge
edge_dt <- data.table(
  edge_id = seq_len(nrow(edge)),
  parent_node = edge[,1],
  child_node  = edge[,2]
)

thresh <- sd(tip_pheno, na.rm = TRUE) * 0.5

subs[, parent_pheno := node_pheno[as.character(parent_node)]]
subs[, child_pheno  := node_pheno[as.character(child_node)]]
subs[, pheno_change := child_pheno - parent_pheno]
hist(subs$pheno_change)
subs[, direction := fcase(
  pheno_change >  thresh, "toward_warm",
  pheno_change < -thresh, "toward_cold",
  default = "neutral"
)]


subs[, .N, by = direction]
subs[, .N, by = .(to_aa, direction)]


# ==== LOAD GWAS RESULTS ====
# UPDATE: Uncomment and set model_files path
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)

# Load partition map
partition_map <- readRDS("raxml_input/partitionMap.rds")

gwas_results <- rbindlist(lapply(model_files, function(f) {
  cat(f)
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      GenePos = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs,
      P_aa_only = m$P_aa_only,
      N = m$N,
      R2_full = m$R2_full,
      effects = list(m$effects),
      residue_counts = list(m$residue_counts)
    )
  }))
}))

summary(gwas_results)
# Add sig_class
gwas_thresh_ctrl <- quantile(gwas_results$P_aa_with_pcs, 0.05, na.rm = TRUE)
gwas_thresh_noctrl <- quantile(gwas_results$P_aa_only, 0.20, na.rm = TRUE)

gwas_results[, sig_class := fcase(
  P_aa_with_pcs < gwas_thresh_ctrl & P_aa_only < gwas_thresh_noctrl, "sig_both",
  P_aa_with_pcs < gwas_thresh_ctrl, "sig_control_only",
  P_aa_only < gwas_thresh_noctrl, "sig_nocontrol_only",
  default = "not_sig"
)]

table(gwas_results$sig_class)

cat("Sig class distribution:\n")
print(gwas_results[, .N, by = sig_class])

# Add GlobalPos via partition_map
gwas_results <- merge(
  gwas_results,
  partition_map[, .(Gene, GenePos, GlobalPos)],
  by = c("Gene", "GenePos"),
  all.x = TRUE
)
stopifnot(sum(is.na(gwas_results$GlobalPos)) == 0)

# How many missing?
cat("Missing GlobalPos:", sum(is.na(gwas_results$GlobalPos)), "of", nrow(gwas_results), "\n")

# Which Gene/GenePos combos aren't in partition_map?
missing <- gwas_results[is.na(GlobalPos), .(Gene, GenePos)]
cat("\nMissing by gene:\n")
print(missing[, .N, by = Gene])

# Check if it's a GenePos range issue
cat("\nExample missing positions:\n")
print(head(missing, 10))

# Compare ranges
for (g in unique(missing$Gene)[1:3]) {
  gwas_range <- range(gwas_results[Gene == g, GenePos])
  pm_range <- range(partition_map[Gene == g, GenePos])
  cat(g, "- GWAS:", gwas_range[1], "-", gwas_range[2], 
      "| partition_map:", pm_range[1], "-", pm_range[2], "\n")
}

# Remove GWAS positions not in supermatrix
gwas_results <- gwas_results[!is.na(GlobalPos)]
cat("GWAS positions retained:", nrow(gwas_results), "\n")

# Now merge with substitutions
subs_gwas <- merge(
  subs,
  gwas_results[, .(position = GlobalPos, Gene, GenePos, P_aa_with_pcs, P_aa_only, 
                   sig_class, effects, residue_counts)],
  by = "position",
  all.x = TRUE
)

# Check: some subs may be at positions not in GWAS (if supermatrix has more positions)
cat("Subs with GWAS match:", sum(!is.na(subs_gwas$sig_class)), "\n")
cat("Subs without GWAS match:", sum(is.na(subs_gwas$sig_class)), "\n")

# For unmatched, assign "not_in_gwas"
subs_gwas[is.na(sig_class), sig_class := "not_in_gwas"]

cat("\nSubstitutions by sig_class:\n")
print(subs_gwas[, .N, by = sig_class])


stopifnot(sum(is.na(subs_gwas$sig_class)) == 0)

cat("\nSubstitutions by sig_class:\n")
print(subs_gwas[, .N, by = sig_class])

# Check merge worked
cat("Substitutions at GWAS hits:", sum(subs_gwas$gwas_hit), "\n")

# Load one gene's GWAS to check
test_gene <- "atpA"
gwas_gene <- readRDS(paste0("results/residue_models_triple/", test_gene, "_effects.rds"))

# Pick a position
test_gwas <- gwas_gene[[5]]
test_global <- partition_map[Gene == test_gene & GenePos == test_gwas$Aligned_Position, GlobalPos]

cat("\n=== Validation: Gene", test_gene, "GenePos", test_gwas$Aligned_Position, 
    "-> GlobalPos", test_global, "===\n")

cat("\nGWAS residue_counts (tip states):\n")
print(unlist(test_gwas$residue_counts))

cat("\nASR 'to_aa' at this position:\n")
print(table(subs[position == test_global, to_aa]))

# Directionality by GWAS hit status
dir_by_gwas <- subs_gwas[, .(
  n_subs = .N,
  n_warm = sum(direction == "toward_warm"),
  n_cold = sum(direction == "toward_cold"),
  prop_warm = mean(direction == "toward_warm")
), by = sig_class]

test <- subs_gwas[sig_class == "sig_both" & !is.na(effects)]
cat("Sig_both subs with effects:", nrow(test), "\n")

barplot(table(subs_gwas$sig_class))

print(dir_by_gwas)

# Test: is directionality biased at GWAS hits vs background?
ct <- subs_gwas[direction != "neutral", table(sig_class, direction)]
print(ct)
print(chisq.test(ct))
# ---- test directionality for significan

# Flatten effects into one row per residue per position
gwas_effects <- rbindlist(lapply(model_files, function(f) {
  cat(f)
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    dt <- copy(m$effects)
    dt[, Gene := m$Gene]
    dt[, GenePos := m$Aligned_Position]
    dt[, residue := gsub("X_aares_factor", "", Residue)]
    dt[, .(Gene, GenePos, residue, Effect, SE, P_value)]
  }))
}))

# Add GlobalPos
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)], 
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]
gwas_effects[, residue := gsub("^X_aa", "", residue)]
setkey(gwas_effects, GlobalPos, residue)

subs_gwas[, c("to_effect", "to_pval") := gwas_effects[.(position, to_aa), .(Effect, P_value), on = c("GlobalPos", "residue")]]
subs_gwas[, c("from_effect", "from_pval") := gwas_effects[.(position, from_aa), .(Effect, P_value), on = c("GlobalPos", "residue")]]

cat("Subs with to_effect:", sum(!is.na(subs_gwas$to_effect)), "\n")

cat("Subs with to_p:", sum(!is.na(subs_gwas$to_pval)), "\n")

subs_gwas$variant_sig_class <- subs_gwas$sig_class
subs_gwas$variant_sig_class[subs_gwas$to_pval>0.05] <- "not_sig"
barplot(table(subs_gwas$variant_sig_class))

# === 1. Basic: Effect sign vs phenotype direction ===
test <- subs_gwas[!is.na(to_effect)]
test$effect_change <- test$from_effect - test$to_effect

test[, effect_sign := fcase(
  effect_change >  thresh,  "hot",
  effect_change < -thresh,  "cold",
  default = "neutral"
)]

table(test$effect_sign)

test[, pheno_sign := fcase(
  pheno_change >  thresh,  "hot",
  pheno_change < -thresh,  "cold",
  default = "neutral"
)]

hist(test$pheno_change)
table(test$pheno_sign)

cat("=== Effect sign vs phenotype direction ===\n")
ct1 <- table(effect_sign = test$effect_sign, pheno_sign = test$pheno_sign)
print(ct1)
print(chisq.test(ct1))

# Proportion concordant
test[, concordant := effect_sign == pheno_sign]
cat("\nOverall concordance:", mean(test$concordant, na.rm = TRUE), "\n")

# === 2. Concordance by sig_class ===
cat("\n=== Concordance by sig_class ===\n")
conc_by_class <- test[, .(
  n = .N,
  concordance = mean(concordant, na.rm = TRUE)
), by = variant_sig_class]
print(conc_by_class[order(-concordance)])

#inflated by neutral

test_eff <- test[effect_sign %in% c("hot","cold")]

ct_eff <- table(
  effect_sign = test_eff$effect_sign,
  pheno_sign  = test_eff$pheno_sign
)

print(ct_eff)
chisq.test(ct_eff)

mean(test_eff$effect_sign == test_eff$pheno_sign)

test_pheno <- test[pheno_sign %in% c("hot","cold")]

ct_pheno <- table(
  effect_sign = test_pheno$effect_sign,
  pheno_sign  = test_pheno$pheno_sign
)

print(ct_pheno)
chisq.test(ct_pheno)

mean(test_pheno$effect_sign == test_pheno$pheno_sign)

test_both <- test[
  effect_sign %in% c("hot","cold") &
    pheno_sign  %in% c("hot","cold")
]

ct_both <- table(
  effect_sign = test_both$effect_sign,
  pheno_sign  = test_both$pheno_sign
)

print(ct_both)
chisq.test(ct_both)

# Subset to directional only
test_both <- test[
  effect_sign %in% c("hot","cold") &
    pheno_sign  %in% c("hot","cold")
]

# Stratified concordance tables
strat_tables <- test_both[, .(
  ct = list(table(effect_sign = effect_sign, pheno_sign = pheno_sign)),
  N  = .N
), by = variant_sig_class]

print(strat_tables)
# Add χ² test and proportion concordant
strat_tables[, `:=`(
  chi2 = lapply(ct, chisq.test),
  concordance = sapply(ct, function(x) sum(diag(x)) / sum(x))
)]

# Print
for (i in seq_len(nrow(strat_tables))) {
  cat("\n=== variant_sig_class:", strat_tables$variant_sig_class[i], "===\n")
  print(strat_tables$ct[[i]])
  print(strat_tables$chi2[[i]])
  cat("Proportion concordant:", strat_tables$concordance[i], "\n")
}

#---- done testing concordance ----

# Mean phenotype of ALL tip descendants vs parent phenotype
get_subtree_tip_pheno <- function(child_node, tip_pheno, tree) {
  desc_tips <- Descendants(tree, child_node, type = "tips")[[1]]
  if (length(desc_tips) == 0) {
    # child_node is itself a tip
    return(tip_pheno[tree$tip.label[child_node]])
  }
  mean(tip_pheno[tree$tip.label[desc_tips]], na.rm = TRUE)
}

# Compute for unique child nodes
unique_children <- unique(subs_gwas$child_node)
cat("Computing subtree tip phenotypes for", length(unique_children), "nodes...\n")

subtree_pheno <- sapply(unique_children, function(cn) {
  get_subtree_tip_pheno(cn, tip_pheno, tree)
})
names(subtree_pheno) <- unique_children

subs_gwas[, subtree_tip_pheno := subtree_pheno[as.character(child_node)]]
subs_gwas[, subtree_shift := subtree_tip_pheno - parent_pheno]

# Now: does subtree_shift correlate with effect_change?
test <- subs_gwas[!is.na(to_effect) & !is.na(from_effect) & !is.na(subtree_shift)]
test[, effect_change := to_effect - from_effect]

strat_cor <- test[, {
  if (.N < 10) {
    list(n = .N, cor = NA_real_, pval = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
  } else {
    ct <- cor.test(effect_change, subtree_shift)
    list(
      n = .N,
      cor = ct$estimate,
      pval = ct$p.value,
      ci_lo = ct$conf.int[1],
      ci_hi = ct$conf.int[2]
    )
  }
}, by = variant_sig_class]

print(strat_cor[order(-cor)])


# Fisher z-transform to compare correlations
fisher_z_test <- function(r1, n1, r2, n2) {
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z_diff <- (z1 - z2) / se
  pval <- 2 * pnorm(-abs(z_diff))
  list(z = z_diff, p = pval)
}

# sig_both vs not_sig
comp <- fisher_z_test(
  strat_cor[variant_sig_class == "sig_both", cor],
  strat_cor[variant_sig_class == "sig_both", n],
  strat_cor[variant_sig_class == "not_sig", cor],
  strat_cor[variant_sig_class == "not_sig", n]
)
cat("sig_both vs not_sig: z =", round(comp$z, 3), ", p =", format(comp$p, digits = 3), "\n")

# Option 2: 2D hex density, faceted
ggplot(test, aes(x = effect_change, y = subtree_shift)) +
  geom_hex(bins = 50) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  facet_wrap(~variant_sig_class, scales = "free") +
  scale_fill_viridis_c(trans = "log10") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Effect change", y = "Subtree shift") +
  theme_minimal()

test[, effect_bin := cut(effect_change, breaks = quantile(effect_change, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]

binned <- test[, .(
  mean_effect = mean(effect_change, na.rm = TRUE),
  mean_shift = mean(subtree_shift, na.rm = TRUE),
  se_shift = sd(subtree_shift, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(variant_sig_class, effect_bin)]

ggplot(binned, aes(x = mean_effect, y = mean_shift, color = variant_sig_class)) +
  geom_pointrange(aes(ymin = mean_shift - 1.96*se_shift, ymax = mean_shift + 1.96*se_shift)) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "GWAS effect change (to - from)", 
       y = "Subtree phenotype shift",
       title = "Binned effect vs downstream phenotype") +
  theme_minimal()


test[, effect_dir := sign(effect_change)]
test[, shift_dir := sign(subtree_shift)]

# Proportion concordant by class (excluding zeros)
test_nonzero <- test[effect_dir != 0 & shift_dir != 0]

concordance <- test_nonzero[, .(
  n = .N,
  prop_concordant = mean(effect_dir == shift_dir),
  se = sqrt(mean(effect_dir == shift_dir) * (1 - mean(effect_dir == shift_dir)) / .N)
), by = variant_sig_class]

concordance[, `:=`(
  ci_lo = prop_concordant - 1.96 * se,
  ci_hi = prop_concordant + 1.96 * se
)]

print(concordance[order(-prop_concordant)])

# Test vs 0.5 (null expectation)
concordance[, pval := binom.test(round(prop_concordant * n), n, p = 0.5)$p.value, by = variant_sig_class]


# Precompute for unique child_node + effect_sign combinations
test_dir <- test[effect_change != 0]
test_dir[, effect_sign := sign(effect_change)]

unique_combos <- unique(test_dir[, .(child_node, parent_pheno, effect_sign)])
cat("Computing proportion concordant for", nrow(unique_combos), "unique combos...\n")

prop_results <- unique_combos[, {
  desc_tips <- Descendants(tree, child_node, type = "tips")[[1]]
  if (length(desc_tips) == 0) desc_tips <- child_node
  
  tip_vals <- tip_pheno[tree$tip.label[desc_tips]]
  shifts <- tip_vals - parent_pheno
  
  prop <- if (effect_sign > 0) {
    mean(shifts > 0, na.rm = TRUE)
  } else {
    mean(shifts < 0, na.rm = TRUE)
  }
  
  .(prop_concordant = prop, n_tips = length(desc_tips))
}, by = .(child_node, parent_pheno, effect_sign)]

# Merge back
test_dir <- merge(test_dir, prop_results, by = c("child_node", "parent_pheno", "effect_sign"), all.x = TRUE)

# Summarize by sig class
prop_summary <- test_dir[, .(
  n_subs = .N,
  mean_prop = mean(prop_concordant, na.rm = TRUE),
  median_prop = median(prop_concordant, na.rm = TRUE),
  se = sd(prop_concordant, na.rm = TRUE) / sqrt(.N)
), by = variant_sig_class]

prop_summary[, `:=`(
  ci_lo = mean_prop - 1.96 * se,
  ci_hi = mean_prop + 1.96 * se
)]

print(prop_summary[order(-mean_prop)])

# Test if mean_prop > 0.5
prop_summary[, tstat := (mean_prop - 0.5) / se]
prop_summary[, pval := 2 * pt(-abs(tstat), df = n_subs - 1)]
print(prop_summary[order(-mean_prop), .(variant_sig_class, n_subs, mean_prop, median_prop, pval)])


test_dir[, depth := node.depth(tree, method = 1)[child_node]]

depth_summary <- test_dir[variant_sig_class == "sig_both", .(
  n = .N,
  mean_prop = mean(prop_concordant, na.rm = TRUE)
), by = cut(depth, breaks = 10)]

print(depth_summary)
test_dir[, effect_magnitude := abs(effect_change)]
test_dir[variant_sig_class == "sig_both", 
         cor.test(effect_magnitude, prop_concordant)]

# ==== VISUALIZATIONS ====

library(ggplot2)

# 1. Correlation bar plot with CIs
cor_plot_data <- data.table(
  sig_class = factor(c("sig_both", "sig_control_only", "sig_nocontrol_only", "not_sig"),
                     levels = c("sig_both", "sig_control_only", "sig_nocontrol_only", "not_sig")),
  cor = c(0.0141, 0.0096, 0.0053, 0.0032),
  ci_lo = c(0.0091, 0.0025, 0.0037, 0.0018),
  ci_hi = c(0.0190, 0.0168, 0.0070, 0.0047)
)

p1 <- ggplot(cor_plot_data, aes(x = sig_class, y = cor)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "GWAS significance class", 
       y = "Correlation (effect change vs subtree shift)",
       title = "GWAS effect-phenotype correlation by significance class") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p1
# 2. Concordance plot (proportion above 0.5)
conc_plot_data <- data.table(
  sig_class = factor(c("sig_both", "sig_control_only", "not_sig", "sig_nocontrol_only"),
                     levels = c("sig_both", "sig_control_only", "not_sig", "sig_nocontrol_only")),
  prop = c(0.5083, 0.5062, 0.5048, 0.5044),
  ci_lo = c(0.5058, 0.5026, 0.5040, 0.5036),
  ci_hi = c(0.5108, 0.5098, 0.5055, 0.5052)
)

p2 <- ggplot(conc_plot_data, aes(x = sig_class, y = prop)) +
  geom_col(fill = "darkorange", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.495, 0.515), labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "GWAS significance class",
       y = "Proportion concordant",
       title = "Sign concordance: subtree shift vs GWAS effect") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p2
# 3. Combined panel
library(patchwork)
p_combined <- p1 + p2 + plot_annotation(tag_levels = 'A')

ggsave("gwas_phylo_validation.pdf", p_combined, width = 12, height = 5)
ggsave("gwas_phylo_validation.png", p_combined, width = 12, height = 5, dpi = 150)

branch_summary <- test[, .(
  total_effect_change = sum(effect_change, na.rm = TRUE),
  n_subs = .N,
  subtree_shift = subtree_shift[1],  # same for all subs on branch
  parent_pheno = parent_pheno[1]
), by = .(edge_idx)]

# ---- branch-wise tests ----
# Now N = number of branches with substitutions
cor.test(branch_summary$total_effect_change, branch_summary$subtree_shift)

plot(branch_summary$total_effect_change, branch_summary$subtree_shift)


set.seed(42)
n_perm <- 1000

# Get unique site effects
site_effects <- unique(test[, .(position, to_aa, to_effect)])

observed_cor <- cor(test$effect_change, test$subtree_shift, use = "complete.obs")

perm_cors <- replicate(n_perm, {
  # Shuffle effect assignments across sites
  shuffled <- site_effects[sample(.N)]
  shuffled[, shuffled_effect := site_effects$to_effect]
  
  test_perm <- merge(test[, .(position, to_aa, subtree_shift, from_effect)], 
                     shuffled[, .(position, to_aa, shuffled_effect)],
                     by = c("position", "to_aa"))
  test_perm[, effect_change_perm := shuffled_effect - from_effect]
  
  cor(test_perm$effect_change_perm, test_perm$subtree_shift, use = "complete.obs")
})

perm_pval <- mean(abs(perm_cors) >= abs(observed_cor))

hist(c(perm_cors, observed_cor))
abline(v=observed_cor, col="seagreen")

# ---- upstream test ---- 

# Get phenotype trend over N ancestors
get_ancestral_trend <- function(node, tree, node_pheno, n_ancestors = 5) {
  edge_dt <- data.table(child = tree$edge[,2], parent = tree$edge[,1])
  setkey(edge_dt, child)
  
  current <- node
  phenos <- numeric(n_ancestors + 1)
  phenos[1] <- node_pheno[as.character(current)]
  
  for (i in seq_len(n_ancestors)) {
    parent <- edge_dt[.(current), parent]
    if (is.na(parent) || length(parent) == 0) {
      phenos[(i+1):(n_ancestors+1)] <- NA
      break
    }
    phenos[i + 1] <- node_pheno[as.character(parent)]
    current <- parent
  }
  
  # Return slope (trend) if enough data
  valid <- which(!is.na(phenos))
  if (length(valid) >= 3) {
    # Positions: 0 = focal node, 1 = parent, 2 = grandparent, etc.
    # Negative slope = ancestors were colder (warming trend toward focal)
    fit <- lm(phenos[valid] ~ valid)
    return(-coef(fit)[2])  # flip sign so positive = warming trend
  }
  NA_real_
}

# Compute for unique parent nodes
unique_parents <- unique(test$parent_node)
cat("Computing ancestral trends for", length(unique_parents), "nodes...\n")

ancestral_trends <- sapply(unique_parents, function(pn) {
  get_ancestral_trend(pn, tree, node_pheno, n_ancestors = 5)
})
names(ancestral_trends) <- unique_parents

test[, preceding_trend := ancestral_trends[as.character(parent_node)]]

# Test
test_trend <- test[!is.na(preceding_trend) & !is.na(effect_change)]
cor.test(test_trend$preceding_trend, test_trend$effect_change)

# Stratified
strat_cor_trend <- test_trend[, {
  if (.N < 10) {
    list(n = .N, cor = NA_real_, pval = NA_real_)
  } else {
    ct <- cor.test(preceding_trend, effect_change)
    list(n = .N, cor = ct$estimate, pval = ct$p.value)
  }
}, by = variant_sig_class]

print(strat_cor_trend[order(-cor)])

comparison <- data.table(
  model = c("consequence", "precedent"),
  cor_sig_both = c(
    strat_cor[variant_sig_class == "sig_both", cor],      # from earlier
    strat_cor_trend[variant_sig_class == "sig_both", cor]
  ),
  pval_sig_both = c(
    strat_cor[variant_sig_class == "sig_both", pval],
    strat_cor_trend[variant_sig_class == "sig_both", pval]
  )
)

print(comparison)
# ---- clean analysis ---- 

# 1. Load ASR probabilities from RAxML
# RAxML outputs marginal probabilities - need to parse them
# Check if there's a .ancestralProbs file or similar
probs_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralProbs"
stopifnot(file.exists(probs_file))

# Parse the probs file - format depends on RAxML output
# Typically: Node\tSite\tAA_probs...
probs_lines <- readLines(probs_file, n = 5)
cat("First 5 lines of probs file:\n")
print(probs_lines)

# Adjust parsing based on actual format - placeholder:
asr_probs <- fread(probs_file)

summary(asr_probs)

# 2. For each substitution, get the ML probability of both parent and child states
# We want: P(parent_aa at parent_node) and P(child_aa at child_node)
# Only keep substitutions where both probabilities > threshold (e.g., 0.8)

ASR_PROB_THRESH <- 0.9

setDT(asr_probs)

# 1. Get max probability (ML state confidence) for each node/site
aa_cols <- paste0("p_", c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"))
asr_probs[, max_prob := do.call(pmax, .SD), .SDcols = aa_cols]

# 2. Convert to long format for easy lookup
asr_probs_long <- melt(asr_probs, 
                       id.vars = c("Node", "Site", "State", "max_prob"),
                       measure.vars = aa_cols,
                       variable.name = "aa_col", 
                       value.name = "prob")
asr_probs_long[, aa := gsub("p_", "", aa_col)]
asr_probs_long <- asr_probs_long[, .(Node, Site, aa, prob, max_prob)]

# Convert Node to numeric (matches ape internal node numbering)
asr_probs_long[, node_num := as.integer(gsub("Node", "", Node)) + n_tips]
setkey(asr_probs_long, node_num, Site, aa)

# 3. For each substitution, get probability of parent and child states
subs_gwas[, parent_state_prob := asr_probs_long[.(parent_node, position, from_aa), prob]]
subs_gwas[, child_state_prob := asr_probs_long[.(child_node, position, to_aa), prob]]

cat("Substitutions with parent_state_prob:", sum(!is.na(subs_gwas$parent_state_prob)), "\n")
cat("Substitutions with child_state_prob:", sum(!is.na(subs_gwas$child_state_prob)), "\n")

# Check distribution of state probabilities
cat("\nParent state prob distribution:\n")
print(summary(subs_gwas$parent_state_prob))
hist(subs_gwas$parent_state_prob)
cat("\nChild state prob distribution:\n")
print(summary(subs_gwas$child_state_prob))
hist(subs_gwas$child_state_prob)
# 4. Identify back-mutations
identify_back_mutations_fast <- function(subs_dt, tree) {
  cat("Building ancestry matrix...\n")
  
  n_nodes <- n_tips + n_internal
  anc_list <- Ancestors(tree, 1:n_nodes, type = "all")
  
  cat("Processing positions...\n")
  
  positions <- unique(subs_dt$position)
  backmut_edges <- vector("list", length(positions))
  
  pb <- txtProgressBar(min = 0, max = length(positions), style = 3)
  
  for (i in seq_along(positions)) {
    pos <- positions[i]
    grp <- subs_dt[position == pos, .(edge_idx, child_node)]
    children <- grp$child_node
    edges <- grp$edge_idx
    n <- length(children)
    
    # For each child, get its ancestors
    # A child has backmut if it's an ancestor of ANY other child in the group
    has_bm <- logical(n)
    
    # Collect all ancestors of all children
    all_ancestors <- unique(unlist(anc_list[children]))
    
    # A child has backmut if it appears in another child's ancestor list
    for (j in seq_len(n)) {
      # Is child_j an ancestor of any other child?
      other_children <- children[-j]
      other_ancestors <- unique(unlist(anc_list[other_children]))
      if (children[j] %in% other_ancestors) {
        has_bm[j] <- TRUE
      }
    }
    
    backmut_edges[[i]] <- edges[has_bm]
    
    if (i %% 500 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  
  backmut_edges <- unique(unlist(backmut_edges))
  cat("\nEdges with back-mutations:", length(backmut_edges), "\n")
  
  subs_dt[, has_backmut := edge_idx %in% backmut_edges]
  subs_dt
}

system.time({
  subs_gwas <- identify_back_mutations_fast(subs_gwas, tree)
})
summary(subs_gwas)
# Check subs per position
subs_per_pos <- subs_gwas[, .N, by = .(position, sig_class)][order(-N)]
cat("Top positions by sub count:\n")
print(head(subs_per_pos, 20))
cat("\nDistribution:\n")
print(summary(subs_per_pos$N))
cat("\nPositions with >100 subs:", sum(subs_per_pos$N > 100), "\n")
cat("Positions with >10 subs:", sum(subs_per_pos$N > 10), "\n")
hist(subs_per_pos$N)
boxplot(N ~ sig_class, subs_per_pos)

cat("Identifying back-mutations...\n")
#subs_gwas <- identify_back_mutations(subs_gwas, tree)
cat("Substitutions with back-mutation downstream:", sum(subs_gwas$has_backmut), 
    "of", nrow(subs_gwas), "(", round(100*mean(subs_gwas$has_backmut), 1), "%)\n")
boxplot(has_backmut ~ sig_class, subs_gwas)

backmut_summary <- subs_gwas[, .(
  total_backmuts = sum(has_backmut, na.rm = TRUE),
  n_muts = .N
), by = .(sig_class)]


backmut_summary$prop_backmut <- backmut_summary$total_backmuts / backmut_summary$n_muts
backmut_summary
barplot(backmut_summary$prop_backmut)
# 5. Define clean dataset
ASR_PROB_THRESH <- 0.8

subs_gwas[, high_conf_parent := parent_state_prob >= ASR_PROB_THRESH]
subs_gwas[, high_conf_child := child_state_prob >= ASR_PROB_THRESH]
subs_gwas[, high_conf_both := high_conf_parent & high_conf_child]

cat("\nFiltering summary:\n")
cat("High-conf parent (>=", ASR_PROB_THRESH, "):", sum(subs_gwas$high_conf_parent, na.rm=TRUE), "\n")
cat("High-conf child (>=", ASR_PROB_THRESH, "):", sum(subs_gwas$high_conf_child, na.rm=TRUE), "\n")
cat("High-conf both:", sum(subs_gwas$high_conf_both, na.rm=TRUE), "\n")
cat("No back-mutation:", sum(!subs_gwas$has_backmut), "\n")
cat("Clean (both criteria):", sum(subs_gwas$high_conf_both & !subs_gwas$has_backmut, na.rm=TRUE), "\n")

# 6. Run analysis on progressively cleaner subsets
run_analysis <- function(subs_subset, label) {
  test_sub <- subs_subset[!is.na(to_effect) & !is.na(from_effect) & !is.na(subtree_shift)]
  test_sub[, effect_change := to_effect - from_effect]
  
  if (nrow(test_sub) < 100) {
    return(data.table(dataset = label, n = nrow(test_sub), cor = NA_real_, pval = NA_real_))
  }
  
  # Overall correlation
  ct <- cor.test(test_sub$effect_change, test_sub$subtree_shift)
  
  # Stratified by sig_class
  strat <- test_sub[, {
    if (.N < 10) list(cor = NA_real_, pval = NA_real_)
    else {
      c <- cor.test(effect_change, subtree_shift)
      list(cor = c$estimate, pval = c$p.value)
    }
  }, by = variant_sig_class]
  
  sig_both_cor <- strat[variant_sig_class == "sig_both", cor]
  sig_both_pval <- strat[variant_sig_class == "sig_both", pval]
  
  data.table(
    dataset = label,
    n = nrow(test_sub),
    n_sig_both = test_sub[variant_sig_class == "sig_both", .N],
    cor_overall = ct$estimate,
    pval_overall = ct$p.value,
    cor_sig_both = sig_both_cor,
    pval_sig_both = sig_both_pval
  )
}

results_comparison <- rbindlist(list(
  run_analysis(subs_gwas, "full"),
  run_analysis(subs_gwas[has_backmut == FALSE], "no_backmut"),
  run_analysis(subs_gwas[high_conf_both == TRUE], "high_conf_ASR"),
  run_analysis(subs_gwas[high_conf_both == TRUE & has_backmut == FALSE], "clean_both")
))

cat("\n=== Comparison across filtering levels ===\n")
print(results_comparison)

# 7. Permutation test on cleanest dataset
subs_clean <- subs_gwas[high_conf_both == TRUE & has_backmut == FALSE]
test_clean <- subs_clean[!is.na(to_effect) & !is.na(from_effect) & !is.na(subtree_shift)]
test_clean[, effect_change := to_effect - from_effect]

if (nrow(test_clean) >= 100) {
  cat("\n=== Permutation test on clean data ===\n")
  
  site_effects_clean <- unique(test_clean[, .(position, to_aa, to_effect)])
  observed_cor_clean <- cor(test_clean$effect_change, test_clean$subtree_shift, use = "complete.obs")
  
  perm_cors_clean <- replicate(1000, {
    shuffled <- site_effects_clean[sample(.N)]
    shuffled[, shuffled_effect := site_effects_clean$to_effect]
    
    test_perm <- merge(test_clean[, .(position, to_aa, subtree_shift, from_effect)],
                       shuffled[, .(position, to_aa, shuffled_effect)],
                       by = c("position", "to_aa"))
    test_perm[, effect_change_perm := shuffled_effect - from_effect]
    
    cor(test_perm$effect_change_perm, test_perm$subtree_shift, use = "complete.obs")
  })
  
  perm_pval_clean <- mean(abs(perm_cors_clean) >= abs(observed_cor_clean))
  
  cat("Observed cor (clean):", round(observed_cor_clean, 5), "\n")
  cat("Permutation null: [", round(min(perm_cors_clean), 5), ",", round(max(perm_cors_clean), 5), "]\n")
  cat("Permutation p:", perm_pval_clean, "\n")
  
  # Histogram
  hist(perm_cors_clean, breaks = 30, main = "Permutation null (clean data)",
       xlab = "Correlation", xlim = range(c(perm_cors_clean, observed_cor_clean)))
  abline(v = observed_cor_clean, col = "seagreen", lwd = 2)
}

# After parsing probs into a usable format (position, node, aa, prob):
# asr_probs_dt <- ...  # data.table with columns: node, position, aa, prob

# Example join logic (adjust once format is known):
# subs_gwas[, parent_state_prob := asr_probs_dt[.(parent_node, position, from_aa), prob]]
# subs_gwas[, child_state_prob := asr_probs_dt[.(child_node, position, to_aa), prob]]
# subs_clean <- subs_gwas[parent_state_prob >= ASR_PROB_THRESH & child_state_prob >= ASR_PROB_THRESH]

# 3. Identify back-mutations: does the site change again in any descendant edge?
identify_back_mutations <- function(subs_dt, tree) {
  # For each substitution, check if any descendant edge also has a sub at same position
  subs_dt[, has_backmut := FALSE]
  
  # Build descendant lookup
  desc_edges <- lapply(1:nrow(tree$edge), function(i) {
    child <- tree$edge[i, 2]
    desc_nodes <- Descendants(tree, child, type = "all")[[1]]
    which(tree$edge[, 1] %in% desc_nodes | tree$edge[, 2] %in% desc_nodes)
  })
  
  # For each position, find all edges with subs at that position
  setkey(subs_dt, position)
  positions_with_subs <- subs_dt[, .(edges = list(unique(edge_idx))), by = position]
  
  # Check each substitution
  for (i in seq_len(nrow(subs_dt))) {
    pos <- subs_dt$position[i]
    edge_i <- subs_dt$edge_idx[i]
    
    # Edges that are descendants of this edge
    desc_edge_idx <- desc_edges[[edge_i]]
    
    # Other edges with subs at this position
    other_edges <- setdiff(positions_with_subs[.(pos), edges[[1]]], edge_i)
    
    # Does any descendant edge have a sub at this position?
    if (any(other_edges %in% desc_edge_idx)) {
      subs_dt[i, has_backmut := TRUE]
    }
  }
  
  subs_dt
}

cat("Identifying back-mutations...\n")
subs_gwas <- identify_back_mutations(subs_gwas, tree)
cat("Substitutions with back-mutation downstream:", sum(subs_gwas$has_backmut), 
    "of", nrow(subs_gwas), "(", round(100*mean(subs_gwas$has_backmut), 1), "%)\n")

# 4. Clean dataset: high-confidence ASR, no back-mutations
# subs_clean <- subs_gwas[parent_state_prob >= ASR_PROB_THRESH & 
#                         child_state_prob >= ASR_PROB_THRESH & 
#                         !has_backmut]

# For now, just filter on back-mutations until we parse ASR probs:
subs_clean <- subs_gwas[has_backmut == FALSE]
cat("Clean substitutions (no back-mut):", nrow(subs_clean), "\n")

# 5. Re-run the main analysis on clean data
test_clean <- subs_clean[!is.na(to_effect) & !is.na(from_effect) & !is.na(subtree_shift)]
test_clean[, effect_change := to_effect - from_effect]

cat("\n=== CLEAN DATA: Correlation by sig_class ===\n")
strat_cor_clean <- test_clean[, {
  if (.N < 10) {
    list(n = .N, cor = NA_real_, pval = NA_real_)
  } else {
    ct <- cor.test(effect_change, subtree_shift)
    list(n = .N, cor = ct$estimate, pval = ct$p.value)
  }
}, by = variant_sig_class]

print(strat_cor_clean[order(-cor)])

# 6. Compare clean vs full
comparison_clean <- data.table(
  dataset = c("full", "no_backmut"),
  n = c(nrow(test), nrow(test_clean)),
  cor_sig_both = c(
    strat_cor[variant_sig_class == "sig_both", cor],
    strat_cor_clean[variant_sig_class == "sig_both", cor]
  ),
  pval_sig_both = c(
    strat_cor[variant_sig_class == "sig_both", pval],
    strat_cor_clean[variant_sig_class == "sig_both", pval]
  )
)

print(comparison_clean)

# 7. Permutation test on clean data
cat("\n=== Permutation test on clean data ===\n")
site_effects_clean <- unique(test_clean[, .(position, to_aa, to_effect)])
observed_cor_clean <- cor(test_clean$effect_change, test_clean$subtree_shift, use = "complete.obs")

perm_cors_clean <- replicate(1000, {
  shuffled <- site_effects_clean[sample(.N)]
  shuffled[, shuffled_effect := site_effects_clean$to_effect]
  
  test_perm <- merge(test_clean[, .(position, to_aa, subtree_shift, from_effect)],
                     shuffled[, .(position, to_aa, shuffled_effect)],
                     by = c("position", "to_aa"))
  test_perm[, effect_change_perm := shuffled_effect - from_effect]
  
  cor(test_perm$effect_change_perm, test_perm$subtree_shift, use = "complete.obs")
})

perm_pval_clean <- mean(abs(perm_cors_clean) >= abs(observed_cor_clean))

cat("Observed cor (clean):", observed_cor_clean, "\n")
cat("Permutation null range:", range(perm_cors_clean), "\n")
cat("Permutation p:", perm_pval_clean, "\n")

# ---- test pheno trajectory ----
# Get all descendants for each node (cached for efficiency)
library(phangorn)
desc_list <- Descendants(tree, 1:(n_tips + n_internal), type = "all")

# Function to get phenotype trajectory after a substitution
get_downstream_pheno <- function(child_node, node_pheno, tree, max_depth = 3) {
  # Get descendants up to max_depth generations
  desc <- Descendants(tree, child_node, type = "children")[[1]]
  
  if (length(desc) == 0) {
    # Terminal - return tip pheno
    return(list(
      immediate = node_pheno[as.character(child_node)],
      descendants = node_pheno[as.character(child_node)],
      n_desc = 1
    ))
  }
  
  # BFS to get nodes within max_depth
  current <- child_node
  all_desc <- c()
  for (d in seq_len(max_depth)) {
    next_gen <- unlist(Descendants(tree, current, type = "children"))
    if (length(next_gen) == 0) break
    all_desc <- c(all_desc, next_gen)
    current <- next_gen
  }
  
  desc_phenos <- node_pheno[as.character(all_desc)]
  
  list(
    immediate = node_pheno[as.character(child_node)],
    descendants = mean(desc_phenos, na.rm = TRUE),
    n_desc = length(all_desc)
  )
}

# For efficiency, vectorize over unique child nodes first
unique_children <- unique(subs_gwas$child_node)
cat("Computing downstream phenotypes for", length(unique_children), "unique child nodes...\n")

downstream_cache <- rbindlist(lapply(unique_children, function(cn) {
  res <- get_downstream_pheno(cn, node_pheno, tree, max_depth = 3)
  data.table(child_node = cn, 
             pheno_immediate = res$immediate,
             pheno_downstream = res$descendants,
             n_descendants = res$n_desc)
}))

# Merge back
subs_gwas <- merge(subs_gwas, downstream_cache, by = "child_node", all.x = TRUE)

# Now compute: does pheno_downstream - parent_pheno align with effect sign?
subs_gwas[, downstream_change := pheno_downstream - parent_pheno]

# Classify
subs_gwas[, downstream_sign := fcase(
  downstream_change >  thresh, "hot",
  downstream_change < -thresh, "cold",
  default = "neutral"
)]

# Test concordance with effect sign
test2 <- subs_gwas[!is.na(to_effect) & !is.na(from_effect)]
test2[, effect_change := to_effect - from_effect]
test2[, effect_sign := fcase(
  effect_change >  thresh, "hot",
  effect_change < -thresh, "cold",
  default = "neutral"
)]

# Filter to directional only
test2_dir <- test2[effect_sign %in% c("hot", "cold") & downstream_sign %in% c("hot", "cold")]

cat("\n=== Downstream concordance (depth=3) ===\n")
ct_downstream <- table(effect = test2_dir$effect_sign, downstream = test2_dir$downstream_sign)
print(ct_downstream)
print(chisq.test(ct_downstream))
cat("Concordance:", sum(diag(ct_downstream)) / sum(ct_downstream), "\n")



mean(test_both$effect_sign == test_both$pheno_sign)

hist(test$pheno_change)
hist(test$to_effect, breaks = 100)

plot(test$effect_change, test$pheno_change, main="Directionality from tip gwas vs node shift",
     xlab="Change in temp effect by substitution", ylab="ASR Temp shift at each node")

# === 3. Restrict to strong effects (low p-value on to_aa) ===
cat("\n=== Concordance for significant residue effects (to_pval < 0.05) ===\n")
test_sig <- test[to_pval < 0.05]
cat("N subs with sig to_effect:", nrow(test_sig), "\n")
cat("Concordance:", mean(test_sig$concordant, na.rm = TRUE), "\n")

conc_sig_by_class <- test_sig[, .(
  n = .N,
  concordance = mean(concordant, na.rm = TRUE)
), by = sig_class]
print(conc_sig_by_class[order(-concordance)])

# === 4. Effect magnitude: do larger effects show stronger concordance? ===
cat("\n=== Concordance by effect magnitude quartile ===\n")
test[, effect_diff := abs(to_effect - from_effect)]
test[, effect_quartile := cut(effect_diff, quantile(effect_diff, 0:4/4, na.rm = TRUE), 
                              labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)]
print(test[, .(n = .N, concordance = mean(concordant, na.rm = TRUE)), by = effect_quartile][order(effect_quartile)])

# === 5. Direction breakdown: warm vs cold separately ===
cat("\n=== Concordance by direction ===\n")
print(test[, .(n = .N, concordance = mean(concordant, na.rm = TRUE)), by = direction])

# === 6. Summary plot data ===
cat("\n=== Summary: sig_both sites with sig effects ===\n")
final_test <- test[sig_class == "sig_both" & to_pval < 0.05]
cat("N:", nrow(final_test), "\n")
cat("Concordance:", mean(final_test$concordant, na.rm = TRUE), "\n")
cat("Expected by chance: 0.5\n")
binom.test(sum(final_test$concordant), nrow(final_test), p = 0.5)

# Variant-level sig_class: demote to not_sig if to_pval >= 0.05
subs_gwas[, variant_sig_class := fcase(
  is.na(to_pval), "no_effect_data",
  to_pval >= 0.05, "not_sig",
  default = sig_class
)]

cat("Variant sig_class distribution:\n")
print(subs_gwas[, .N, by = variant_sig_class][order(-N)])

# Rerun concordance by variant_sig_class
test <- subs_gwas[!is.na(to_effect) & !is.na(from_effect) & direction != "neutral"]
test[, effect_sign := sign(to_effect - from_effect)]
test[, pheno_sign := sign(pheno_change)]
test[, concordant := effect_sign == pheno_sign]

cat("\n=== Concordance by variant_sig_class ===\n")
conc_by_var <- test[, .(
  n = .N,
  concordance = mean(concordant, na.rm = TRUE),
  se = sd(concordant, na.rm = TRUE) / sqrt(.N)
), by = variant_sig_class]
print(conc_by_var[order(-concordance)])

# Binom test for each class vs 0.5
cat("\n=== Binomial tests ===\n")
for (cls in conc_by_var$variant_sig_class) {
  sub <- test[variant_sig_class == cls & !is.na(concordant)]
  bt <- binom.test(sum(sub$concordant), nrow(sub), p = 0.5)
  cat(cls, ": ", round(bt$estimate, 4), " [", round(bt$conf.int[1], 4), "-", 
      round(bt$conf.int[2], 4), "] p=", format.pval(bt$p.value, digits = 3), "\n", sep = "")
}

# Use full subs_gwas, not filtered
test <- subs_gwas[!is.na(to_effect) & !is.na(from_effect)]

# Effect sign with neutral zone (effect_diff < 1 SD)
effect_diff <- test$to_effect - test$from_effect
effect_sd <- sd(abs(effect_diff), na.rm = TRUE)
cat("Effect SD:", effect_sd, "\n")

test[, effect_sign := fcase(
  (to_effect - from_effect) > effect_sd, 1L,
  (to_effect - from_effect) < -effect_sd, -1L,
  default = 0L
)]

# Pheno sign with neutral zone (pheno_change < 1 SD)
pheno_sd <- sd(abs(test$pheno_change), na.rm = TRUE)
cat("Pheno SD:", pheno_sd, "\n")

test[, pheno_sign := fcase(
  pheno_change > pheno_sd, 1L,
  pheno_change < -pheno_sd, -1L,
  default = 0L
)]

cat("\nEffect sign distribution:\n")
print(table(test$effect_sign))

cat("\nPheno sign distribution:\n")
print(table(test$pheno_sign))

# Concordance: only among non-neutral on both axes
test[, concordant := fcase(
  effect_sign == 0L | pheno_sign == 0L, NA,
  effect_sign == pheno_sign, TRUE,
  default = FALSE
)]

cat("\n=== Concordance by variant_sig_class (both non-neutral) ===\n")
conc_by_var <- test[!is.na(concordant), .(
  n = .N,
  concordance = mean(concordant),
  se = sd(concordant) / sqrt(.N)
), by = variant_sig_class]
print(conc_by_var[order(-concordance)])

# Cross-tab: effect_sign x pheno_sign x variant_sig_class
cat("\n=== Full contingency: sig_both ===\n")
print(test[variant_sig_class == "sig_both", table(effect_sign, pheno_sign)])

cat("\n=== Full contingency: not_sig ===\n")
print(test[variant_sig_class == "not_sig", table(effect_sign, pheno_sign)])

# Binomial tests for each class
cat("=== Binomial tests vs 0.5 ===\n")
for (cls in conc_by_var[order(-concordance), variant_sig_class]) {
  sub <- test[variant_sig_class == cls & !is.na(concordant)]
  n_conc <- sum(sub$concordant)
  n_tot <- nrow(sub)
  bt <- binom.test(n_conc, n_tot, p = 0.5)
  cat(sprintf("%s: %.3f [%.3f-%.3f] n=%d p=%s\n", 
              cls, bt$estimate, bt$conf.int[1], bt$conf.int[2], n_tot, format.pval(bt$p.value, digits = 3)))
}

# Chi-square: is sig_both different from not_sig?
cat("\n=== Chi-square: sig_both vs not_sig ===\n")
ct_compare <- test[variant_sig_class %in% c("sig_both", "not_sig") & !is.na(concordant), 
                   table(variant_sig_class, concordant)]
print(ct_compare)
print(chisq.test(ct_compare))

# Odds ratio
cat("\n=== Odds ratio: sig_both vs not_sig ===\n")
or <- (ct_compare["sig_both", "TRUE"] / ct_compare["sig_both", "FALSE"]) / 
  (ct_compare["not_sig", "TRUE"] / ct_compare["not_sig", "FALSE"])
cat("OR:", round(or, 3), "\n")

# Fisher's exact for sig_both alone
cat("\n=== Fisher exact: sig_both contingency ===\n")
ct_sig <- test[variant_sig_class == "sig_both", table(effect_sign, pheno_sign)]
ct_sig_2x2 <- ct_sig[c("-1", "1"), c("-1", "1")]
print(ct_sig_2x2)
print(fisher.test(ct_sig_2x2))


# ==== AUTOCORRELATION ANALYSIS ====
# Test if substitutions of same type cluster on the tree

library(ape)

# Get edge adjacency matrix (edges sharing a node are adjacent)
n_edges <- nrow(tree$edge)

# Build edge adjacency: edges share a node if one's child is another's parent
edge_adj <- matrix(0L, n_edges, n_edges)
for (i in 1:n_edges) {
  child_i <- tree$edge[i, 2]
  # Find edges where this child is the parent
  neighbors <- which(tree$edge[, 1] == child_i)
  edge_adj[i, neighbors] <- 1L
  edge_adj[neighbors, i] <- 1L
}

cat("Edge adjacency built:", sum(edge_adj)/2, "adjacent pairs\n")

# Summarize subs per edge by variant_sig_class and direction
edge_subs <- test[, .(
  n_sig_both = sum(variant_sig_class == "sig_both"),
  n_sig_both_warm = sum(variant_sig_class == "sig_both" & pheno_sign == 1),
  n_sig_both_cold = sum(variant_sig_class == "sig_both" & pheno_sign == -1),
  n_sig_both_concordant = sum(variant_sig_class == "sig_both" & concordant == TRUE, na.rm = TRUE),
  n_not_sig = sum(variant_sig_class == "not_sig"),
  n_not_sig_concordant = sum(variant_sig_class == "not_sig" & concordant == TRUE, na.rm = TRUE)
), by = edge_idx]

# Fill missing edges with zeros
all_edges <- data.table(edge_idx = 1:n_edges)
edge_subs <- merge(all_edges, edge_subs, by = "edge_idx", all.x = TRUE)
for (col in names(edge_subs)[-1]) {
  edge_subs[is.na(get(col)), (col) := 0L]
}

cat("Edges with sig_both subs:", sum(edge_subs$n_sig_both > 0), "\n")

# Moran's I function
morans_i <- function(x, W) {
  n <- length(x)
  x_bar <- mean(x)
  x_dev <- x - x_bar
  
  num <- sum(W * outer(x_dev, x_dev))
  denom <- sum(x_dev^2)
  W_sum <- sum(W)
  
  I <- (n / W_sum) * (num / denom)
  
  # Expected value under null
  E_I <- -1 / (n - 1)
  
  # Permutation test
  n_perm <- 99
  I_perm <- replicate(n_perm, {
    x_shuf <- sample(x)
    x_dev_shuf <- x_shuf - mean(x_shuf)
    num_shuf <- sum(W * outer(x_dev_shuf, x_dev_shuf))
    (n / W_sum) * (num_shuf / sum(x_dev_shuf^2))
  })
  
  p_val <- (sum(I_perm >= I) + 1) / (n_perm + 1)
  
  list(I = I, E_I = E_I, p = p_val)
}

# Run Moran's I for different metrics
cat("\n=== Moran's I: spatial autocorrelation of subs on tree ===\n")

cat("\nn_sig_both per edge:\n")
res <- morans_i(edge_subs$n_sig_both, edge_adj)
cat(sprintf("I=%.4f, E[I]=%.4f, p=%.4f\n", res$I, res$E_I, res$p))

cat("\nn_sig_both_warm per edge:\n")
res2 <- morans_i(edge_subs$n_sig_both_warm, edge_adj)
cat(sprintf("I=%.4f, E[I]=%.4f, p=%.4f\n", res2$I, res2$E_I, res2$p))

cat("\nn_sig_both_cold per edge:\n")
res3 <- morans_i(edge_subs$n_sig_both_cold, edge_adj)
cat(sprintf("I=%.4f, E[I]=%.4f, p=%.4f\n", res3$I, res3$E_I, res3$p))

cat("\nn_sig_both_concordant per edge:\n")
res4 <- morans_i(edge_subs$n_sig_both_concordant, edge_adj)
cat(sprintf("I=%.4f, E[I]=%.4f, p=%.4f\n", res4$I, res4$E_I, res4$p))
