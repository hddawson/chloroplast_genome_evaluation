#!/usr/bin/env Rscript
# Tree-based analysis of GWAS-selected sites
# Depth, Convergence, and Correlation

library(ape)
library(data.table)
library(Biostrings)
library(phangorn)  # for ancestral state reconstruction
library(arrow)
# PART 0: TREE QUALITY CHECKS
# Run these checks before proceeding with analysis

cat("\n=== TREE QUALITY CHECKS ===\n")

# --- Load tree ---
# Terminal: ensure tree file exists
# ls -la raxml_input/aa_tree_fiveParsimony.raxml.bestTree

data <- read_parquet("data/processed_data.parquet")

tree_file <- "raxml_input/aa_tree_fiveParsimony.raxml.bestTree"
stopifnot(file.exists(tree_file))

tree <- read.tree(tree_file)
stopifnot(inherits(tree, "phylo"))

# --- Basic stats ---
cat("\nBasic tree statistics:\n")
cat("  N tips:", Ntip(tree), "\n")
cat("  N internal nodes:", Nnode(tree), "\n")
cat("  Rooted:", is.rooted(tree), "\n")
cat("  Binary:", is.binary(tree), "\n")
cat("  Ultrametric:", is.ultrametric(tree), "\n")

# --- Branch length diagnostics ---
bl <- tree$edge.length
hist(bl)
cat("\nBranch length distribution:\n")
cat("  Min:", min(bl), "\n")
cat("  Max:", max(bl), "\n")
cat("  Median:", median(bl), "\n")
cat("  Mean:", mean(bl), "\n")
cat("  N near-zero (<1e-6):", sum(bl < 1e-6), "\n")
cat("  N very long (>1):", sum(bl > 1), "\n")

# RAxML reported 122 near-zero branches - verify
n_nearzero <- sum(bl < 1e-6)
#stopifnot(n_nearzero > 0)  # expect some given log warning

# --- Check for polytomies ---
if (!is.binary(tree)) {
  cat("\nWARNING: Tree has polytomies. Consider:\n")
  cat("  - Using multi2di() to resolve randomly\n")
  cat("  - This may affect ancestral state reconstruction\n")
  # Resolve for downstream analysis
  tree <- multi2di(tree, random = FALSE)
  cat("  - Resolved polytomies with zero-length branches\n")
}

# --- Root  ---
# For ASR, tree should ideally be rooted
# If unrooted, you need an outgroup or midpoint rooting
if (!is.rooted(tree)) {
  cat("\nWARNING: Tree is unrooted.\n")
  cat("Options for rooting:\n")
  cat("  1. Midpoint rooting: tree <- midpoint(tree)  # requires phangorn\n")
  cat("  2. Outgroup rooting: tree <- root(tree, outgroup='TAXON_ID')\n")
  cat("  3. Proceed unrooted (some analyses may fail)\n")
  
  # Midpoint root by default
  pinales <- intersect(tree$tip.label, data$ID[data$Order == "Pinales"])
  is.monophyletic(tree, pinales)
  
  tree <- root(tree, outgroup = pinales, resolve.root = TRUE)
  
  cat("Applied outgroup rooting w pinales\n")
}

is.rooted(tree)
# --- Check tip labels match data ---
# Load GWAS data to verify overlap
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

# Get tip labels from one GWAS result
sample_gwas <- readRDS(model_files[1])
gwas_ids <- sample_gwas[[1]]$IDs
tree_tips <- tree$tip.label

overlap <- length(intersect(gwas_ids, tree_tips))
cat("\nTip-GWAS overlap:\n")
cat("  GWAS samples:", length(gwas_ids), "\n")
cat("  Tree tips:", length(tree_tips), "\n")
cat("  Overlap:", overlap, "\n")

stopifnot(overlap > 0)  # must have some overlap

# --- Visual QC (save to file) ---
# Terminal: view histogram of branch lengths
par(mfrow = c(2, 2))

hist(log10(bl + 1e-10), breaks = 50, main = "Branch lengths (log10)",
     xlab = "log10(branch length)", col = "steelblue")
abline(v = log10(1e-6), col = "red", lty = 2)

hist(bl[bl > 1e-6], breaks = 50, main = "Non-zero branch lengths",
     xlab = "Branch length", col = "steelblue")

# Tip-to-root distances (depth)
tip_depths <- node.depth.edgelength(tree)[1:Ntip(tree)]
hist(tip_depths, breaks = 50, main = "Tip-to-root distances",
     xlab = "Distance", col = "steelblue")

# Node degree distribution
node_degrees <- tabulate(tree$edge[, 1])
barplot(table(node_degrees), main = "Node degree distribution",
        xlab = "N children", col = "steelblue")

par(mfrow=c(1,1))

# Save cleaned tree
saveRDS(tree, "results/tree_cleaned.rds")
cat("Saved cleaned tree to: results/tree_cleaned.rds\n")

# PART 1: LOAD AND PREPARE DATA

cat("\n=== LOADING DATA ===\n")

# --- Load cleaned tree ---
tree <- readRDS("results/tree_cleaned.rds")
stopifnot(inherits(tree, "phylo"))

# --- Load GWAS results ---
gwas_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    dt <- data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs,
      P_aa_only = m$P_aa_only,
      N = m$N,
      R2_full = m$R2_full
    )
    # Add effect sizes
    if (!is.null(m$effects)) {
      dt[, effects := list(list(m$effects))]
      dt[, residue_counts := list(list(m$residue_counts))]
    }
    dt
  }))
}))

cat("Loaded", nrow(gwas_results), "GWAS positions\n")

# --- Define foreground sites ---
# Using your existing threshold approach
gwas_results[, sig_class := {
  thresh_control <- quantile(P_aa_with_pcs, 0.05, na.rm = TRUE)
  thresh_nocontrol <- quantile(P_aa_only, 0.20, na.rm = TRUE)
  
  sig_ctrl <- P_aa_with_pcs < thresh_control
  sig_noctrl <- P_aa_only < thresh_nocontrol
  
  fcase(
    sig_ctrl & sig_noctrl, "sig_both",
    sig_ctrl & !sig_noctrl, "sig_control",
    !sig_ctrl & sig_noctrl, "sig_nocontrol",
    default = "not_sig"
  )
}]

thresh_control <- quantile(gwas_results$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(gwas_results$P_aa_only, 0.20, na.rm = TRUE)

# Colors for each class
cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

# Base plot
plot(-log10(gwas_results$P_aa_only), 
     -log10(gwas_results$P_aa_with_pcs),
     col = cols[gwas_results$sig_class],
     pch = 16, cex = 0.6,
     xlab = "-log10(P_aa_only)", 
     ylab = "-log10(P_aa_with_pcs)",
     main="Site classification")

# Threshold lines
abline(h = -log10(thresh_control), lty = 2, col = "darkorange")
abline(v = -log10(thresh_nocontrol), lty = 2, col = "steelblue")

# Counts annotation
counts <- table(gwas_results$sig_class)
legend_text <- paste0(names(counts), ": ", counts)
legend("topright", legend = legend_text, 
       fill = cols[names(counts)], bty = "n", cex = 0.8)


# Check
stopifnot(all(gwas_results$sig_class %in% c("not_sig", "sig_nocontrol", "sig_control", "sig_both")))
table(gwas_results$sig_class)
cat("Foreground sites (top 5%):", sum(gwas_results$is_foreground), "\n")
cat("Background sites:", sum(!gwas_results$is_foreground), "\n")

# PART 2: IDENTIFY VARIANTS OF INTEREST AT EACH SITE
# Not all variants at foreground sites matter - use effect sizes and p-values

cat("\n=== IDENTIFYING VARIANTS OF INTEREST ===\n")

# Function to identify significant variants at a site
identify_significant_variants <- function(effects_dt, p_thresh = 0.05, 
                                          effect_thresh = NULL) {
  if (is.null(effects_dt) || nrow(effects_dt) == 0) return(character(0))
  
  # Filter by p-value
  sig <- effects_dt[P_value < p_thresh]
  
  # Optionally filter by effect size
  if (!is.null(effect_thresh)) {
    sig <- sig[abs(Effect) > effect_thresh]
  }
  
  # Return residue names (clean up "X_aares_factor" prefix)
  residues <- gsub("^X_aa", "", sig$Residue)
  return(residues)
}

# Apply to foreground sites
foreground <- gwas_results[is_foreground == TRUE]

foreground[, significant_variants := lapply(effects, function(e) {
  identify_significant_variants(e, p_thresh = 0.05)
})]

foreground[, n_sig_variants := sapply(significant_variants, length)]

cat("Sites with significant variants:\n")
print(table(foreground$n_sig_variants))

# ---- raxml asr ---- 

# Read as tab-delimited, not fasta
anc_raw <- fread("raxml_input/aa_tree_ASR.raxml.ancestralStates", 
                 header = FALSE, sep = "\t", col.names = c("node", "seq"))

#load the partition map for the supermatrix
partitionMap <- readRDS("raxml_input/partitionMap.rds")
head(partitionMap)
geneBounds <- partitionMap[, .(start = min(GlobalPos), end = max(GlobalPos)), by = Gene]
stopifnot(nrow(geneBounds) > 0)


partLines <- geneBounds[, paste0("LG+G4, ", Gene, " = ", start, "-", end)]
writeLines(partLines, "raxml_input/partitions.txt")

# Extract positions 420 and 439 from each sequence
anc_raw[, site_420 := substr(seq, 7420, 7420)]
anc_raw[, site_439 := substr(seq, 7426, 7426)]

# Load the node-labeled tree
anc_tree <- read.tree("raxml_input/aa_tree_ASR.raxml.ancestralTree")

# Scaled tree-based analysis using RAxML ASR
# Requires: pre-loaded anc_raw, anc_tree, partitionMap, gwas_results
#!/usr/bin/env Rscript

# Scaled tree-based analysis using RAxML ASR
# Depth, Convergence, Correlation for all foreground sites


library(data.table)
library(ape)


# LOAD DATA


cat("\n=== LOADING DATA ===\n")

anc_raw <- fread("raxml_input/aa_tree_ASR.raxml.ancestralStates", 
                 header = FALSE, sep = "\t", col.names = c("node", "seq"))
anc_tree <- read.tree("raxml_input/aa_tree_ASR.raxml.ancestralTree")
partitionMap <- readRDS("raxml_input/partitionMap.rds")

stopifnot(nrow(anc_raw) > 0)
stopifnot(inherits(anc_tree, "phylo"))
stopifnot(nrow(partitionMap) > 0)

cat("Loaded:", nrow(anc_raw), "ancestral states\n")
cat("Tree tips:", Ntip(anc_tree), ", internal nodes:", Nnode(anc_tree), "\n")
cat("Partition map rows:", nrow(partitionMap), "\n")

# Node labels for tree
node_labels <- c(anc_tree$tip.label, anc_tree$node.label)
n_tips <- Ntip(anc_tree)


# IDENTIFY FOREGROUND SITES


cat("\n=== FOREGROUND SITES ===\n")

foreground <- partitionMap[gwas_hit == TRUE]
background <- partitionMap[gwas_hit == FALSE]

cat("Foreground sites:", nrow(foreground), "\n")
cat("Background sites:", nrow(background), "\n")

table(foreground$Gene)


# DEPTH ANALYSIS


cat("\n=== DEPTH ANALYSIS ===\n")

# Precompute tree depths once
topo_depths <- node.depth(anc_tree)
hist(topo_depths)
bl_depths <- node.depth.edgelength(anc_tree)
hist(bl_depths)

# Check what node names look like
head(anc_raw$node, 20)
tail(anc_raw$node, 20)

# Are tip labels in there?
sum(anc_tree$tip.label %in% anc_raw$node)

# What patterns exist?
head(grep("^Node", anc_raw$node, value = TRUE))
head(grep("^[^N]", anc_raw$node, value = TRUE))


depth_results <- list()

# Check if supermatrix exists
file.exists("raxml_input/superaa_collapsed.fasta")  # or whatever your AA supermatrix is called

# If so, load it
library(Biostrings)
supermat <- readAAStringSet("raxml_input/superaa_collapsed.fasta")  # adjust filename
names(supermat) <- sub("\\|.*", "", names(supermat))

# Check overlap with tree tips
sum(names(supermat) %in% anc_tree$tip.label)

# Convert supermatrix to character matrix for fast column access
supermat_mat <- as.matrix(supermat)
dim(supermat_mat)  # rows = tips, cols = positions

# Verify alignment length matches partition map
max(partitionMap$GlobalPos)
ncol(supermat_mat)
stopifnot(max(partitionMap$GlobalPos) <= ncol(supermat_mat))

# Test with a single position
gpos <- foreground$GlobalPos[1]

# Tip states from supermatrix
tip_states <- supermat_mat[, gpos]
names(tip_states) <- rownames(supermat_mat)

# Internal node states from anc_raw
node_states <- substr(anc_raw$seq, gpos, gpos)
names(node_states) <- anc_raw$node

# Combine
all_states <- c(tip_states, node_states)

# Check
length(all_states)
sum(names(all_states) %in% anc_tree$tip.label)  # should be 6094
sum(names(all_states) %in% anc_tree$node.label)  # should be ~6092

table(tip_states)
table(node_states)


depth_results <- list()

for (i in seq_len(nrow(foreground))) {
  site <- foreground[i]
  gpos <- site$GlobalPos
  
  # Tip states from supermatrix
  tip_states <- supermat_mat[, gpos]
  names(tip_states) <- rownames(supermat_mat)
  
  # Get variants (exclude gaps)
  variants <- unique(tip_states[tip_states != "-" & !is.na(tip_states)])
  if (length(variants) == 0) next
  
  for (v in variants) {
    tips_with <- names(tip_states)[tip_states == v & !is.na(tip_states)]
    tips_with <- intersect(tips_with, anc_tree$tip.label)
    
    if (length(tips_with) < 2) next
    
    mrca_node <- getMRCA(anc_tree, tips_with)
    if (is.null(mrca_node)) next
    
    # Crown depth
    subtree <- tryCatch(extract.clade(anc_tree, mrca_node), error = function(e) NULL)
    crown <- if (!is.null(subtree)) max(node.depth.edgelength(subtree)) else NA_real_
    
    depth_results[[length(depth_results) + 1]] <- data.table(
      Gene = site$Gene,
      Position = site$GenePos,
      GlobalPos = gpos,
      variant = v,
      n_tips = length(tips_with),
      mrca_node = mrca_node,
      topo_depth = topo_depths[mrca_node],
      bl_depth = bl_depths[mrca_node],
      crown_depth = crown
    )
  }
  
  if (i %% 50 == 0) cat("  Depth:", i, "/", nrow(foreground), "\n")
}

depth_dt <- rbindlist(depth_results)
cat("Depth results:", nrow(depth_dt), "rows\n")

# Quick look
summary(depth_dt$bl_depth)
hist(depth_dt$bl_depth, breaks = 50, main = "Branch length depth of variant MRCAs")
boxplot(topo_depth ~ Gene, depth_dt, las = 2, main = "Depth by gene")
hist(depth_dt$topo_depth)
hist(depth_dt$crown_depth)


#now w sampling background 

depth_results <- list()

bg_sample <- background[sample(.N, nrow(foreground))]  # match foreground size

all_sites <- rbind(
  foreground[, .(Gene, GenePos, GlobalPos, class = "foreground")],
  bg_sample[, .(Gene, GenePos, GlobalPos, class = "background")]
)


for (i in seq_len(nrow(all_sites))) {
  site <- all_sites[i]
  gpos <- site$GlobalPos
  
  tip_states <- supermat_mat[, gpos]
  names(tip_states) <- rownames(supermat_mat)
  
  variants <- unique(tip_states[tip_states != "-" & !is.na(tip_states)])
  if (length(variants) == 0) next
  
  for (v in variants) {
    tips_with <- names(tip_states)[tip_states == v & !is.na(tip_states)]
    tips_with <- intersect(tips_with, anc_tree$tip.label)
    
    if (length(tips_with) < 2) next
    
    mrca_node <- getMRCA(anc_tree, tips_with)
    if (is.null(mrca_node)) next
    
    subtree <- tryCatch(extract.clade(anc_tree, mrca_node), error = function(e) NULL)
    crown <- if (!is.null(subtree)) max(node.depth.edgelength(subtree)) else NA_real_
    
    depth_results[[length(depth_results) + 1]] <- data.table(
      Gene = site$Gene,
      Position = site$GenePos,
      GlobalPos = gpos,
      class = site$class,
      variant = v,
      n_tips = length(tips_with),
      mrca_node = mrca_node,
      topo_depth = topo_depths[mrca_node],
      bl_depth = bl_depths[mrca_node],
      crown_depth = crown
    )
  }
  
  if (i %% 500 == 0) cat("  Depth:", i, "/", nrow(all_sites), "\n")
}

depth_dt <- rbindlist(depth_results)
cat("Depth results:", nrow(depth_dt), "rows\n")
table(depth_dt$class)

boxplot(bl_depth ~ class, depth_dt, main = "BL depth: foreground vs background")
wilcox.test(bl_depth ~ class, depth_dt)

# ---- CONVERGENCE ANALYSIS ----

# Build edge table once
edges <- anc_tree$edge

conv_results <- list()

# Combine foreground and background
all_sites <- rbind(
  foreground[, .(Gene, GenePos, GlobalPos, class = "foreground")],
  background[, .(Gene, GenePos, GlobalPos, class = "background")]
)

# Build node_labels lookup once
node_labels <- c(anc_tree$tip.label, anc_tree$node.label)
edges <- anc_tree$edge

for (i in seq_len(nrow(all_sites))) {
  site <- all_sites[i]
  gpos <- site$GlobalPos
  
  # Tip states from supermatrix
  tip_states <- supermat_mat[, gpos]
  names(tip_states) <- rownames(supermat_mat)
  
  # Internal node states from anc_raw
  node_states <- substr(anc_raw$seq, gpos, gpos)
  names(node_states) <- anc_raw$node
  
  # Combine
  all_states <- c(tip_states, node_states)
  
  # Get parent/child states for each edge
  parent_states <- all_states[node_labels[edges[, 1]]]
  child_states <- all_states[node_labels[edges[, 2]]]
  
  # Variants from tips only (exclude gaps)
  variants <- unique(tip_states[tip_states != "-" & !is.na(tip_states)])
  if (length(variants) == 0) next
  
  for (v in variants) {
    n_tips <- sum(tip_states == v, na.rm = TRUE)
    # Count transitions TO this state
    n_origins <- sum(child_states == v & parent_states != v, na.rm = TRUE)
    
    conv_results[[length(conv_results) + 1]] <- data.table(
      Gene = site$Gene,
      Position = site$GenePos,
      GlobalPos = gpos,
      class = site$class,
      variant = v,
      n_tips = n_tips,
      n_origins = n_origins,
      convergence_ratio = if (n_tips > 0) n_origins / n_tips else NA_real_
    )
  }
  
  if (i %% 500 == 0) cat("  Convergence:", i, "/", nrow(all_sites), "\n")
}

conv_dt <- rbindlist(conv_results)
cat("Convergence results:", nrow(conv_dt), "rows\n")
table(conv_dt$class)

# Quick look
summary(conv_dt$n_origins)
hist(conv_dt$n_origins)
hist(conv_dt$n_origins, breaks = 50, main = "Number of independent origins per variant")

par(mfrow=c(1,2))
# Compare foreground vs background
boxplot(n_origins ~ class, conv_dt,
        main = "Independent origins: foreground vs background",
        ylim=c(0,50))
boxplot(n_origins ~ class, conv_dt,
        main = "Independent origins: foreground vs background")
wilcox.test(n_origins ~ class, conv_dt)
par(mfrow=c(1,1))

boxplot(convergence_ratio ~ class, conv_dt, main = "Convergence ratio: foreground vs background")
wilcox.test(convergence_ratio ~ class, conv_dt)

# By gene (foreground only)
boxplot(n_origins ~ Gene, conv_dt[class == "foreground"], las = 2, main = "Origins by gene (foreground)")

# High convergence variants
conv_dt[n_origins > 10][order(-n_origins)]


site_conv <- conv_dt[, .(
  total_origins = sum(n_origins),
  total_tips = sum(n_tips),
  n_variants = .N,
  max_origins = max(n_origins),
  weighted_conv_ratio = sum(n_origins) / sum(n_tips)
), by = .(Gene, Position, GlobalPos, class)]

# Control for number of variants and total tips
site_conv[, origins_per_variant := total_origins / n_variants]

# Compare
boxplot(origins_per_variant ~ class, site_conv, 
        main = "Origins per variant (controls for polymorphism level)")
wilcox.test(origins_per_variant ~ class, site_conv)

boxplot(weighted_conv_ratio ~ class, site_conv,
        main = "Weighted convergence ratio")
wilcox.test(weighted_conv_ratio ~ class, site_conv)

# Logistic regression: does convergence predict foreground status?
site_conv[, is_foreground := as.integer(class == "foreground")]
glm_conv <- glm(is_foreground ~ origins_per_variant + n_variants + total_tips, 
                data = site_conv, family = binomial)
summary(glm_conv)

# ---- pPAIRWISE SITE CORRELATION (HAPLOTYPE STRUCTURE) ----
# PAIRWISE SITE CORRELATION (HAPLOTYPE STRUCTURE)

cat("\n=== PAIRWISE SITE CORRELATION ===\n")

# Cramér's V for categorical correlation
cramers_v <- function(x, y) {
  valid <- !is.na(x) & !is.na(y) & x != "-" & y != "-"
  if (sum(valid) < 100) return(NA_real_)
  
  tab <- table(x[valid], y[valid])
  chi2 <- tryCatch(chisq.test(tab)$statistic, error = function(e) NA)
  if (is.na(chi2)) return(NA_real_)
  
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  if (k <= 1) return(NA_real_)
  
  sqrt(chi2 / (n * (k - 1)))
}

# --- WITHIN FOREGROUND: find haplotype blocks ---
cat("\n--- Foreground site correlations ---\n")

fg_positions <- foreground$GlobalPos
n_fg <- length(fg_positions)
cat("Foreground sites:", n_fg, "\n")

# Sample pairs if too many
set.seed(123)
n_pairs <- min(50000, choose(n_fg, 2))
pair_idx <- t(combn(n_fg, 2))
if (nrow(pair_idx) > n_pairs) {
  pair_idx <- pair_idx[sample(nrow(pair_idx), n_pairs), ]
}

cat("Computing", nrow(pair_idx), "foreground pairs...\n")

fg_cor_results <- list()
for (p in seq_len(nrow(pair_idx))) {
  i <- pair_idx[p, 1]
  j <- pair_idx[p, 2]
  
  x <- supermat_mat[, fg_positions[i]]
  y <- supermat_mat[, fg_positions[j]]
  
  v <- cramers_v(x, y)
  
  if (!is.na(v) && v >= 0.3) {
    fg_cor_results[[length(fg_cor_results) + 1]] <- data.table(
      pos1 = fg_positions[i],
      pos2 = fg_positions[j],
      cramers_v = v,
      class = "foreground"
    )
  }
  
  if (p %% 10000 == 0) cat("  Pairs:", p, "/", nrow(pair_idx), "\n")
}

fg_cor_dt <- rbindlist(fg_cor_results)
cat("Foreground correlated pairs (V >= 0.3):", nrow(fg_cor_dt), "\n")

# --- WITHIN BACKGROUND (sample to match) ---
cat("\n--- Background site correlations ---\n")

bg_positions <- background$GlobalPos
n_bg <- length(bg_positions)

# Sample same number of pairs as foreground
set.seed(456)
bg_sample_idx <- sample(n_bg, min(n_fg, n_bg))
bg_positions_sample <- bg_positions[bg_sample_idx]

pair_idx_bg <- t(combn(length(bg_positions_sample), 2))
if (nrow(pair_idx_bg) > n_pairs) {
  pair_idx_bg <- pair_idx_bg[sample(nrow(pair_idx_bg), n_pairs), ]
}

cat("Computing", nrow(pair_idx_bg), "background pairs...\n")

bg_cor_results <- list()
for (p in seq_len(nrow(pair_idx_bg))) {
  i <- pair_idx_bg[p, 1]
  j <- pair_idx_bg[p, 2]
  
  x <- supermat_mat[, bg_positions_sample[i]]
  y <- supermat_mat[, bg_positions_sample[j]]
  
  v <- cramers_v(x, y)
  
  if (!is.na(v) && v >= 0.3) {
    bg_cor_results[[length(bg_cor_results) + 1]] <- data.table(
      pos1 = bg_positions_sample[i],
      pos2 = bg_positions_sample[j],
      cramers_v = v,
      class = "background"
    )
  }
  
  if (p %% 10000 == 0) cat("  Pairs:", p, "/", nrow(pair_idx_bg), "\n")
}

bg_cor_dt <- rbindlist(bg_cor_results)
cat("Background correlated pairs (V >= 0.3):", nrow(bg_cor_dt), "\n")

# --- COMPARE ---
cor_dt <- rbind(fg_cor_dt, bg_cor_dt)

cat("\n--- Comparison ---\n")
cat("Foreground pairs V >= 0.3:", nrow(fg_cor_dt), "\n")
cat("Background pairs V >= 0.3:", nrow(bg_cor_dt), "\n")

# Rate of strong correlations
fg_rate <- nrow(fg_cor_dt) / nrow(pair_idx)
bg_rate <- nrow(bg_cor_dt) / nrow(pair_idx_bg)
cat("Foreground correlation rate:", round(fg_rate, 4), "\n")
cat("Background correlation rate:", round(bg_rate, 4), "\n")

# Distribution of V values
if (nrow(cor_dt) > 0) {
  boxplot(cramers_v ~ class, cor_dt, 
          main = "Cramér's V: foreground vs background (V >= 0.3 only)")
  
  hist(fg_cor_dt$cramers_v, breaks = 30, col = rgb(1, 0, 0, 0.5),
       main = "Distribution of Cramér's V", xlab = "Cramér's V")
  hist(bg_cor_dt$cramers_v, breaks = 30, col = rgb(0, 0, 1, 0.5), add = TRUE)
  legend("topright", c("Foreground", "Background"), 
         fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))
}

# Top correlated foreground pairs
if (nrow(fg_cor_dt) > 0) {
  setorder(fg_cor_dt, -cramers_v)
  cat("\nTop 20 foreground correlated pairs:\n")
  print(head(fg_cor_dt, 20))
}




# ---- phytools asr ---- 

# PART 3: ANCESTRAL STATE RECONSTRUCTION (ASR)#
# This is computationally intensive for large trees

cat("\n=== ANCESTRAL STATE RECONSTRUCTION ===\n")

# --- Load alignment data ---
# You'll need the actual alignment to get character states

perform_asr <- function(tree, gene, position, aln_file, method = "parsimony") {
  # Read alignment
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*", "", names(aln))
  
  # Get character states at this position
  aln_mat <- as.matrix(aln)
  
  stopifnot(position <= ncol(aln_mat))
  
  states <- aln_mat[, position]
  names(states) <- names(aln)
  
  # Prune tree to match alignment
  common_tips <- intersect(tree$tip.label, names(states))
  if (length(common_tips) < 100) {
    warning("Too few overlapping tips for ASR")
    return(NULL)
  }
  
  pruned_tree <- keep.tip(tree, common_tips)
  states <- states[common_tips]
  
  # Convert to phyDat format
  states_factor <- factor(states)
  phy_dat <- phyDat(as.matrix(states_factor), type = "USER", 
                    levels = levels(states_factor))
  
  # Perform ASR
  if (method == "ML") {
    # Maximum parsimony (faster) or ML
    anc <- ancestral.pars(pruned_tree, phy_dat, type = "MPR")
  } else if (method == "parsimony") {
    anc <- ancestral.pars(pruned_tree, phy_dat, type = "ACCTRAN")
  }
  
  return(list(
    tree = pruned_tree,
    tip_states = states,
    ancestral = anc,
    levels = levels(states_factor)
  ))
}

# Example usage (do not run on all sites - too slow)
asr_result <- perform_asr(tree, "rbcL", 439, 
                          "data/tmp/alignedGenes/rbcL_AA_aligned.fasta")

# Extract states directly from the phyDat object
anc <- asr_result$ancestral
pruned_tree <- asr_result$tree

# Get state labels
state_labels <- attr(anc, "levels")

# Extract as list, then convert
anc_list <- lapply(anc, function(x) {
  # Each element is a matrix of probabilities (or 0/1 for parsimony)
  # Rows = sites (just 1), Cols = states
  if (is.matrix(x)) {
    state_labels[which.max(x[1, ])]
  } else {
    state_labels[which.max(x)]
  }
})
anc_states <- unlist(anc_list)
names(anc_states) <- names(anc)

# Split into tips and internal nodes
n_tips <- Ntip(pruned_tree)
tip_states <- anc_states[1:n_tips]
node_states <- anc_states[(n_tips + 1):length(anc_states)]

# Quick summary
cat("Tip states:\n")
print(table(tip_states))

cat("\nInternal node states:\n")
print(table(node_states))

# Count transitions TO each state (approximate origins)
edges <- pruned_tree$edge
parent_states <- anc_states[edges[, 1]]
child_states <- anc_states[edges[, 2]]

transitions <- data.table(
  parent = parent_states,
  child = child_states
)
transitions[, is_change := parent != child]

cat("\nState transitions:\n")
print(transitions[is_change == TRUE, .N, by = .(from = parent, to = child)])

# Example usage (do not run on all sites - too slow)
asr_result <- perform_asr(tree, "rbcL", 420, 
                          "data/tmp/alignedGenes/rbcL_AA_aligned.fasta")

# Extract states directly from the phyDat object
anc <- asr_result$ancestral
pruned_tree <- asr_result$tree

# Get state labels
state_labels <- attr(anc, "levels")

# Extract as list, then convert
anc_list <- lapply(anc, function(x) {
  # Each element is a matrix of probabilities (or 0/1 for parsimony)
  # Rows = sites (just 1), Cols = states
  if (is.matrix(x)) {
    state_labels[which.max(x[1, ])]
  } else {
    state_labels[which.max(x)]
  }
})
anc_states <- unlist(anc_list)
names(anc_states) <- names(anc)

# Split into tips and internal nodes
n_tips <- Ntip(pruned_tree)
tip_states <- anc_states[1:n_tips]
node_states <- anc_states[(n_tips + 1):length(anc_states)]

# Quick summary
cat("Tip states:\n")
print(table(tip_states))

cat("\nInternal node states:\n")
print(table(node_states))

# Count transitions TO each state (approximate origins)
edges <- pruned_tree$edge
parent_states <- anc_states[edges[, 1]]
child_states <- anc_states[edges[, 2]]

transitions <- data.table(
  parent = parent_states,
  child = child_states
)
transitions[, is_change := parent != child]

cat("\nState transitions:\n")
print(transitions[is_change == TRUE, .N, by = .(from = parent, to = child)])

# Get tip states for both sites
asr_439 <- perform_asr(tree, "rbcL", 439, "data/tmp/alignedGenes/rbcL_AA_aligned.fasta")
asr_420 <- perform_asr(tree, "rbcL", 420, "data/tmp/alignedGenes/rbcL_AA_aligned.fasta")

# Extract tip states using the workaround
extract_tip_states <- function(asr_result) {
  anc <- asr_result$ancestral
  state_labels <- attr(anc, "levels")
  anc_list <- lapply(anc, function(x) {
    if (is.matrix(x)) state_labels[which.max(x[1, ])] else state_labels[which.max(x)]
  })
  states <- unlist(anc_list)
  names(states) <- names(anc)
  states[1:Ntip(asr_result$tree)]
}

tips_439 <- extract_tip_states(asr_439)
tips_420 <- extract_tip_states(asr_420)

# Align to common tips
common <- intersect(names(tips_439), names(tips_420))
stopifnot(length(common) > 1000)

# Contingency table
tab <- table(site_439 = tips_439[common], site_420 = tips_420[common])
print(tab)

# Chi-squared test
chisq.test(tab)

# Cramér's V for effect size
n <- sum(tab)
chi2 <- chisq.test(tab)$statistic
k <- min(nrow(tab), ncol(tab))
cramers_v <- sqrt(chi2 / (n * (k - 1)))
cat("\nCramér's V:", round(cramers_v, 3), "\n")

# Get edge-level transitions for both sites
get_edge_transitions <- function(asr_result) {
  anc <- asr_result$ancestral
  pruned_tree <- asr_result$tree
  state_labels <- attr(anc, "levels")
  
  anc_list <- lapply(anc, function(x) {
    if (is.matrix(x)) state_labels[which.max(x[1, ])] else state_labels[which.max(x)]
  })
  anc_states <- unlist(anc_list)
  names(anc_states) <- names(anc)
  
  edges <- pruned_tree$edge
  data.table(
    edge_id = seq_len(nrow(edges)),
    parent_node = edges[, 1],
    child_node = edges[, 2],
    parent_state = anc_states[edges[, 1]],
    child_state = anc_states[edges[, 2]]
  )
}

trans_439 <- get_edge_transitions(asr_439)
trans_420 <- get_edge_transitions(asr_420)

# Flag transitions
trans_439[, change_439 := paste0(parent_state, "->", child_state)]
trans_439[, has_change_439 := parent_state != child_state]

trans_420[, change_420 := paste0(parent_state, "->", child_state)]
trans_420[, has_change_420 := parent_state != child_state]

# Merge on edge_id
merged <- merge(trans_439[, .(edge_id, change_439, has_change_439)],
                trans_420[, .(edge_id, change_420, has_change_420)],
                by = "edge_id")

# Do transitions co-occur on same branches?
table(merged$has_change_439, merged$has_change_420, 
      dnn = c("change_328", "change_309"))

# Which specific transitions co-occur?
merged[has_change_439 & has_change_420, .N, by = .(change_439, change_420)]

# ---- variant analysis ---- 
# Build comprehensive variant-level results table
# Requires: anc_tree, supermat_mat, anc_raw, partitionMap, gwas_results already loaded

cat("\n=== BUILDING VARIANT-LEVEL TABLE ===\n")

# Precompute tree info
topo_depths <- node.depth(anc_tree)
bl_depths <- node.depth.edgelength(anc_tree)
root_depths <- max(bl_depths) - bl_depths  # depth from root
node_labels <- c(anc_tree$tip.label, anc_tree$node.label)
edges <- anc_tree$edge
n_tips <- Ntip(anc_tree)

summary(topo_depths)
summary(bl_depths)
hist(bl_depths)
summary(root_depths)
hist(root_depths)
summary(edges)
summary(n_tips)


# Merge partitionMap with gwas_results to get p-values and effects
stopifnot("GlobalPos" %in% names(partitionMap))
stopifnot("Position" %in% names(gwas_results))

# Create key for merging
gwas_results[, merge_key := paste(Gene, Position, sep = "_")]
partitionMap[, merge_key := paste(Gene, GenePos, sep = "_")]

site_info <- merge(
  partitionMap[, .(GlobalPos, Gene, GenePos, merge_key)],
  gwas_results[, .(merge_key, P_aa_only, P_aa_with_pcs, sig_class, effects, residue_counts)],
  by = "merge_key", all.x = TRUE
)

stopifnot(nrow(site_info) > 0)
cat("Sites with GWAS info:", sum(!is.na(site_info$P_aa_only)), "/", nrow(site_info), "\n")

# Main loop over all sites with GWAS data
variant_results <- list()

sites_to_process <- site_info[!is.na(P_aa_only)]
cat("Processing", nrow(sites_to_process), "sites...\n")

for (i in seq_len(nrow(sites_to_process))) {
  site <- sites_to_process[i]
  gpos <- site$GlobalPos
  
  # Tip states
  tip_states <- supermat_mat[, gpos]
  names(tip_states) <- rownames(supermat_mat)
  
  # Internal node states
  node_states <- substr(anc_raw$seq, gpos, gpos)
  names(node_states) <- anc_raw$node
  
  all_states <- c(tip_states, node_states)
  
  # Parent/child states for edges
  parent_states <- all_states[node_labels[edges[, 1]]]
  child_states <- all_states[node_labels[edges[, 2]]]
  
  # Get variants (exclude gaps)
  variants <- unique(tip_states[tip_states != "-" & !is.na(tip_states)])
  if (length(variants) == 0) next
  
  # Extract per-residue p-values and effects
  eff_dt <- site$effects[[1]]
  
  for (v in variants) {
    tips_with <- names(tip_states)[tip_states == v & !is.na(tip_states)]
    tips_with <- intersect(tips_with, anc_tree$tip.label)
    n_tip <- length(tips_with)
    
    if (n_tip < 1) next
    
    # Count independent origins
    n_origins <- sum(child_states == v & parent_states != v, na.rm = TRUE)
    
    # Identify parent alleles (what states transition TO this variant)
    parent_alleles <- unique(parent_states[child_states == v & parent_states != v])
    parent_alleles <- parent_alleles[!is.na(parent_alleles) & parent_alleles != "-"]
    parent_allele_str <- if (length(parent_alleles) > 0) paste(sort(parent_alleles), collapse = ",") else NA_character_
    
    # MRCA and depths (only if >= 2 tips)
    if (n_tip >= 2) {
      mrca_node <- getMRCA(anc_tree, tips_with)
      topo_d <- topo_depths[mrca_node]
      bl_d <- bl_depths[mrca_node]
      root_d <- root_depths[mrca_node]
    } else {
      mrca_node <- NA_integer_
      topo_d <- NA_real_
      bl_d <- NA_real_
      root_d <- NA_real_
    }
    
    # Per-residue stats from effects table
    p_residue <- NA_real_
    residue_effect <- NA_real_
    if (!is.null(eff_dt) && nrow(eff_dt) > 0) {
      # Match residue name (effects table has "X_aa<residue>" format)
      residue_row <- eff_dt[grepl(paste0(v, "$"), Residue)]
      if (nrow(residue_row) == 1) {
        p_residue <- residue_row$P_value
        residue_effect <- residue_row$Effect
      }
    }
    
    variant_results[[length(variant_results) + 1]] <- data.table(
      Gene = site$Gene,
      Position = site$GenePos,
      GlobalPos = gpos,
      variant = v,
      parent_alleles = parent_allele_str,
      n_tips = n_tip,
      n_origins = n_origins,
      convergence_ratio = n_origins / n_tip,
      mrca_node = mrca_node,
      topo_depth = topo_d,
      bl_depth = bl_d,
      root_depth = root_d,
      P_aa_only = site$P_aa_only,
      P_aa_with_pcs = site$P_aa_with_pcs,
      P_residue = p_residue,
      residue_effect = residue_effect,
      sig_class = site$sig_class
    )
  }
  
  if (i %% 500 == 0) cat("  Processed:", i, "/", nrow(sites_to_process), "\n")
}




variant_dt <- rbindlist(variant_results)
cat("\nFinal table:", nrow(variant_dt), "rows\n")

# Basic checks
stopifnot(nrow(variant_dt) > 0)
stopifnot(all(c("variant", "parent_alleles", "n_origins", "sig_class") %in% names(variant_dt)))

# Summary
cat("\nBy sig_class:\n")
print(variant_dt[, .N, by = sig_class])

cat("\nVariants with residue-level p-values:\n")
cat(sum(!is.na(variant_dt$P_residue)), "/", nrow(variant_dt), "\n")

clean_var_dt <- variant_dt[!is.na(variant_dt$P_residue),]

qqnorm(-log10(clean_var_dt$P_residue), main = "QQ plot of variant p-values")
qqline(-log10(clean_var_dt$P_residue))
abline(v=-log10(0.05))
clean_var_dt$variant_sig_class <- clean_var_dt$sig_class
table(clean_var_dt$variant_sig_class)
clean_var_dt$variant_sig_class[clean_var_dt$P_residue > 0.05] <- "not_sig"
table(clean_var_dt$variant_sig_class)

summary(clean_var_dt)
plot(clean_var_dt$n_origins, clean_var_dt$convergence_ratio)
boxplot(n_origins ~ variant_sig_class, clean_var_dt)
boxplot(convergence_ratio ~ sig_class, clean_var_dt)

hist(clean_var_dt$n_tips)
clean_var_dt <- clean_var_dt[n_tips >= 10]

boxplot(bl_depth ~ sig_class, clean_var_dt)

# ---- convergence by gene class ---- 

# Colors for each class
cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

# Base plot
plot(-log10(clean_var_dt$P_aa_only), 
     -log10(clean_var_dt$P_aa_with_pcs),
     col = cols[clean_var_dt$sig_class],
     pch = 16, cex = 0.6,
     xlab = "-log10(P_aa_only)", 
     ylab = "-log10(P_aa_with_pcs)",
     main="Variant classification")

# Threshold lines
abline(h = -log10(thresh_control), lty = 2, col = "darkorange")
abline(v = -log10(thresh_nocontrol), lty = 2, col = "steelblue")

# Counts annotation
counts <- table(clean_var_dt$sig_class)
legend_text <- paste0(names(counts), ": ", counts)
legend("topright", legend = legend_text, 
       fill = cols[names(counts)], bty = "n", cex = 0.8)



# Boxplot with significance annotations
cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

# Get levels in consistent order
class_order <- c("not_sig", "sig_nocontrol", "sig_control", "sig_both")
clean_var_dt[, sig_class := factor(sig_class, levels = class_order)]

# Boxplot
boxplot(n_origins ~ sig_class, clean_var_dt, 
        col = cols[class_order],
        main = "Variant origins by significance class",
        ylab = "Variant origins", xlab = "")

library(ggplot2)
library(ggsignif)

cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

class_order <- c("not_sig", "sig_nocontrol", "sig_control", "sig_both")
clean_var_dt[, variant_sig_class := factor(variant_sig_class, levels = class_order)]



all_comparisons <- combn(class_order, 2, simplify = FALSE)

ggplot(clean_var_dt, aes(x = variant_sig_class, y = n_origins, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = cols) +
  geom_signif(
    comparisons = all_comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.1
  ) +
  labs(x = "", y = "Number of origins", title = "Independent origins by significance class") +
  theme_bw() +
  theme(legend.position = "none")

# Option 2: violin + boxplot combo (shows distribution shape)
ggplot(clean_var_dt, aes(x = variant_sig_class, y = n_origins, fill = variant_sig_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  #coord_cartesian(ylim = c(0, quantile(clean_var_dt$n_origins, 0.95))) +
  scale_fill_manual(values = cols) +
  geom_signif(
    comparisons = all_comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.08,
    y_position = 300 * c(1.05, 1.15, 1.25, 1.35, 1.45, 1.55)
  ) +
  labs(x = "", y = "Number of origins ", 
       title = "Independent origins by significance class") +
  theme_bw() +
  theme(legend.position = "none")

# Option 2: violin + boxplot combo (shows distribution shape)
ggplot(clean_var_dt, aes(x = variant_sig_class, y = convergence_ratio, fill = variant_sig_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_signif(
    comparisons = all_comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.08,
    y_position = quantile(clean_var_dt$convergence_ratio, 0.95) * c(1.05, 1.15, 1.25, 1.35, 1.45, 1.55)
  ) +
  labs(x = "", y = "Convergence Ratio (95th percentile view)", 
       title = "Convergence Ratio by significance class") +
  theme_bw() +
  theme(legend.position = "none")

# Option 2: violin + boxplot combo (shows distribution shape)
ggplot(clean_var_dt, aes(x = variant_sig_class, y = convergence_ratio, fill = variant_sig_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(clean_var_dt$convergence_ratio, 0.95))) +
  scale_fill_manual(values = cols) +
  geom_signif(
    comparisons = all_comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.08,
    y_position = quantile(clean_var_dt$convergence_ratio, 0.95) * c(1.05, 1.15, 1.25, 1.35, 1.45, 1.55)
  ) +
  labs(x = "", y = "Convergence Ratio", 
       title = "Convergence Ratio by significance class") +
  theme_bw() +
  theme(legend.position = "none")

plot(clean_var_dt$bl_depth, clean_var_dt$n_origins)
plot(clean_var_dt$bl_depth, clean_var_dt$convergence_ratio)

# ---- depth vs convergence ---- 
# Plot 1: bl_depth vs n_origins
ggplot(clean_var_dt, aes(x = bl_depth, y = n_origins, color = variant_sig_class)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = cols) +
  facet_wrap(~variant_sig_class) +
  labs(x = "Branch length depth", y = "Number of origins",
       title = "Convergence vs depth by significance class") +
  theme_bw() +
  theme(legend.position = "none")

# Plot 2: bl_depth vs convergence_ratio
ggplot(clean_var_dt, aes(x = bl_depth, y = convergence_ratio, color = variant_sig_class)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = cols) +
  #scale_y_sqrt(breaks = c(0, 0.5, 1, 2, 5, 10)) +
  facet_wrap(~variant_sig_class) +
  labs(x = "Branch length depth", y = "Convergence ratio",
       title = "Convergence ratio vs depth by significance class") +
  theme_bw() +
  theme(legend.position = "none")

hist(clean_var_dt$residue_effect)
boxplot(residue_effect ~ variant_sig_class, clean_var_dt)

plot(clean_var_dt$residue_effect, clean_var_dt$n_origins)

# ---- variant pca ----
pca_data <- clean_var_dt[variant_sig_class != "not_sig", c("n_origins", "n_tips", "convergence_ratio", "bl_depth", "residue_effect")]

pca <- prcomp(pca_data, scale.=T, center=T)
pvar <- pca$sdev^2 / sum(pca$sdev^2)
pca
barplot(pvar)
plot(pca$x[,1],pca$x[,2])
plot(pca$x[,1],pca$x[,3])
pca_var_dt <- clean_var_dt[variant_sig_class != "not_sig",]
# Add PCA scores to data
pca_var_dt[, `:=`(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3])]

# Gene prefix (first 3 letters)
pca_var_dt[, gene_prefix := substr(Gene, 1, 3)]

# --- 1. By significance class ---
p_sig <- ggplot(pca_var_dt, aes(x = PC1, y = PC2, color = sig_class)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = cols) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA by significance class") +
  theme_bw()
p_sig
# --- 2. By amino acid ---
# Count per AA to find common ones
aa_counts <- pca_var_dt[, .N, by = variant][order(-N)]
top_aa <- aa_counts[N > 50, variant]  # adjust threshold as needed

p_aa <- ggplot(pca_var_dt[variant %in% top_aa], 
               aes(x = PC1, y = PC2, color = variant)) +
  geom_point(alpha = 0.5, size = 1) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA by amino acid (common only)") +
  theme_bw()

p_aa

ggplot(pca_var_dt[variant %in% top_aa], 
       aes(x = PC1, y = PC3, color = variant)) +
  geom_point(alpha = 0.5, size = 1) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA by amino acid (common only)") +
  theme_bw()

pca_var_dt$AminoAcid <- as.factor(pca_var_dt$variant)
boxplot(residue_effect ~ AminoAcid, pca_var_dt[variant_sig_class != "not_sig"],
        ylim = c(-10,10))
# AA summary stats
aa_summary <- pca_var_dt[, .(
  mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
  mean_effect = mean(residue_effect, na.rm = TRUE),
  n = .N
), by = variant][order(-abs(mean_PC1))]
print("Amino acids by mean PC1:")
print(head(aa_summary, 10))

# --- 3. By effect direction (hot/cold) ---
pca_var_dt[, effect_dir := fcase(
  residue_effect > 0, "hot",
  residue_effect < 0, "cold",
  default = "neutral"
)]

p_effect <- ggplot(pca_var_dt[effect_dir != "neutral"], 
                   aes(x = PC1, y = PC2, color = effect_dir)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c(cold = "blue", hot = "red")) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA by effect direction") +
  theme_bw()
p_effect

# --- 4. By gene prefix ---
gene_counts <- pca_var_dt[, .N, by = gene_prefix][order(-N)]
top_genes <- gene_counts[N > 100, gene_prefix]  # adjust threshold

p_gene <- ggplot(pca_var_dt[gene_prefix %in% top_genes], 
                 aes(x = PC1, y = PC2, color = gene_prefix)) +
  geom_point(alpha = 0.5, size = 1) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA by gene (common only)") +
  theme_bw()

p_gene
# Gene summary stats
gene_summary <- pca_var_dt[, .(
  mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
  mean_origins = mean(n_origins),
  mean_depth = mean(bl_depth, na.rm = TRUE),
  n = .N
), by = gene_prefix][order(-abs(mean_PC1))]
print("Genes by mean PC1:")
print(head(gene_summary, 10))

# --- 5. Combined facet view ---
p_facet <- ggplot(pca_var_dt, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sig_class), alpha = 0.3, size = 0.8) +
  scale_color_manual(values = cols) +
  facet_wrap(~gene_prefix) +
  labs(title = "PCA faceted by gene") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_sig)
print(p_aa)
print(p_effect)
print(p_gene)
print(p_facet)

# --- Statistical tests ---
cat("\n--- ANOVA: PC1 ~ sig_class ---\n")
print(summary(aov(PC1 ~ sig_class, data = pca_var_dt)))

cat("\n--- ANOVA: PC1 ~ gene_prefix ---\n")
print(summary(aov(PC1 ~ gene_prefix, data = pca_var_dt)))

cat("\n--- t-test: PC1 by effect direction ---\n")
print(t.test(PC1 ~ effect_dir, data = pca_var_dt[effect_dir %in% c("hot", "cold")]))


table(pca_var_dt$gene_prefix)

heatmap_dt <- clean_var_dt[, .N, by = .(gene_prefix, sig_class)]

# Heatmap of counts
p1 <- ggplot(heatmap_dt, aes(x = sig_class, y = gene_prefix, fill = N)) +
  geom_tile() +
  geom_text(aes(label = N), size = 3) +
  scale_fill_viridis_c() +
  labs(title = "Variant counts by gene and significance class", 
       x = "", y = "", fill = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Proportion within gene (more informative - controls for gene size)
heatmap_dt[, prop := N / sum(N), by = gene_prefix]

p2 <- ggplot(heatmap_dt, aes(x = sig_class, y = gene_prefix, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", prop)), size = 3, color="white") +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Proportion of variants by significance class (within gene)", 
       x = "", y = "", fill = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Chi-squared test for independence
contingency <- dcast(heatmap_dt, gene_prefix ~ sig_class, value.var = "N", fill = 0)
mat <- as.matrix(contingency[, -1])
rownames(mat) <- contingency$gene_prefix
cat("\nChi-squared test (gene vs sig_class):\n")
print(chisq.test(mat))

print(p1)
print(p2)

summary(heatmap_dt)
heatmap_dt

# ---- transition analysis ---- 
# Build transition-level table instead of tip-level
# Each row = one edge where a substitution occurred

cat("\n=== BUILDING TRANSITION-LEVEL TABLE ===\n")

edges <- anc_tree$edge
node_labels <- c(anc_tree$tip.label, anc_tree$node.label)

# Get branch lengths for each edge
edge_lengths <- anc_tree$edge.length

# Precompute depths at child node of each edge
bl_depths <- node.depth.edgelength(anc_tree)
topo_depths <- node.depth(anc_tree)

transition_results <- list()

sites_to_process <- site_info[!is.na(P_aa_only)]

for (i in seq_len(nrow(sites_to_process))) {
  site <- sites_to_process[i]
  gpos <- site$GlobalPos
  
  # Tip states
  tip_states <- supermat_mat[, gpos]
  names(tip_states) <- rownames(supermat_mat)
  
  # Internal node states
  node_states <- substr(anc_raw$seq, gpos, gpos)
  names(node_states) <- anc_raw$node
  
  all_states <- c(tip_states, node_states)
  
  # Parent/child for each edge
  parent_state <- all_states[node_labels[edges[, 1]]]
  child_state <- all_states[node_labels[edges[, 2]]]
  
  # Find edges with substitutions (excluding gaps)
  is_sub <- parent_state != child_state & 
    parent_state != "-" & child_state != "-" &
    !is.na(parent_state) & !is.na(child_state)
  
  if (sum(is_sub) == 0) next
  
  # Get per-residue effects
  eff_dt <- site$effects[[1]]
  
  # For each substitution edge
  sub_idx <- which(is_sub)
  for (j in sub_idx) {
    from_aa <- parent_state[j]
    to_aa <- child_state[j]
    child_node <- edges[j, 2]
    
    # Effect of derived state
    p_residue <- NA_real_
    residue_effect <- NA_real_
    if (!is.null(eff_dt) && nrow(eff_dt) > 0) {
      residue_row <- eff_dt[grepl(paste0(to_aa, "$"), Residue)]
      if (nrow(residue_row) == 1) {
        p_residue <- residue_row$P_value
        residue_effect <- residue_row$Effect
      }
    }
    
    transition_results[[length(transition_results) + 1]] <- data.table(
      Gene = site$Gene,
      Position = site$GenePos,
      GlobalPos = gpos,
      from_aa = from_aa,
      to_aa = to_aa,
      substitution = paste0(from_aa, "→", to_aa),
      edge_id = j,
      child_node = child_node,
      is_tip_edge = child_node <= Ntip(anc_tree),
      edge_length = edge_lengths[j],
      child_bl_depth = bl_depths[child_node],
      child_topo_depth = topo_depths[child_node],
      P_aa_only = site$P_aa_only,
      P_aa_with_pcs = site$P_aa_with_pcs,
      P_residue = p_residue,
      residue_effect = residue_effect,
      sig_class = site$sig_class
    )
  }
  
  if (i %% 500 == 0) cat("  Processed:", i, "/", nrow(sites_to_process), "\n")
}

trans_dt <- rbindlist(transition_results)
cat("Transitions:", nrow(trans_dt), "\n")
trans_dt$variant_sig_class <- trans_dt$sig_class
trans_dt$variant_sig_class[trans_dt$P_residue > 0.05] <- "not_sig"

# Basic checks
stopifnot(nrow(trans_dt) > 0)
table(trans_dt$sig_class)
table(trans_dt$variant_sig_class)
# 1. Are substitutions at sig sites happening on different parts of tree?
boxplot(child_bl_depth ~ variant_sig_class, trans_dt, 
        main = "Depth of substitutions by sig class")

# 2. Are sig sites substituting on shorter/longer branches?
#    (short branches = recent, long = old or fast-evolving)
boxplot(edge_length ~ variant_sig_class, trans_dt,
        main = "Branch length of substitution edges")

# 3. Substitution matrix by class - are certain AA changes enriched?
sub_counts <- trans_dt[, .N, by = .(variant_sig_class, substitution)]
sub_wide <- dcast(sub_counts, substitution ~ variant_sig_class, value.var = "N", fill = 0)

# Top substitutions in sig_both vs not_sig
trans_dt[, .N, by = .(variant_sig_class, from_aa, to_aa)][order(-N)][variant_sig_class == "sig_both"][1:20]

# 4. Ratio of substitutions per site (transition density)
site_trans_rate <- trans_dt[, .(
  n_subs = .N,
  mean_edge_len = mean(edge_length)
), by = .(Gene, Position, GlobalPos, variant_sig_class)]

boxplot(n_subs ~ variant_sig_class, site_trans_rate,
        main = "Substitutions per site by sig class")
boxplot(n_subs ~ Gene, site_trans_rate[variant_sig_class=="sig_both"],
        main = "Substitutions per gene in both_sig",
        las=0.45)
# Build transition matrices per sig_class and test for differences

library(data.table)

# Count transitions by class
trans_counts <- trans_dt[, .N, by = .(variant_sig_class, from_aa, to_aa)]

# Get all amino acids observed
all_aa <- sort(unique(c(trans_dt$from_aa, trans_dt$to_aa)))

# Function to build matrix from counts
build_trans_matrix <- function(counts_dt, aa_levels) {
  mat <- matrix(0, nrow = length(aa_levels), ncol = length(aa_levels),
                dimnames = list(from = aa_levels, to = aa_levels))
  for (i in seq_len(nrow(counts_dt))) {
    from <- counts_dt$from_aa[i]
    to <- counts_dt$to_aa[i]
    if (from %in% aa_levels && to %in% aa_levels) {
      mat[from, to] <- counts_dt$N[i]
    }
  }
  mat
}

# Build matrix for each class
classes <- c("not_sig", "sig_nocontrol", "sig_control", "sig_both")
trans_mats <- lapply(setNames(classes, classes), function(cl) {
  build_trans_matrix(trans_counts[variant_sig_class == cl], all_aa)
})

# Print totals
cat("Total transitions per class:\n")
sapply(trans_mats, sum)

# --- Test 1: Chi-squared on flattened matrices (pairwise) ---
cat("\n=== PAIRWISE CHI-SQUARED TESTS ===\n")

compare_matrices <- function(mat1, mat2, name1, name2) {
  # Flatten to vectors, combine as 2-row contingency table
  v1 <- as.vector(mat1)
  v2 <- as.vector(mat2)
  
  # Remove cells where both are zero
  keep <- (v1 + v2) > 0
  v1 <- v1[keep]
  v2 <- v2[keep]
  
  # Chi-squared test
  contingency <- rbind(v1, v2)
  rownames(contingency) <- c(name1, name2)
  
  test <- chisq.test(contingency, simulate.p.value = TRUE, B = 2000)
  list(
    comparison = paste(name1, "vs", name2),
    chi_sq = test$statistic,
    p_value = test$p.value
  )
}

# All pairwise comparisons
pairs <- combn(classes, 2, simplify = FALSE)
pairwise_results <- rbindlist(lapply(pairs, function(p) {
  compare_matrices(trans_mats[[p[1]]], trans_mats[[p[2]]], p[1], p[2])
}))
print(pairwise_results)

par(mfrow=c(2,2))

heatmap(trans_mats[[1]])
par(mfrow=c(1,1))

# --- Test 2: Which specific transitions are enriched? ---
cat("\n=== ENRICHED TRANSITIONS (sig_both vs not_sig) ===\n")

# Normalize to proportions
prop_mats <- lapply(trans_mats, function(m) m / sum(m))

# Difference: sig_both - not_sig
diff_mat <- prop_mats[["sig_both"]] - prop_mats[["not_sig"]]

# Find biggest differences
diff_dt <- as.data.table(as.table(diff_mat))
setnames(diff_dt, c("from_aa", "to_aa", "diff"))
diff_dt <- diff_dt[order(-abs(diff))]

cat("Top enriched in sig_both:\n")
print(head(diff_dt[diff > 0], 15))

cat("\nTop depleted in sig_both:\n")
print(head(diff_dt[diff < 0], 15))

# --- Test 3: Log-odds ratio for specific transitions ---
cat("\n=== LOG-ODDS RATIOS ===\n")

calc_log_odds <- function(mat1, mat2, pseudocount = 0.5) {
  # Add pseudocount, normalize
  p1 <- (mat1 + pseudocount) / sum(mat1 + pseudocount)
  p2 <- (mat2 + pseudocount) / sum(mat2 + pseudocount)
  log2(p1 / p2)
}

lor_mat <- calc_log_odds(trans_mats[["sig_both"]], trans_mats[["not_sig"]])

lor_dt <- as.data.table(as.table(lor_mat))
setnames(lor_dt, c("from_aa", "to_aa", "log2_OR"))

# Add counts for context
lor_dt <- merge(lor_dt, 
                trans_counts[sig_class == "sig_both", .(from_aa, to_aa, N_sigboth = N)],
                by = c("from_aa", "to_aa"), all.x = TRUE)
lor_dt[is.na(N_sigboth), N_sigboth := 0]

cat("Transitions most enriched in sig_both (log2 OR):\n")
print(lor_dt[N_sigboth >= 5][order(-log2_OR)][1:20])

cat("\nTransitions most depleted in sig_both:\n")
print(lor_dt[N_sigboth >= 5][order(log2_OR)][1:20])

# --- Visualization: Heatmap of log-odds ---

# Filter to common amino acids for cleaner heatmap
common_aa <- names(which(rowSums(trans_mats[["not_sig"]]) > 20))
lor_plot <- lor_dt[from_aa %in% common_aa & to_aa %in% common_aa]

ggplot(lor_plot, aes(x = to_aa, y = from_aa, fill = log2_OR)) +
  geom_tile() +
  geom_text(aes(label = N_sigboth), size = 2.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
  labs(title = "Log2 odds ratio: sig_both vs not_sig",
       subtitle = "Numbers = count in sig_both",
       x = "To (derived)", y = "From (ancestral)", fill = "log2 OR") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0))

library(ggplot2)
library(data.table)

# Normalize matrices to proportions
prop_mats <- lapply(trans_mats, function(m) m / sum(m))

# Convert to long format for ggplot
mat_to_dt <- function(mat, class_name) {
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("from_aa", "to_aa", "prop"))
  dt[, sig_class := class_name]
  dt
}

plot_dt <- rbindlist(lapply(names(prop_mats), function(cl) {
  mat_to_dt(prop_mats[[cl]], cl)
}))

# Order sig_class factor
plot_dt[, sig_class := factor(sig_class, levels = c("not_sig", "sig_nocontrol", 
                                                    "sig_control", "sig_both"))]

# Filter to common amino acids
common_aa <- names(which(rowSums(trans_mats[["not_sig"]]) > 20))
plot_dt <- plot_dt[from_aa %in% common_aa & to_aa %in% common_aa]

# Order amino acids by property (hydrophobic -> polar -> charged)
aa_order <- c("A", "V", "L", "I", "M", "F", "W", "P", "G", 
              "S", "T", "N", "Q", "Y", "C", 
              "K", "R", "H", "D", "E")
aa_order <- intersect(aa_order, common_aa)  # keep only those present

plot_dt[, from_aa := factor(from_aa, levels = aa_order)]
plot_dt[, to_aa := factor(to_aa, levels = aa_order)]
sum(is.na(plot_dt))
plot_dt <- na.omit(plot_dt)

ggplot(plot_dt, aes(x = to_aa, y = from_aa, fill = prop)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "sqrt",
                       labels = scales::percent) +
  facet_wrap(~sig_class, ncol = 2) +
  labs(title = "Transition matrices by significance class",
       x = "To (derived)", y = "From (ancestral)", fill = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"))

# Classify substitutions by biochemical impact
conservative <- c("V-I", "I-V", "L-M", "M-L", "L-F", "F-L", "K-R", "R-K", "E-D", "D-E", "S-T", "T-S")
trans_dt[, sub_pair := paste(from_aa, to_aa, sep = "-")]
trans_dt[, is_conservative := sub_pair %in% conservative]

# Test
table(trans_dt$sig_class, trans_dt$is_conservative)
chisq.test(table(trans_dt$sig_class, trans_dt$is_conservative))

# Proportion conservative by class
trans_dt[, .(prop_conservative = mean(is_conservative), n = .N), by = sig_class]

# Log-odds ratio vs not_sig for each class
calc_log_odds <- function(mat1, mat2, pseudocount = 0.5) {
  p1 <- (mat1 + pseudocount) / sum(mat1 + pseudocount)
  p2 <- (mat2 + pseudocount) / sum(mat2 + pseudocount)
  log2(p1 / p2)
}

lor_mats <- lapply(c("sig_nocontrol", "sig_control", "sig_both"), function(cl) {
  calc_log_odds(trans_mats[[cl]], trans_mats[["not_sig"]])
})
names(lor_mats) <- c("sig_nocontrol", "sig_control", "sig_both")

# Convert to long format
mat_to_dt <- function(mat, class_name) {
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("from_aa", "to_aa", "log2_OR"))
  dt[, sig_class := class_name]
  dt
}

lor_mats_clean <- lapply(lor_mats, function(m) {
  diag(m) <- NA
  m
})

# Rebuild plot_dt
plot_dt <- rbindlist(lapply(names(lor_mats_clean), function(cl) {
  mat_to_dt(lor_mats_clean[[cl]], cl)
}))



# Add counts from each class for filtering
for (cl in names(lor_mats)) {
  counts <- as.data.table(as.table(trans_mats[[cl]]))
  setnames(counts, c("from_aa", "to_aa", "N"))
  plot_dt[sig_class == cl, N := counts$N[match(paste(from_aa, to_aa), 
                                               paste(counts$from_aa, counts$to_aa))]]
}

# Filter to amino acids with sufficient data
common_aa <- names(which(rowSums(trans_mats[["not_sig"]]) > 50))

# Order by biochemical property
aa_order <- c("A", "V", "L", "I", "M", "F", "W", "P", "G", 
              "S", "T", "N", "Q", "Y", "C", 
              "K", "R", "H", "D", "E")
aa_order <- intersect(aa_order, common_aa)

plot_dt <- plot_dt[from_aa %in% aa_order & to_aa %in% aa_order]
plot_dt[, from_aa := factor(from_aa, levels = aa_order)]
plot_dt[, to_aa := factor(to_aa, levels = aa_order)]
plot_dt[, sig_class := factor(sig_class, levels = c("sig_nocontrol", "sig_control", "sig_both"))]

# Main plot
ggplot(plot_dt, aes(x = to_aa, y = from_aa, fill = log2_OR)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = "log2 OR\nvs not_sig") +
  facet_wrap(~sig_class, ncol = 3) +
  labs(title = "Transition enrichment relative to not_sig",
       subtitle = "Red = enriched, Blue = depleted",
       x = "To (derived)", y = "From (ancestral)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"))

diag(trans_mats[["not_sig"]])
diag(trans_mats[["sig_both"]])
trans_dt[from_aa == to_aa, .N, by = sig_class]

# Calculate log-odds with significance test per cell
calc_lor_with_sig <- function(mat1, mat2, pseudocount = 0.5) {
  # Counts
  n1 <- sum(mat1)
  n2 <- sum(mat2)
  
  results <- list()
  aas <- rownames(mat1)
  
  for (from in aas) {
    for (to in aas) {
      if (from == to) next
      
      a <- mat1[from, to]  # count in class 1
      b <- mat2[from, to]  # count in class 2
      
      # Log odds ratio with pseudocount
      p1 <- (a + pseudocount) / (n1 + pseudocount * length(aas)^2)
      p2 <- (b + pseudocount) / (n2 + pseudocount * length(aas)^2)
      lor <- log2(p1 / p2)
      
      # Fisher's exact test (or chi-squared for speed)
      # 2x2 table: this_transition vs all_others, class1 vs class2
      cont <- matrix(c(a, n1 - a, b, n2 - b), nrow = 2)
      pval <- tryCatch(
        fisher.test(cont)$p.value,
        error = function(e) chisq.test(cont)$p.value
      )
      
      results[[length(results) + 1]] <- data.table(
        from_aa = from, to_aa = to, 
        n_test = a, n_ref = b,
        log2_OR = lor, p_value = pval
      )
    }
  }
  rbindlist(results)
}

# Calculate for each class vs not_sig
sig_tests <- lapply(c("sig_nocontrol", "sig_control", "sig_both"), function(cl) {
  dt <- calc_lor_with_sig(trans_mats[[cl]], trans_mats[["not_sig"]])
  dt[, sig_class := cl]
  dt
})
sig_dt <- rbindlist(sig_tests)

# FDR correction
sig_dt[, p_adj := p.adjust(p_value, method = "BH")]

# Filter: significant and meaningful effect size
sig_dt[, show := p_adj < 0.01 & abs(log2_OR) > 0.5]

cat("Significant transitions per class:\n")
sig_dt[show == TRUE, .N, by = sig_class]

# Filter to common AAs
aa_order <- c("A", "V", "L", "I", "M", "F", "W", "P", "G", 
              "S", "T", "N", "Q", "Y", "C", 
              "K", "R", "H", "D", "E")
plot_dt <- sig_dt[from_aa %in% aa_order & to_aa %in% aa_order]
plot_dt[, from_aa := factor(from_aa, levels = aa_order)]
plot_dt[, to_aa := factor(to_aa, levels = aa_order)]
plot_dt[, sig_class := factor(sig_class, levels = c("sig_nocontrol", "sig_control", "sig_both"))]

# Set non-significant to NA (will appear grey)
plot_dt[show == FALSE, log2_OR := NA]

ggplot(plot_dt, aes(x = to_aa, y = from_aa, fill = log2_OR)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = "log2 OR", na.value = "grey90") +
  facet_wrap(~sig_class, ncol = 3) +
  labs(title = "Transition enrichment relative to not_sig",
       subtitle = "Colored = FDR < 0.01 & |log2 OR| > 0.5; Grey = NS",
       x = "To (derived)", y = "From (ancestral)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"))


library(Peptides)

# Count derived alleles by class
to_counts <- trans_dt[, .N, by = .(sig_class, to_aa)]
to_counts[, prop := N / sum(N), by = sig_class]

# Wide format for comparison
to_wide <- dcast(to_counts, to_aa ~ sig_class, value.var = "prop", fill = 0)

# Log-odds vs not_sig
to_wide[, lor_sigboth := log2((sig_both + 0.001) / (not_sig + 0.001))]
to_wide[, lor_sigcontrol := log2((sig_control + 0.001) / (not_sig + 0.001))]
to_wide[, lor_signocontrol := log2((sig_nocontrol + 0.001) / (not_sig + 0.001))]

# Sort by enrichment in sig_both
setorder(to_wide, -lor_sigboth)
print(to_wide[to_aa %in% aa_order, .(to_aa, not_sig, sig_both, lor_sigboth)])

# Simple bar plot
to_plot <- melt(to_counts[to_aa %in% aa_order], 
                id.vars = c("sig_class", "to_aa"), 
                measure.vars = "prop")

ggplot(to_plot, aes(x = to_aa, y = value, fill = sig_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(not_sig = "grey70", sig_nocontrol = "steelblue",
                               sig_control = "darkorange", sig_both = "firebrick")) +
  labs(x = "Derived amino acid", y = "Proportion", 
       title = "Derived allele preference by significance class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Chi-squared test: is derived AA distribution different across classes?
cont_table <- dcast(to_counts[to_aa %in% aa_order], to_aa ~ sig_class, value.var = "N", fill = 0)
mat <- as.matrix(cont_table[, -1])
rownames(mat) <- cont_table$to_aa
chisq.test(mat)

# ---- aa propoertie preference ----

library(Peptides)
data(AAdata)


# Get AA properties
hydro <- AAdata$Hydrophobicity$Tanford
kidera <- AAdata$kideraFactors

# Build property table for each AA
aa_props <- data.table(
  aa = names(hydro),
  hydrophobicity = as.numeric(hydro)
)

# Add Kidera factors
for (kf in names(kidera)) {
  aa_props[, (kf) := as.numeric(kidera[[kf]][aa])]
}

# --- Q1: Is there a preference for hydrophobic derived alleles? ---
cat("\n=== HYDROPHOBICITY OF DERIVED ALLELES ===\n")

# Add hydrophobicity to trans_dt
trans_dt[, to_hydro := hydro[to_aa]]
trans_dt[, from_hydro := hydro[from_aa]]
trans_dt[, delta_hydro := to_hydro - from_hydro]

# Mean hydrophobicity of derived allele by class
hydro_summary <- trans_dt[, .(
  mean_to_hydro = mean(to_hydro, na.rm = TRUE),
  mean_from_hydro = mean(from_hydro, na.rm = TRUE),
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = variant_sig_class]
print(hydro_summary)

# Test
cat("\nKruskal-Wallis: derived hydrophobicity ~ sig_class\n")
kruskal.test(to_hydro ~ variant_sig_class, data = trans_dt)

cat("\nKruskal-Wallis: delta hydrophobicity ~ sig_class\n")
kruskal.test(delta_hydro ~ variant_sig_class, data = trans_dt)

# Visualize
boxplot(to_hydro ~ variant_sig_class, trans_dt, 
        col = c("grey70", "steelblue", "darkorange", "firebrick")[
          match(levels(factor(trans_dt$variant_sig_class)), 
                c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))],
        main = "Hydrophobicity of derived allele", ylab = "Tanford hydrophobicity")

boxplot(delta_hydro ~ variant_sig_class, trans_dt,
        col = c("grey70", "steelblue", "darkorange", "firebrick")[
          match(levels(factor(trans_dt$variant_sig_class)), 
                c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))],
        main = "Change in hydrophobicity (to - from)", ylab = "Delta hydrophobicity")
abline(h = 0, lty = 2)

# --- Q2: Non-conservative changes via Kidera distance ---
cat("\n=== KIDERA FACTOR DISTANCE (BIOCHEMICAL DISSIMILARITY) ===\n")

# Euclidean distance in Kidera space between from and to
kf_names <- paste0("KF", 1:10)

# Build matrix of Kidera factors
kf_mat <- as.matrix(aa_props[, ..kf_names])
rownames(kf_mat) <- aa_props$aa

# Calculate pairwise Kidera distance for each transition
calc_kidera_dist <- function(from, to, kf_mat) {
  if (is.na(from) || is.na(to) || !from %in% rownames(kf_mat) || !to %in% rownames(kf_mat)) {
    return(NA_real_)
  }
  sqrt(sum((kf_mat[from, ] - kf_mat[to, ])^2))
}

trans_dt[, kidera_dist := mapply(calc_kidera_dist, from_aa, to_aa, MoreArgs = list(kf_mat = kf_mat))]

# Summary by class
kidera_summary <- trans_dt[, .(
  mean_kidera_dist = mean(kidera_dist, na.rm = TRUE),
  median_kidera_dist = median(kidera_dist, na.rm = TRUE),
  n = .N
), by = variant_sig_class]
print(kidera_summary)

# Test
cat("\nKruskal-Wallis: Kidera distance ~ sig_class\n")
kruskal.test(kidera_dist ~ variant_sig_class, data = trans_dt)

# Pairwise vs not_sig
cat("\nWilcoxon: sig_both vs not_sig\n")
wilcox.test(kidera_dist ~ variant_sig_class, 
            data = trans_dt[variant_sig_class %in% c("sig_both", "not_sig")])

# Visualize
cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

trans_dt[, variant_sig_class := factor(variant_sig_class, levels = names(cols))]

ggplot(trans_dt[!is.na(kidera_dist)], aes(x = variant_sig_class, y = kidera_dist, fill = variant_sig_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "Kidera distance (biochemical dissimilarity)",
       title = "Substitution magnitude by significance class",
       subtitle = "Higher = more non-conservative") +
  theme_bw() +
  theme(legend.position = "none")

conservative_pairs <- c("V-I", "I-V", "L-M", "M-L", "L-F", "F-L", "K-R", "R-K", "E-D", "D-E", "S-T", "T-S")

trans_dt[, .(
  mean_kidera = mean(kidera_dist, na.rm = TRUE),
  n = .N
), by = is_conservative]

# Bin Kidera distance and look at distribution by class
trans_dt[, kidera_bin := cut(kidera_dist, 
                             breaks = c(0, 3, 3.5, 4, 4.5, 5, Inf),
                             labels = c("<3", "3-3.5", "3.5-4", "4-4.5", "4.5-5", ">5"))]

# Proportion in each bin by class
kidera_dist_table <- trans_dt[, .N, by = .(variant_sig_class, kidera_bin)]
kidera_dist_table[, prop := N / sum(N), by = variant_sig_class]

dcast(kidera_dist_table, kidera_bin ~ variant_sig_class, value.var = "prop")

# Pairwise Kidera distance heatmap between all amino acids

# Calculate all pairwise distances
aa_order <- c("A", "V", "L", "I", "M", "F", "W", "P", "G", 
              "S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E")

# Build pairwise distance matrix
kidera_pairwise <- expand.grid(from_aa = aa_order, to_aa = aa_order, stringsAsFactors = FALSE)
kidera_pairwise <- as.data.table(kidera_pairwise)

kidera_pairwise[, kidera_dist := mapply(calc_kidera_dist, from_aa, to_aa, 
                                        MoreArgs = list(kf_mat = kf_mat))]

# Set diagonal to NA
kidera_pairwise[from_aa == to_aa, kidera_dist := NA]

kidera_pairwise[, from_aa := factor(from_aa, levels = aa_order)]
kidera_pairwise[, to_aa := factor(to_aa, levels = aa_order)]

ggplot(kidera_pairwise, aes(x = to_aa, y = from_aa, fill = kidera_dist)) +
  geom_tile() +
  geom_text(aes(label = round(kidera_dist, 1)), size = 2.5) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "white",
                       name = "Kidera\ndistance") +
  labs(title = "Pairwise Kidera distance between amino acids",
       subtitle = "Lower = more biochemically similar (conservative substitution)",
       x = "To", y = "From") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 9),
        axis.text.y = element_text(size = 9))

# Extract upper triangle of pairwise distances
kidera_mat <- dcast(kidera_pairwise, from_aa ~ to_aa, value.var = "kidera_dist")
kidera_mat_m <- as.matrix(kidera_mat[, -1])
rownames(kidera_mat_m) <- kidera_mat$from_aa

# Get upper triangle values (excluding diagonal)
upper_vals <- kidera_mat_m[upper.tri(kidera_mat_m)]

# Histogram
hist(upper_vals, breaks = 30, col = "steelblue", border = "white",
     main = "Distribution of pairwise Kidera distances",
     xlab = "Kidera distance", ylab = "Frequency")
abline(v = mean(upper_vals, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
abline(v = median(upper_vals, na.rm = TRUE), col = "darkred", lty = 1, lwd = 2)
legend("topright", c(paste("Mean:", round(mean(upper_vals, na.rm = TRUE), 2)),
                     paste("Median:", round(median(upper_vals, na.rm = TRUE), 2))),
       col = c("red", "darkred"), lty = c(2, 1), lwd = 2, bty = "n")

# Add reference lines for conservative threshold
abline(v = 3.5, col = "orange", lty = 2, lwd = 2)
text(3.5, par("usr")[4] * 0.9, "Conservative\nthreshold", pos = 2, cex = 0.8)



# ---- stratify by hot.cold ----
# Stratify by effect direction
trans_dt[, effect_dir := fcase(
  residue_effect > 0, "hot",
  residue_effect < 0, "cold",
  default = NA_character_
)]

cat("Effect direction counts:\n")
table(trans_dt$effect_dir, trans_dt$variant_sig_class, useNA = "ifany")

# Filter to transitions with known effect direction
trans_eff <- trans_dt[!is.na(effect_dir)]
cat("\nTransitions with effect direction:", nrow(trans_eff), "\n")

# --- HYDROPHOBICITY BY EFFECT DIRECTION ---
cat("\n=== HYDROPHOBICITY BY EFFECT DIRECTION ===\n")

hydro_by_effect <- trans_eff[, .(
  mean_to_hydro = mean(to_hydro, na.rm = TRUE),
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = .(variant_sig_class, effect_dir)]
print(dcast(hydro_by_effect, variant_sig_class ~ effect_dir, value.var = "mean_to_hydro"))
print(dcast(hydro_by_effect, variant_sig_class ~ effect_dir, value.var = "mean_delta_hydro"))

# Test interaction
cat("\nTwo-way ANOVA: to_hydro ~ sig_class * effect_dir\n")
summary(aov(mean_delta_ ~ variant_sig_class * effect_dir, data = trans_eff))

# Visualize
cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

ggplot(trans_eff, aes(x = variant_sig_class, y = delta_hydro, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = cols) +
  facet_wrap(~effect_dir) +
  labs(x = "", y = "Delta hydrophobicity",
       title = "Hydrophobicity Change of transition by effect direction") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

all_comparisons <- combn(class_order, 2, simplify = FALSE)

ggplot(trans_eff, aes(x = variant_sig_class, y = delta_hydro, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = cols) +
  geom_signif(
    comparisons = all_comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.08,
    size = 0.3,
    textsize = 2.5
  ) +
  facet_wrap(~effect_dir) +
  labs(x = "", y = "Delta hydrophobicity",
       title = "Hydrophobicity change of transition by effect direction") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


# --- KIDERA DISTANCE BY EFFECT DIRECTION ---
cat("\n=== KIDERA DISTANCE BY EFFECT DIRECTION ===\n")

kidera_by_effect <- trans_eff[, .(
  mean_kidera = mean(kidera_dist, na.rm = TRUE),
  n = .N
), by = .(variant_sig_class, effect_dir)]
print(dcast(kidera_by_effect, variant_sig_class ~ effect_dir, value.var = "mean_kidera"))

cat("\nTwo-way ANOVA: kidera_dist ~ sig_class * effect_dir\n")
summary(aov(kidera_dist ~ sig_class * effect_dir, data = trans_eff))


library(ggplot2)
library(ggsignif)

cols <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")

class_order <- c("not_sig", "sig_nocontrol", "sig_control", "sig_both")
trans_eff[, variant_sig_class := factor(variant_sig_class, levels = class_order)]

ggplot(trans_eff[!is.na(kidera_dist)], aes(x = variant_sig_class, y = kidera_dist, fill = variant_sig_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  facet_wrap(~effect_dir) +
  geom_signif(
    comparisons = list(c("not_sig", "sig_both")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = max(trans_eff$kidera_dist, na.rm = TRUE) * 0.95
  ) +
  labs(x = "", y = "Kidera distance (biochemical dissimilarity)",
       title = "Substitution magnitude by significance class and effect direction",
       subtitle = "Higher = more non-conservative") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(trans_dt[variant_sig_class %in% c("sig_both", "not_sig") & !is.na(kidera_dist)], 
       aes(x = kidera_dist, fill = variant_sig_class)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(not_sig = "grey70", sig_both = "firebrick")) +
  labs(x = "Kidera distance", y = "Density",
       title = "Kidera Distance Density: sig_both vs not_sig",
       fill = "") +
  theme_bw() +
  theme(legend.position = "top")

ggplot(trans_eff[variant_sig_class == "sig_both" & !is.na(kidera_dist)], 
       aes(x = effect_dir, y = kidera_dist, fill = effect_dir)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c(cold = "blue", hot = "red")) +
  geom_signif(
    comparisons = list(c("cold", "hot")),
    test = "wilcox.test",
    map_signif_level = TRUE
  ) +
  labs(x = "", y = "Kidera distance (biochemical dissimilarity)",
       title = "sig_both: Substitution magnitude by thermal direction",
       subtitle = "Higher = more non-conservative") +
  theme_bw() +
  theme(legend.position = "none")

# --- CONSERVATIVE SUBSTITUTIONS BY EFFECT DIRECTION ---
cat("\n=== CONSERVATIVE SUBSTITUTIONS BY EFFECT DIRECTION ===\n")

cons_by_effect <- trans_eff[, .(
  prop_conservative = mean(is_conservative),
  n = .N
), by = .(variant_sig_class, effect_dir)]
print(dcast(cons_by_effect, variant_sig_class ~ effect_dir, value.var = "prop_conservative"))

# Chi-squared for interaction
cat("\nChi-squared: conservative ~ sig_class, stratified by effect_dir\n")
cat("\nHOT:\n")
print(chisq.test(table(trans_eff[effect_dir == "hot", .(variant_sig_class, is_conservative)])))
cat("\nCOLD:\n
")
print(chisq.test(table(trans_eff[effect_dir == "cold", .(variant_sig_class, is_conservative)])))

# --- DERIVED AA PREFERENCE BY EFFECT DIRECTION ---
cat("\n=== DERIVED AA ENRICHMENT: HOT vs COLD ===\n")

# Within sig_both, compare hot vs cold
to_by_effect <- trans_eff[variant_sig_class == "sig_both", .N, by = .(effect_dir, to_aa)]
to_by_effect[, prop := N / sum(N), by = effect_dir]

to_wide <- dcast(to_by_effect, to_aa ~ effect_dir, value.var = "prop", fill = 0)
to_wide[, lor_hot_vs_cold := log2((hot + 0.001) / (cold + 0.001))]
setorder(to_wide, -lor_hot_vs_cold)

cat("\nDerived AA preference in sig_both (hot vs cold):\n")
print(to_wide)
aa_order <- Peptides::aaList()
# Visualize
ggplot(to_by_effect[to_aa %in% aa_order], aes(x = to_aa, y = prop, fill = effect_dir)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(cold = "blue", hot = "red")) +
  labs(x = "Derived amino acid", y = "Proportion",
       title = "Derived AA preference in sig_both: hot vs cold") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- KIDERA BIN DISTRIBUTION BY EFFECT ---
cat("\n=== KIDERA BIN DISTRIBUTION BY EFFECT ===\n")

kidera_bin_effect <- trans_eff[variant_sig_class == "sig_both", .N, by = .(effect_dir, kidera_bin)]
kidera_bin_effect[, prop := N / sum(N), by = effect_dir]

print(dcast(kidera_bin_effect, kidera_bin ~ effect_dir, value.var = "prop"))

#stratify by hot/cold 
ggplot(trans_eff[variant_sig_class == "sig_both"], aes(x = effect_dir, y = delta_hydro, fill = effect_dir)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c(cold = "blue", hot = "red")) +
  labs(title = "sig_both: Hydrophobicity change by thermal direction",
       y = "Delta hydrophobicity", x = "") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(trans_eff[variant_sig_class == "sig_both"], aes(x = effect_dir, y = delta_hydro, fill = effect_dir)) +
  geom_violin() +
  geom_boxplot(width=0.15) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c(cold = "steelblue", hot = "tomato")) +
  geom_signif(
    comparisons = list(c("cold", "hot")),
    test = "wilcox.test",
    map_signif_level = TRUE
  ) +
  labs(title = "sig_both: Hydrophobicity change by thermal direction",
       y = "Delta hydrophobicity", x = "") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(trans_eff[variant_sig_class == "sig_both"], aes(x = delta_hydro, fill = effect_dir)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_fill_manual(values = c(cold = "blue", hot = "red")) +
  labs(title = "sig_both: Hydrophobicity change by thermal direction",
       x = "Delta hydrophobicity", y = "Count", fill = "") +
  theme_bw() +
  theme(legend.position = "top")

# Transition counts by effect direction within sig_both
trans_by_effect <- trans_eff[variant_sig_class == "sig_both", .N, by = .(effect_dir, from_aa, to_aa)]
trans_by_effect[, prop := N / sum(N), by = effect_dir]
table(trans_by_effect$N)

# Wide format for comparison
trans_wide <- dcast(trans_by_effect, from_aa + to_aa ~ effect_dir, 
                    value.var = c("N", "prop"), fill = 0)

# Log odds ratio hot vs cold
trans_wide[, lor := log2((prop_hot + 0.0001) / (prop_cold + 0.0001))]

# Add delta hydrophobicity for each transition
hydro <- AAdata$Hydrophobicity$Tanford
trans_wide[, delta_hydro := hydro[to_aa] - hydro[from_aa]]

# Total count for filtering
trans_wide[, N_total := N_hot + N_cold]

# Filter to reasonably common transitions
trans_sig <- trans_wide[N_total >= 5]


head(trans_sig)
cold_vals <- rep(trans_sig$delta_hydro, trans_sig$N_cold)
hot_vals <- rep(trans_sig$delta_hydro, trans_sig$N_hot)

stopifnot(length(cold_vals) == sum(trans_sig$N_cold))
stopifnot(length(hot_vals) == sum(trans_sig$N_hot))

# Compare distributions
wilcox.test(cold_vals, hot_vals)
t.test(cold_vals, hot_vals)

hist(cold_vals)
hist(hot_vals)

boxplot(
  list(Cold = cold_vals, Hot = hot_vals),
  ylab = "delta_hydro",
  col = c("steelblue", "tomato")
)


cat("=== TRANSITIONS ENRICHED IN HOT (more hydrophobic?) ===\n")
print(trans_sig[order(-lor)][1:20, .(from_aa, to_aa, N_cold, N_hot, lor, delta_hydro)])

cat("\n=== TRANSITIONS ENRICHED IN COLD (more hydrophilic?) ===\n")
print(trans_sig[order(lor)][1:20, .(from_aa, to_aa, N_cold, N_hot, lor, delta_hydro)])

# Visualize: lor vs delta_hydro (expect positive correlation)
ggplot(trans_sig, aes(x = delta_hydro, y = lor, size = N_total)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = trans_sig[abs(lor) > 1 | N_total > 200], 
            aes(label = paste0(from_aa, "→", to_aa)), 
            size = 3, nudge_y = 0.1) +
  labs(x = "Delta hydrophobicity (to - from)",
       y = "Log2 OR (hot / cold)",
       title = "Transitions driving thermal adaptation",
       subtitle = "Upper right = hot-enriched & hydrophobic gain; Lower left = cold-enriched & hydrophilic gain",
       size = "N transitions") +
  theme_bw()

# Correlation test
cat("\nCorrelation: lor vs delta_hydro\n")
cor.test(trans_sig$lor, trans_sig$delta_hydro, method="spearman")

# Heatmap of log-odds by transition (filtered)
aa_order <- c("A", "V", "L", "I", "M", "F", "W", "P", "G", 
              "S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E")

plot_trans <- trans_sig[from_aa %in% aa_order & to_aa %in% aa_order]
plot_trans[, from_aa := factor(from_aa, levels = aa_order)]
plot_trans[, to_aa := factor(to_aa, levels = aa_order)]

ggplot(plot_trans, aes(x = to_aa, y = from_aa, fill = lor)) +
  geom_tile() +
  geom_text(aes(label = N_total), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = "log2 OR\nhot/cold") +
  labs(title = "Transition preference: hot vs cold adaptation",
       subtitle = "Red = enriched in hot; Blue = enriched in cold; Numbers = total count",
       x = "To (derived)", y = "From (ancestral)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add all Kidera factor deltas to transitions
kidera <- AAdata$kideraFactors

for (kf in names(kidera)) {
  kf_vec <- kidera[[kf]]
  trans_sig[, paste0("delta_", kf) := kf_vec[to_aa] - kf_vec[from_aa]]
}

hydro <- AAdata$Hydrophobicity$Tanford
trans_sig[, delta_hydro := hydro[to_aa] - hydro[from_aa]]

sapply(trans_sig[, .(delta_hydro, delta_KF1, delta_KF2, lor)], class)


# Correlate lor with each property
properties <- c("delta_hydro", paste0("delta_KF", 1:10))

cor_results <- rbindlist(lapply(properties, function(p) {
  x <- as.numeric(trans_sig[[p]])
  y <- trans_sig$lor
  ct <- cor.test(x, y, use = "complete.obs")
  data.table(property = p, r = ct$estimate, p_value = ct$p.value)
}))


cor_dt <- data.table(
  property = properties,
  r = cor_results["r.cor", ],
  p_value = cor_results["p", ]
)


cor_results[, p_adj := p.adjust(p_value, method = "BH")]
setorder(cor_results, -r)


cat("Correlation of transition enrichment (hot/cold) with biochemical properties:\n")
print(cor_results)

# Kidera factor key:
# KF1: Helix/bend preference
# KF2: Side-chain size
# KF3: Extended structure preference
# KF4: Hydrophobicity (Kidera version)
# KF5: Double-bend preference
# KF6: Partial specific volume
# KF7: Flat extended preference
# KF8: Occurrence in alpha region
# KF9: pK-C
# KF10: Surrounding hydrophobicity

# Multiple regression: which properties jointly predict hot vs cold?
trans_sig_complete <- trans_sig[complete.cases(trans_sig[, ..properties])]

fit <- lm(lor ~ delta_hydro + delta_KF1 + delta_KF2 + delta_KF3 + delta_KF4 + 
            delta_KF5 + delta_KF6 + delta_KF7 + delta_KF8 + delta_KF9 + delta_KF10,
          data = trans_sig_complete, weights = N_total)

cat("\n=== MULTIPLE REGRESSION: lor ~ biochemical properties ===\n")
summary(fit)

# Visualize the top 2
par(mfrow = c(1, 2))

plot(trans_sig$delta_hydro, trans_sig$lor, 
     xlab = "Delta hydro", ylab = "log2 OR (hot/cold)",
     main = "delta hydro vs hot/cold enrichment", pch = 16, col = rgb(0,0,0,0.5))
abline(lm(lor ~ delta_hydro, trans_sig), col = "red")
abline(h = 0, v = 0, lty = 2)

plot(trans_sig$delta_KF4, trans_sig$lor,
     xlab = "Delta KF4 (kidera hydro)", ylab = "log2 OR (hot/cold)",
     main = "delta kidera hydro vs hot/cold enrichment", pch = 16, col = rgb(0,0,0,0.5))
abline(lm(lor ~ delta_KF4, trans_sig), col = "red")
abline(h = 0, v = 0, lty = 2)

par(mfrow = c(1, 1))

trans_sig

# High KF4 gain + low KF3 gain = most "hot-like"
trans_sig[, hot_score := delta_KF4 - delta_KF3]
trans_sig[order(-hot_score)][1:10, .(from_aa, to_aa, lor, delta_KF4, delta_KF3, N_total)]
trans_sig[order(hot_score)][1:10, .(from_aa, to_aa, lor, delta_KF4, delta_KF3, N_total)]

trans_sig[from_aa == "P" | to_aa == "P", 
          .(from_aa, to_aa, lor, delta_hydro, delta_KF4, N_total)][order(-lor)]

# Add proline indicator to model
trans_sig[, loses_P := as.integer(from_aa == "P")]
trans_sig[, gains_P := as.integer(to_aa == "P")]

fit2 <- lm(lor ~ delta_hydro + delta_KF3 + delta_KF4 + loses_P + gains_P, 
           data = trans_sig, weights = N_total)
summary(fit2)

# ---- PART 2!!! ----

#!/usr/bin/env Rscript
# Phenotypic Ancestral State Reconstruction Analysis
# Compares trait-based transition classification with sequence-based GWAS classification
# Tests whether hydrophobicity patterns differ when transitions are classified by phenotype shifts

library(data.table)
library(ape)
library(phytools)
library(Peptides)

# PART 1: LOAD REQUIRED DATA

cat("\n=== LOADING DATA ===\n")

# Tree and ancestral states
anc_tree <- read.tree("raxml_input/aa_tree_ASR.raxml.ancestralTree")
anc_raw <- fread("raxml_input/aa_tree_ASR.raxml.ancestralStates", 
                 header = FALSE, sep = "\t", col.names = c("node", "seq"))
partitionMap <- readRDS("raxml_input/partitionMap.rds")
supermat <- readAAStringSet("raxml_input/superaa_collapsed.fasta")
names(supermat) <- sub("\\|.*", "", names(supermat))
supermat_mat <- as.matrix(supermat)

# Phenotype data
pheno_data <- read_parquet("data/processed_data.parquet")
setDT(pheno_data)

stopifnot(nrow(anc_raw) > 0)
stopifnot(inherits(anc_tree, "phylo"))

# Identify temperature column
temp_col <- "pheno_wc2.1_2.5m_bio_8_p50"
stopifnot(!is.na(temp_col))
cat("Using temperature column:", temp_col, "\n")
hist(pheno_data$pheno_wc2.1_2.5m_bio_8_p50)
# Match phenotype to tree tips
pheno_data <- pheno_data[ID %in% anc_tree$tip.label]
cat("Tips with phenotype data:", nrow(pheno_data), "/", Ntip(anc_tree), "\n")

# PART 2: PHENOTYPIC ASR - RECONSTRUCT ANCESTRAL TEMPERATURES

cat("\n=== PHENOTYPIC ASR ===\n")

# Prepare tip values
tip_temps <- setNames(pheno_data[[temp_col]], pheno_data$ID)
tip_temps <- tip_temps[anc_tree$tip.label]

# Remove tips with missing phenotype
valid_tips <- names(tip_temps)[!is.na(tip_temps)]
cat("Tips with valid temperature:", length(valid_tips), "\n")

# Prune tree to tips with phenotype data
pruned_tree <- keep.tip(anc_tree, valid_tips)
tip_temps <- tip_temps[valid_tips]

stopifnot(length(tip_temps) == Ntip(pruned_tree))
stopifnot(all(names(tip_temps) == pruned_tree$tip.label))

# Perform continuous ASR using fastAnc (ML estimation)
cat("Running fastAnc for temperature reconstruction...\n")
anc_temps <- fastAnc(pruned_tree, tip_temps, vars = TRUE, CI = TRUE)
hist(anc_temps$ace)
hist(anc_temps$var)
cat("ASR complete. Internal nodes:", length(anc_temps$ace), "\n")
cat("Temperature range (tips):", round(range(tip_temps), 2), "\n")
cat("Temperature range (internal):", round(range(anc_temps$ace), 2), "\n")

# PART 3: CLASSIFY EDGES BY PHENOTYPE SHIFT

cat("\n=== CLASSIFYING EDGES BY PHENOTYPE SHIFT ===\n")

# Combine tip and ancestral temperatures
n_tips <- Ntip(pruned_tree)
all_temps <- c(tip_temps, anc_temps$ace)
names(all_temps) <- c(pruned_tree$tip.label, (n_tips + 1):(n_tips + Nnode(pruned_tree)))

# Build edge temperature table
edges <- pruned_tree$edge
edge_temps <- data.table(
  edge_id = seq_len(nrow(edges)),
  parent_node = edges[, 1],
  child_node = edges[, 2],
  parent_temp = all_temps[as.character(edges[, 1])],
  child_temp = all_temps[as.character(edges[, 2])],
  edge_length = pruned_tree$edge.length
)

edge_temps[, delta_temp := child_temp - parent_temp]
edge_temps[, is_tip_edge := child_node <= n_tips]
par(mfrow=c(1,1))
hist(edge_temps$delta_temp)
# Classify edges by temperature change direction
# Using a threshold to avoid noise from small changes
temp_threshold <- sd(edge_temps$delta_temp, na.rm = TRUE) * 0.5

edge_temps[, pheno_class := fcase(
  delta_temp > temp_threshold, "hot",
  delta_temp < -temp_threshold, "cold",
  default = "stable"
)]

cat("\nEdge classification by temperature shift:\n")
print(table(edge_temps$pheno_class))

# PART 4: LINK SUBSTITUTIONS TO PHENOTYPE-CLASSIFIED EDGES

cat("\n=== LINKING SUBSTITUTIONS TO PHENOTYPE SHIFTS ===\n")

# Load GWAS results for site info
gwas_results <- rbindlist(lapply(
  list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE),
  function(f) {
    models <- readRDS(f)
    rbindlist(lapply(models, function(m) {
      data.table(
        Gene = m$Gene,
        Position = m$Aligned_Position,
        P_aa_with_pcs = m$P_aa_with_pcs,
        P_aa_only = m$P_aa_only,
        effects = list(m$effects)
      )
    }))
  }
))

# Merge with partitionMap
gwas_results[, merge_key := paste(Gene, Position, sep = "_")]
partitionMap[, merge_key := paste(Gene, GenePos, sep = "_")]

site_info <- merge(
  partitionMap[, .(GlobalPos, Gene, GenePos, merge_key)],
  gwas_results[, .(merge_key, P_aa_only, P_aa_with_pcs, effects)],
  by = "merge_key", all.x = TRUE
)

# Define GWAS significance classes
thresh_control <- quantile(site_info$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(site_info$P_aa_only, 0.20, na.rm = TRUE)

site_info[, sig_class := fcase(
  P_aa_with_pcs < thresh_control & P_aa_only < thresh_nocontrol, "sig_both",
  P_aa_with_pcs < thresh_control, "sig_control",
  P_aa_only < thresh_nocontrol, "sig_nocontrol",
  default = "not_sig"
)]

# Build node label lookup for pruned tree
node_labels_pruned <- c(pruned_tree$tip.label, as.character((n_tips + 1):(n_tips + Nnode(pruned_tree))))

# Match supermatrix rows to pruned tree tips
common_tips <- intersect(rownames(supermat_mat), pruned_tree$tip.label)
supermat_pruned <- supermat_mat[common_tips, ]

# Map internal node labels between trees
# anc_raw uses "Node1", "Node2" etc format
internal_mapping <- data.table(
  anc_node = anc_raw$node[!anc_raw$node %in% anc_tree$tip.label],
  seq_idx = which(!anc_raw$node %in% anc_tree$tip.label)
)

cat("Processing substitutions on phenotype-classified edges...\n")

pheno_trans_results <- list()
# ---- linking subs to asr pheno shifts ----
cat("\n=== LINKING SUBSTITUTIONS TO PHENOTYPE SHIFTS (FAST) ===\n")

# Step 1: Build state matrices for ALL positions at once
# Rows = edges, Cols = positions

n_edges <- nrow(edge_temps)
n_sites <- nrow(sites_to_process)

# Precompute node label lookups
tip_labels <- pruned_tree$tip.label
n_tips <- length(tip_labels)

# Parent states: all internal nodes
# Child states: tips or internal nodes

# For tips: use supermat_pruned directly
# For internal: use anc_raw

# Build edge-to-label mapping once
edge_temps[, parent_label := paste0("Node", parent_node - n_tips)]
edge_temps[, child_label := fifelse(
  child_node <= n_tips,
  tip_labels[child_node],
  paste0("Node", child_node - n_tips)
)]

# Check which internal node labels exist in anc_raw
valid_internal <- intersect(edge_temps$parent_label, anc_raw$node)
valid_internal_child <- intersect(edge_temps[child_node > n_tips]$child_label, anc_raw$node)

cat("Valid parent nodes:", length(valid_internal), "/", n_edges, "\n")

# Step 2: For each site, vectorized extraction
cat("Processing", n_sites, "sites...\n")

# Pre-extract all ancestral sequences into a matrix (much faster than repeated substr)
anc_seq_vec <- anc_raw$seq
names(anc_seq_vec) <- anc_raw$node

# Get all global positions
all_gpos <- sites_to_process$GlobalPos

pheno_trans_results <- vector("list", n_sites)

for (i in seq_len(n_sites)) {
  site <- sites_to_process[i]
  gpos <- site$GlobalPos
  
  
  # Tip states
  tip_state_vec <- supermat_pruned[, gpos]
  names(tip_state_vec) <- rownames(supermat_pruned)
  
  # Internal node states
  node_state_vec <- substr(anc_seq_vec, gpos, gpos)
  
  # Vectorized: get parent/child states for ALL edges at once
  parent_states <- node_state_vec[edge_temps$parent_label]
  child_states <- ifelse(
    edge_temps$child_node <= n_tips,
    tip_state_vec[edge_temps$child_label],
    node_state_vec[edge_temps$child_label]
  )
  
  # Find substitutions (vectorized)
  is_sub <- !is.na(parent_states) & !is.na(child_states) &
    parent_states != child_states &
    parent_states != "-" & child_states != "-"
  
  if (sum(is_sub) == 0) next
  
  # Extract only substitution edges
  sub_idx <- which(is_sub)
  
  # Get residue effect for derived alleles
  eff_dt <- site$effects[[1]]
  
  # Build result for this site
  result <- data.table(
    Gene = site$Gene,
    Position = site$GenePos,
    GlobalPos = gpos,
    edge_id = sub_idx,
    from_aa = parent_states[sub_idx],
    to_aa = child_states[sub_idx],
    delta_temp = edge_temps$delta_temp[sub_idx],
    pheno_class = edge_temps$pheno_class[sub_idx],
    P_aa_only = site$P_aa_only,
    P_aa_with_pcs = site$P_aa_with_pcs,
    sig_class = site$sig_class
  )
  
  # Add residue-level effects (vectorized lookup)
  if (!is.null(eff_dt) && nrow(eff_dt) > 0) {
    eff_dt[, aa := sub(".*([A-Z])$", "\\1", Residue)]
    result <- merge(result, eff_dt[, .(aa, P_residue = P_value, residue_effect = Effect)],
                    by.x = "to_aa", by.y = "aa", all.x = TRUE)
  } else {
    result[, `:=`(P_residue = NA_real_, residue_effect = NA_real_)]
  }
  
  pheno_trans_results[[i]] <- result
  
  if (i %% 500 == 0) cat("  Processed:", i, "/", n_sites, "\n")
}

pheno_trans_dt <- rbindlist(pheno_trans_results, fill = TRUE)
cat("\nTotal transitions:", nrow(pheno_trans_dt), "\n")

# Add variant-level significance classification
pheno_trans_dt[, variant_sig_class := sig_class]
pheno_trans_dt[P_residue > 0.05 | is.na(P_residue), variant_sig_class := "not_sig"]

cat("\nTransitions by phenotype class:\n")
print(table(pheno_trans_dt$pheno_class))

cat("\nTransitions by GWAS sig class:\n")
print(table(pheno_trans_dt$variant_sig_class))

# PART 5: COMPARE CLASSIFICATIONS

cat("\n=== COMPARING CLASSIFICATIONS ===\n")

# Q1: Do GWAS-significant transitions occur preferentially on warming/cooling edges?
cat("\n--- Q1: GWAS sig class vs phenotype class ---\n")

class_table <- table(pheno_trans_dt$sig_class, pheno_trans_dt$pheno_class)
print(class_table)

cat("\nChi-squared test:\n")
print(chisq.test(class_table))

# Proportions within each GWAS class
prop_by_gwas <- pheno_trans_dt[, .(
  prop_warming = mean(pheno_class == "hot"),
  prop_cooling = mean(pheno_class == "cold"),
  prop_stable = mean(pheno_class == "stable"),
  n = .N
), by = variant_sig_class]
print(prop_by_gwas)

# Q2: Do hot/cold effect variants match warming/cooling edges?
cat("\n--- Q2: Effect direction vs phenotype direction ---\n")

hist(pheno_trans_dt$residue_effect)
sd(na.omit(pheno_trans_dt$residue_effect)) * 0.5

pheno_trans_dt[, effect_dir := fcase(
  residue_effect > 5, "hot",
  residue_effect < -5, "cold",
  default = NA_character_
)]
pheno_trans_dt_full <- pheno_trans_dt
pheno_trans_dt$

table(pheno_trans_dt$effect_dir)
# Cross-tabulation
effect_pheno_table <- table(
  pheno_trans_dt[!is.na(effect_dir)]$effect_dir,
  pheno_trans_dt[!is.na(effect_dir)]$pheno_class
)
print(effect_pheno_table)

cat("\nChi-squared test (effect direction vs phenotype class):\n")
print(chisq.test(effect_pheno_table))

# Agreement rate: hot effects on warming edges, cold on cooling
agreement <- pheno_trans_dt[!is.na(effect_dir) & pheno_class != "stable", .(
  agree = sum((effect_dir == "hot" & pheno_class == "hot") |
                (effect_dir == "cold" & pheno_class == "cold")),
  disagree = sum((effect_dir == "hot" & pheno_class == "cold") |
                   (effect_dir == "cold" & pheno_class == "hot")),
  n = .N
)]
cat("\nAgreement (hot-hot or cold-cold):", agreement$agree, "\n")
cat("Disagreement (hot-cold or cold-cold):", agreement$disagree, "\n")
cat("Agreement ratio:", round(agreement$agree / (agreement$agree + agreement$disagree), 3), "\n")

agreement <- pheno_trans_dt[!is.na(effect_dir) & pheno_class != "stable" & variant_sig_class!="not_sig", .(
  agree = sum((effect_dir == "hot" & pheno_class == "hot") |
                (effect_dir == "cold" & pheno_class == "cold")),
  disagree = sum((effect_dir == "hot" & pheno_class == "cold") |
                   (effect_dir == "cold" & pheno_class == "hot")),
  n = .N
)]
cat("\nAgreement (hot-hot or cold-cold):", agreement$agree, "\n")
cat("Disagreement (hot-cold or cold-cold):", agreement$disagree, "\n")
cat("Agreement ratio:", round(agreement$agree / (agreement$agree + agreement$disagree), 3), "\n")

# PART 6: HYDROPHOBICITY ANALYSIS BY PHENOTYPE CLASS


cat("\n=== HYDROPHOBICITY BY PHENOTYPE CLASS ===\n")

# Add hydrophobicity
data(AAdata)
hydro <- AAdata$Hydrophobicity$Tanford

pheno_trans_dt[, to_hydro := hydro[to_aa]]
pheno_trans_dt[, from_hydro := hydro[from_aa]]
pheno_trans_dt[, delta_hydro := to_hydro - from_hydro]

# Mean delta hydrophobicity by phenotype class
hydro_by_pheno <- pheno_trans_dt[, .(
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  se_delta_hydro = sd(delta_hydro, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = pheno_class]
print(hydro_by_pheno)

cat("\nKruskal-Wallis: delta_hydro ~ pheno_class\n")
print(kruskal.test(delta_hydro ~ pheno_class, data = pheno_trans_dt))

# Pairwise tests
cat("\nWilcoxon: warming vs cooling\n")
print(wilcox.test(delta_hydro ~ pheno_class, 
                  data = pheno_trans_dt[pheno_class %in% c("hot", "cold")]))

# Stratify by GWAS significance
cat("\n--- Hydrophobicity by pheno class, stratified by GWAS class ---\n")
hydro_stratified <- pheno_trans_dt[, .(
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = .(variant_sig_class, pheno_class)]
print(dcast(hydro_stratified, variant_sig_class ~ pheno_class, value.var = "mean_delta_hydro"))

cols_gwas <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
               sig_control = "darkorange", sig_both = "firebrick")

# Filter out stable, order factors
plot_dt <- pheno_trans_dt[pheno_class != "stable" & !is.na(delta_hydro)]
plot_dt[, sig_class := factor(sig_class, levels = c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))]
plot_dt[, pheno_class := factor(pheno_class, levels = c("cold", "hot"))]

ggplot(plot_dt, aes(x = pheno_class, y = delta_hydro, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = cols_gwas) +
  facet_wrap(~variant_sig_class, nrow = 1) +
  labs(x = "Phenotype shift direction", y = "Delta hydrophobicity",
       title = "Hydrophobicity change by phenotype shift, stratified by GWAS class") +
  theme_bw() +
  theme(legend.position = "none")

# Compute p-values per facet
pvals <- plot_dt[, .(
  p = wilcox.test(delta_hydro ~ pheno_class)$p.value
), by = variant_sig_class]
pvals[, label := fifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p))]

ggplot(plot_dt, aes(x = pheno_class, y = delta_hydro, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = cols_gwas) +
  geom_signif(comparisons = list(c("cold", "hot")), 
              test = "wilcox.test", 
              map_signif_level = FALSE,
              textsize = 3) +
  facet_wrap(~variant_sig_class, nrow = 1) +
  labs(x = "Phenotype shift direction", y = "Delta hydrophobicity",
       title = "Hydrophobicity change by phenotype shift, stratified by variant significance") +
  theme_bw() +
  theme(legend.position = "none")


# Count transitions per site
site_counts <- pheno_trans_dt[, .N, by = .(GlobalPos, from_aa, to_aa)]

# Keep only transitions that occur >= 5 times at that site
pheno_trans_dt_filtered <- merge(
  pheno_trans_dt,
  site_counts[N >= 10, .(GlobalPos, from_aa, to_aa)],
  by = c("GlobalPos", "from_aa", "to_aa")
)

cat("Before filter:", nrow(pheno_trans_dt), "\n")
cat("After filter:", nrow(pheno_trans_dt_filtered), "\n")

# Now plot with filtered data
plot_dt <- pheno_trans_dt_filtered[pheno_class != "stable" & !is.na(delta_hydro)]
plot_dt[, variant_sig_class := factor(variant_sig_class, levels = c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))]
plot_dt[, pheno_class := factor(pheno_class, levels = c("cold", "hot"))]

pvals <- plot_dt[, .(
  p = wilcox.test(delta_hydro ~ pheno_class)$p.value,
  y = max(delta_hydro, na.rm = TRUE) + 0.5
), by = variant_sig_class]
pvals[, label := fifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p))]

ggplot(plot_dt, aes(x = pheno_class, y = delta_hydro, fill = variant_sig_class)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = cols_gwas) +
  geom_text(data = pvals, aes(x = 1.5, y = y, label = label), 
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~variant_sig_class, nrow = 1) +
  labs(x = "Phenotype shift direction", y = "Delta hydrophobicity",
       title = "Hydrophobicity change by phenotype shift") +
  theme_bw() +
  theme(legend.position = "none")

table(plot_dt$variant_sig_class)

# PART 7: TRANSITION MATRICES BY PHENOTYPE CLASS

cat("\n=== TRANSITION MATRICES BY PHENOTYPE CLASS ===\n")

# Build transition matrices
all_aa <- sort(unique(c(pheno_trans_dt_filtered$from_aa, pheno_trans_dt$to_aa)))
all_aa <- all_aa[all_aa != "-" & !is.na(all_aa)]

build_trans_matrix <- function(counts_dt, aa_levels) {
  mat <- matrix(0, nrow = length(aa_levels), ncol = length(aa_levels),
                dimnames = list(from = aa_levels, to = aa_levels))
  for (i in seq_len(nrow(counts_dt))) {
    from <- counts_dt$from_aa[i]
    to <- counts_dt$to_aa[i]
    if (from %in% aa_levels && to %in% aa_levels) {
      mat[from, to] <- counts_dt$N[i]
    }
  }
  mat
}

trans_counts <- pheno_trans_dt_filtered[, .N, by = .(pheno_class, from_aa, to_aa)]

pheno_trans_mats <- lapply(c("hot", "cold", "stable"), function(cl) {
  build_trans_matrix(trans_counts[pheno_class == cl], all_aa)
})
names(pheno_trans_mats) <- c("hot", "cold", "stable")

# Log-odds: warming vs cooling
calc_log_odds <- function(mat1, mat2, pseudocount = 0.5) {
  p1 <- (mat1 + pseudocount) / sum(mat1 + pseudocount)
  p2 <- (mat2 + pseudocount) / sum(mat2 + pseudocount)
  log2(p1 / p2)
}

lor_warm_cool <- calc_log_odds(pheno_trans_mats[["hot"]], pheno_trans_mats[["cold"]])

lor_dt <- as.data.table(as.table(lor_warm_cool))
setnames(lor_dt, c("from_aa", "to_aa", "log2_OR"))

# Add counts
warm_counts <- as.data.table(as.table(pheno_trans_mats[["hot"]]))
setnames(warm_counts, c("from_aa", "to_aa", "N_warming"))
cool_counts <- as.data.table(as.table(pheno_trans_mats[["cold"]]))
setnames(cool_counts, c("from_aa", "to_aa", "N_cooling"))

lor_dt <- merge(lor_dt, warm_counts, by = c("from_aa", "to_aa"))
lor_dt <- merge(lor_dt, cool_counts, by = c("from_aa", "to_aa"))
lor_dt[, N_total := N_warming + N_cooling]

# Add delta hydrophobicity
lor_dt[, delta_hydro := hydro[to_aa] - hydro[from_aa]]

cat("\nTop transitions enriched in WARMING edges:\n")
print(lor_dt[N_total >= 10][order(-log2_OR)][1:15, 
                                             .(from_aa, to_aa, log2_OR, delta_hydro, N_warming, N_cooling)])

cat("\nTop transitions enriched in COOLING edges:\n")
print(lor_dt[N_total >= 10][order(log2_OR)][1:15, 
                                            .(from_aa, to_aa, log2_OR, delta_hydro, N_warming, N_cooling)])

# Correlation: log-odds vs delta_hydro
cat("\nCorrelation: warming/cooling log-odds vs delta_hydro\n")
print(cor.test(lor_dt[N_total >= 5]$log2_OR, lor_dt[N_total >= 5]$delta_hydro, method = "spearman"))

# PART 8: COMPARE WITH GWAS-BASED HOT/COLD CLASSIFICATION

cat("\n=== COMPARING PHENOTYPE-BASED VS GWAS-BASED CLASSIFICATIONS ===\n")

# Within sig_both sites, compare hydrophobicity patterns
sig_both_trans <- pheno_trans_dt_filtered[variant_sig_class == "sig_both" & !is.na(effect_dir)]
summary(sig_both_trans)
cat("\nWithin sig_both sites:\n")

# GWAS-based classification
gwas_hydro <- sig_both_trans[, .(
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = effect_dir]
cat("By GWAS effect direction:\n")
print(gwas_hydro)

# Phenotype-based classification
pheno_hydro <- sig_both_trans[pheno_class != "stable", .(
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = pheno_class]
cat("\nBy phenotype shift direction:\n")
print(pheno_hydro)

# Combined classification
combined_hydro <- sig_both_trans[pheno_class != "stable", .(
  mean_delta_hydro = mean(delta_hydro, na.rm = TRUE),
  n = .N
), by = .(effect_dir, pheno_class)]
cat("\nCombined (effect_dir x pheno_class):\n")
print(dcast(combined_hydro, effect_dir ~ pheno_class, value.var = "mean_delta_hydro"))

# Test interaction
if (nrow(sig_both_trans[pheno_class != "stable"]) > 50) {
  cat("\nTwo-way ANOVA: delta_hydro ~ effect_dir * pheno_class\n")
  print(summary(aov(delta_hydro ~ effect_dir * pheno_class, 
                    data = sig_both_trans[pheno_class != "stable"])))
}

# ============================================================================
# PART 9: SAVE RESULTS
# ============================================================================

cat("\n=== SAVING RESULTS ===\n")

saveRDS(list(
  pheno_trans_dt = pheno_trans_dt,
  edge_temps = edge_temps,
  anc_temps = anc_temps,
  lor_dt = lor_dt,
  pheno_trans_mats = pheno_trans_mats
), "results/phenotype_asr_results.rds")

cat("Results saved to: results/phenotype_asr_results.rds\n")

# ============================================================================
# PART 10: VISUALIZATION
# ============================================================================

cat("\n=== GENERATING PLOTS ===\n")

library(ggplot2)

cols_pheno <- c(warming = "tomato", cooling = "steelblue", stable = "grey70")
cols_gwas <- c(not_sig = "grey70", sig_nocontrol = "steelblue", 
               sig_control = "darkorange", sig_both = "firebrick")

# Plot 1: Delta hydrophobicity by phenotype class
p1 <- ggplot(pheno_trans_dt_filtered[variant_sig_class=="sig_both"], 
             aes(x = pheno_class, y = delta_hydro, fill = pheno_class)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Delta hydrophobicity",
       title = "Hydrophobicity change by phenotype shift direction") +
  theme_bw() +
  theme(legend.position = "none")

# Plot 2: Combined effect_dir x pheno_class (sig_both only)
p2 <- ggplot(sig_both_trans[pheno_class != "stable"], 
             aes(x = pheno_class, y = delta_hydro, fill = effect_dir)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c(cold = "steelblue", hot = "tomato")) +
  labs(x = "Phenotype shift", y = "Delta hydrophobicity",
       title = "sig_both sites: GWAS effect vs phenotype shift",
       fill = "GWAS effect") +
  theme_bw()

# Plot 3: Log-odds vs delta_hydro
p3 <- ggplot(lor_dt[N_total >= 5], aes(x = delta_hydro, y = log2_OR, size = N_total)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  labs(x = "Delta hydrophobicity (to - from)",
       y = "Log2 OR (warming / cooling)",
       title = "Phenotype-based transition enrichment vs hydrophobicity",
       size = "N transitions") +
  theme_bw()

# Save plots
pdf("results/phenotype_asr_plots.pdf", width = 10, height = 8)
print(p1)
print(p2)
print(p3)
dev.off()

cat("Plots saved to: results/phenotype_asr_plots.pdf\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# ---- pheno ASR analysis ----

