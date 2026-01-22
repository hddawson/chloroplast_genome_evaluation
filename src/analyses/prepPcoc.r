#!/usr/bin/env Rscript
# ==============================================================================
# Prepare Data for PCOC Convergence Analysis
#
# PCOC requires:
# 1. AA alignment (fasta)
# 2. Tree with branch lengths (newick)
# 3. Scenario string defining convergent clades
#
# This script identifies warm-adapted clades and generates the scenario string
# ==============================================================================

library(ape)
library(data.table)
library(phangorn)
library(arrow)

# ==== CONFIGURATION ====
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
output_dir <- "pcoc_input"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==== LOAD DATA ====

cat("Loading tree...\n")
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

# --- EDA: Phenotype distribution ---
par(mfrow = c(1, 2))
hist(tip_pheno, breaks = 50, main = "Tip phenotypes", 
     xlab = "Mean temp wettest quarter", col = "steelblue")

# Define "warm" threshold - use upper quartile for stringent definition
warm_thresh <- quantile(tip_pheno, 0.25, na.rm = TRUE)
abline(v = warm_thresh, col = "red", lwd = 2)
legend("topright", paste("Warm threshold:", round(warm_thresh, 1)), bty = "n")

tip_is_warm <- !(tip_pheno > warm_thresh & !is.na(tip_pheno))
barplot(table(tip_is_warm), main = "Warm vs Cold tips",
        names.arg = c("Cold", "Warm"), col = c("dodgerblue", "tomato"))
par(mfrow = c(1, 1))

cat("\nWarm threshold (75th percentile):", round(warm_thresh, 2), "\n")
cat("Warm tips:", sum(tip_is_warm, na.rm = TRUE), "\n")
cat("Cold tips:", sum(!tip_is_warm, na.rm = TRUE), "\n")
#STUPID INCOMING I SWAPPED WARM COLD SO COLD IS THE CONVERGENT SCENARIO
# BUT IN THIS IT IS WARM - SWAP WARM AND COLD WHENEVER U READ EM 

# ==== IDENTIFY CONVERGENT WARM CLADES ====

#' Find monophyletic clades where all tips are warm
#' Returns the MRCA node of each such clade

find_warm_clades <- function(tree, tip_is_warm, min_clade_size = 3) {
  
  warm_tips <- which(tip_is_warm)
  warm_tip_names <- tree$tip.label[warm_tips]
  
  cat("Finding warm clades from", length(warm_tips), "warm tips...\n")
  
  # Strategy: find connected components of warm tips
  # A "warm clade" is a monophyletic group where ALL tips are warm
  
  # Get all internal nodes
  internal_nodes <- (n_tips + 1):(n_tips + n_internal)
  
  warm_clades <- list()
  
  for (node in internal_nodes) {
    # Get all descendant tips
    desc_tips <- Descendants(tree, node, type = "tips")[[1]]
    
    if (length(desc_tips) < min_clade_size) next
    
    # Check if ALL descendants are warm
    all_warm <- all(tip_is_warm[desc_tips])
    
    if (all_warm) {
      # Check this isn't a subset of a larger warm clade we already found
      # (we want maximal clades)
      is_subset <- FALSE
      for (existing in warm_clades) {
        if (all(desc_tips %in% existing$tips)) {
          is_subset <- TRUE
          break
        }
      }
      
      if (!is_subset) {
        # Remove any existing clades that are subsets of this one
        warm_clades <- Filter(function(x) !all(x$tips %in% desc_tips), warm_clades)
        
        warm_clades[[length(warm_clades) + 1]] <- list(
          mrca = node,
          tips = desc_tips,
          n_tips = length(desc_tips)
        )
      }
    }
  }
  
  cat("Found", length(warm_clades), "maximal warm clades\n")
  warm_clades
}

warm_clades <- find_warm_clades(tree, tip_is_warm, min_clade_size = 3)

# --- EDA: Clade size distribution ---
clade_sizes <- sapply(warm_clades, function(x) x$n_tips)
hist(clade_sizes, breaks = 20, main = "Warm clade sizes",
     xlab = "Number of tips", col = "tomato")

cat("\nClade size summary:\n")
print(summary(clade_sizes))

# ==== FILTER TO INDEPENDENT CLADES ====

# PCOC wants truly independent convergent events
# Remove nested clades and keep only well-separated ones

filter_independent_clades <- function(warm_clades, tree, min_separation = 5) {
  
  if (length(warm_clades) <= 1) return(warm_clades)
  
  # Sort by size (largest first)
  warm_clades <- warm_clades[order(-sapply(warm_clades, function(x) x$n_tips))]
  
  # Check pairwise: are any clades nested or too close?
  keep <- rep(TRUE, length(warm_clades))
  
  for (i in 1:(length(warm_clades) - 1)) {
    if (!keep[i]) next
    
    for (j in (i + 1):length(warm_clades)) {
      if (!keep[j]) next
      
      # Check if j is descendant of i or vice versa
      node_i <- warm_clades[[i]]$mrca
      node_j <- warm_clades[[j]]$mrca
      
      # Get ancestors
      anc_i <- Ancestors(tree, node_i, type = "all")
      anc_j <- Ancestors(tree, node_j, type = "all")
      
      if (node_j %in% anc_i || node_i %in% anc_j) {
        # Nested - remove the smaller one
        keep[j] <- FALSE
      }
    }
  }
  
  warm_clades[keep]
}

independent_clades <- filter_independent_clades(warm_clades, tree)
cat("\nIndependent warm clades:", length(independent_clades), "\n")

# --- EDA: Show independent clades ---
cat("\nIndependent clade details:\n")
for (i in seq_along(independent_clades)) {
  clade <- independent_clades[[i]]
  cat("Clade", i, ": MRCA node", clade$mrca, ",", clade$n_tips, "tips\n")
  cat("  Example tips:", paste(head(tree$tip.label[clade$tips], 3), collapse = ", "), "\n")
}

# ==== BUILD PCOC SCENARIO STRING ====

#' Convert ape node numbers to PCOC format
#' PCOC uses 0-indexed tip numbers, then internal nodes
#' ape uses 1-indexed: tips 1:n_tips, internals (n_tips+1):(n_tips+n_internal)

build_pcoc_scenario <- function(independent_clades, tree) {
  
  n_tips <- Ntip(tree)
  
  # PCOC numbering: tips are 0:(n_tips-1), internals follow
  # Need to map ape node IDs to PCOC node IDs
  
  ape_to_pcoc <- function(ape_node) {
    if (ape_node <= n_tips) {
      # It's a tip: PCOC uses 0-indexing
      return(ape_node - 1)
    } else {
      # It's internal: PCOC numbers after tips
      # But we need to figure out PCOC's internal numbering...
      # This is tricky - PCOC numbers by tree traversal
      # For now, assume same ordering (may need adjustment)
      return(ape_node - 1)
    }
  }
  
  scenario_parts <- list()
  
  for (i in seq_along(independent_clades)) {
    clade <- independent_clades[[i]]
    
    # Get all nodes in this clade (tips + internal)
    mrca <- clade$mrca
    all_desc <- Descendants(tree, mrca, type = "all")[[1]]
    all_nodes <- c(mrca, all_desc)
    
    # Convert to PCOC numbering
    pcoc_nodes <- sapply(all_nodes, ape_to_pcoc)
    
    # MRCA must be first
    mrca_pcoc <- ape_to_pcoc(mrca)
    other_pcoc <- setdiff(pcoc_nodes, mrca_pcoc)
    
    scenario_parts[[i]] <- paste(c(mrca_pcoc, other_pcoc), collapse = ",")
  }
  
  paste(scenario_parts, collapse = "/")
}

scenario_string <- build_pcoc_scenario(independent_clades, tree)
cat("\n=== PCOC Scenario String ===\n")
cat("(First 500 chars):\n")
cat(substr(scenario_string, 1, 500), "...\n")

# ==== ALTERNATIVE: USE PCOC_NUM_TREE TO GET CORRECT NUMBERING ====

# The safest approach is to:
# 1. Export tree
# 2. Run pcoc_num_tree.py to get numbered tree
# 3. Manually verify/adjust scenario string

# For now, let's also create a simpler scenario with just the largest clades

# Take top N independent clades by size
top_n <- length(independent_clades)
top_clades <- independent_clades[order(-sapply(independent_clades, function(x) x$n_tips))][1:top_n]

cat("\n=== Top", top_n, "clades for PCOC ===\n")
for (i in seq_along(top_clades)) {
  clade <- top_clades[[i]]
  cat("Clade", i, ":", clade$n_tips, "tips, MRCA =", clade$mrca, "\n")
}

scenario_top <- build_pcoc_scenario(top_clades, tree)

# ==== PREPARE OUTPUT FILES ====

# 1. Tree file (newick)
tree_out <- file.path(output_dir, "tree.nw")
outtree <- tree
outtree$node.label <- NULL
outtree$edge.label <- NULL
write.tree(outtree, tree_out)
cat("\nTree written to:", tree_out, "\n")

# 2. Alignment file - PCOC needs fasta, you already have this
# Just copy or symlink
library(seqinr)

aln <- read.fasta(aln_file, seqtype="AA", as.string=FALSE)
aln <- lapply(aln, function(x) { x[x %in% c("Z","J","B")] <- "-"; x })
write.fasta(aln, names(aln), aln_out)
m <- do.call(rbind, aln)
gap_per_site <- colSums(m == "-")

hist(gap_per_site, breaks=50, xlab="Fraction gaps per site")

# 3. Scenario string
scenario_out <- file.path(output_dir, "scenario.txt")
writeLines(scenario_string, scenario_out)
cat("Full scenario written to:", scenario_out, "\n")

scenario_top_out <- file.path(output_dir, "scenario_top10.txt")
writeLines(scenario_top, scenario_top_out)
cat("Top-10 scenario written to:", scenario_top_out, "\n")

# 4. Warm tips list (for verification)
warm_tips_out <- file.path(output_dir, "warm_tips.txt")
writeLines(tree$tip.label[tip_is_warm], warm_tips_out)
cat("Warm tips list written to:", warm_tips_out, "\n")

# 5. Clade info table
clade_info <- rbindlist(lapply(seq_along(independent_clades), function(i) {
  clade <- independent_clades[[i]]
  data.table(
    clade_id = i,
    mrca_ape = clade$mrca,
    n_tips = clade$n_tips,
    example_tips = paste(head(tree$tip.label[clade$tips], 5), collapse = ";")
  )
}))
clade_info_out <- file.path(output_dir, "clade_info.tsv")
fwrite(clade_info, clade_info_out, sep = "\t")
cat("Clade info written to:", clade_info_out, "\n")

# ==== GENERATE PCOC COMMANDS ====

cat("\n")
cat("=======================================================\n")
cat("PCOC COMMANDS\n")
cat("=======================================================\n")
cat("\n")
cat("# 1. First, number the tree to verify node IDs:\n")
cat("singularity exec pcoc_latest.sif pcoc_num_tree.py \\\n")
"
singularity exec pcoc_latest.sif pcoc_num_tree.py \
  -t pcoc_input/tree.nw \
  -o pcoc_input/numbered_tree.pdf 
"
cat("  -t", tree_out, "\\\n")
cat("  -o", file.path(output_dir, "numbered_tree.pdf"), "\n")
cat("\n")
cat("# 2. Visualize scenario (check colored branches):\n")
cat("singularity exec pcoc_latest.sif pcoc_num_tree.py \\\n")
cat("  -t", tree_out, "\\\n")
cat("  -o", file.path(output_dir, "scenario_tree.pdf"), "\\\n")
cat("  -m \"$(cat", scenario_top_out, ")\"\n")
cat("\n")
cat("# 3. Run PCOC detection:\n")
cat("singularity exec pcoc_latest.sif pcoc_det.py \\\n")
cat("  -t", tree_out, "\\\n")
cat("  -aa", aln_out, "\\\n")
cat("  -o", file.path(output_dir, "pcoc_results"), "\\\n")
cat("  -m \"$(cat", scenario_top_out, ")\" \\\n")
cat("  -f 0.8 \\\n")
cat("  --plot --plot_complete_ali\n")
cat("\n")
cat("=======================================================\n")

# ==== VISUALIZE WARM CLADES ON TREE ====

cat("\nGenerating tree visualization...\n")

# Color tips by warm/cold
tip_colors <- ifelse(tip_is_warm, "tomato", "steelblue")
tip_colors[is.na(tip_is_warm)] <- "gray"

# For large trees, plot a subset or use clade highlighting
if (n_tips > 200) {
  cat("Tree too large for full visualization, plotting clade structure...\n")
  
  # Plot showing which clades are warm
  par(mfrow = c(1, 1))
  
  # Create a simplified view: barplot of clade sizes
  barplot(sort(clade_sizes, decreasing = TRUE), 
          main = paste("Warm clade sizes (n =", length(clade_sizes), "clades)"),
          xlab = "Clade rank", ylab = "Number of tips",
          col = "tomato", border = NA)
  
} else {
  # Small enough to plot
  plot(tree, type = "fan", show.tip.label = FALSE, 
       edge.width = 0.5, main = "Warm (red) vs Cold (blue) tips")
  tiplabels(pch = 16, col = tip_colors, cex = 0.5)
}

# ==== SANITY CHECKS ====

cat("\n=== Sanity Checks ===\n")

# Check tree and alignment have same tips
aln <- read.FASTA(aln_file, type = "AA")
aln_names <- names(aln)
tree_names <- tree$tip.label

in_both <- intersect(aln_names, tree_names)
only_tree <- setdiff(tree_names, aln_names)
only_aln <- setdiff(aln_names, tree_names)

cat("Tips in both:", length(in_both), "\n")
cat("Only in tree:", length(only_tree), "\n")
cat("Only in alignment:", length(only_aln), "\n")

stopifnot(length(only_tree) == 0)
stopifnot(length(only_aln) == 0)

cat("\n=== Setup Complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Run the PCOC commands above to perform detection.\n")
cat("\nNOTE: You may need to verify the scenario string against\n")
cat("the numbered tree output from pcoc_num_tree.py\n")