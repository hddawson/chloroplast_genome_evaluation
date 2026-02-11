# ============================================================================
# CHLOROPLAST GWAS 3D STRUCTURE ANALYSIS - EXAMPLE WORKFLOW
# ============================================================================
# 
# This script demonstrates the analysis workflow for mapping temperature
# adaptation GWAS signals onto protein 3D structures.
#
# Prerequisites:
#   - sites_df: GWAS results (Gene, Position, P_aa_only, P_aa_with_pcs, N)
#   - struct_var: Secondary structure predictions (Gene, Position, consensus, ...)
#   - Aligned FASTA files in data/tmp/alignedGenes/
#   - data object with ID and Organism columns for reference lookup
#
# ============================================================================

# Source the analysis functions
source("src/analyses/structureAnalysis.r")

library(arrow)
data <- read_parquet("data/processed_data.parquet")

model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_list <- lapply(model_files, function(f) {
  cat(f)
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_only = m$P_aa_only,
      P_aa_with_pcs = m$P_aa_with_pcs,
      N = m$N
    )
  }))
})
sites_df <- rbindlist(gwas_list)


# ---- SETUP: Define thresholds and paths ----
# Paths
ALIGNMENT_DIR <- "data/tmp/alignedGenes"

# Reference species pattern - adjust based on your data
# For structures from spinach, we need to find a close relative in our alignment
# Here we use Arabidopsis as it's commonly present
REFERENCE_PATTERN <- data$ID[grep("thaliana", data$Organism)]

# Define significance thresholds (genome-wide quantiles)
# These should be computed ONCE from full sites_df
message("Computing genome-wide significance thresholds...")
THRESH_AA_ONLY <- quantile(sites_df$P_aa_only, 0.25, na.rm = TRUE)
THRESH_AA_PCS <- quantile(sites_df$P_aa_with_pcs, 0.05, na.rm = TRUE)

message("Thresholds:")
message("  P_aa_only < ", signif(THRESH_AA_ONLY, 3), " (25th percentile)")
message("  P_aa_with_pcs < ", signif(THRESH_AA_PCS, 3), " (5th percentile)")

# Quick summary of what we're working with
message("\nGWAS data summary:")
message("  Genes: ", length(unique(sites_df$Gene)))
message("  Total sites: ", nrow(sites_df))
message("  Sites P_aa_only < threshold: ", sum(sites_df$P_aa_only < THRESH_AA_ONLY, na.rm = TRUE))
message("  Sites P_aa_with_pcs < threshold: ", sum(sites_df$P_aa_with_pcs < THRESH_AA_PCS, na.rm = TRUE))


# ============================================================================
# EXAMPLE 1: ATP SYNTHASE
# ============================================================================

message("\n\n")
message(paste(rep("*", 70), collapse = ""))
message("*  EXAMPLE 1: ATP SYNTHASE")
message(paste(rep("*", 70), collapse = ""))

# --- Step 1: Load the structure ---
# Download PDB manually then load
atp_structure <- load_pdb_structure("6FKF")
# Examine what we got
message("\nStructure summary:")
print(atp_structure$chain_summary)

message("\nLigands present:")
print(atp_structure$ligands[, .(ligand_type, n = .N), by = ligand_type])

# Visualize the complex
p_atp_overview <- ggplot(atp_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() +
  labs(title = "ATP synthase (6FKF) - XZ projection",
       subtitle = "Colored by chain") +
  theme_classic() +
  theme(legend.position = "none")

p_atp_overview_xy <- ggplot(atp_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() +
  labs(title = "ATP synthase (6FKF) - XZ projection",
       subtitle = "Colored by chain") +
  theme_classic() +
  theme(legend.position = "none")

p_atp_overview + p_atp_overview_xy


# --- Step 2: Analyze one gene manually (atpB) to understand the workflow ---

message("\n--- Detailed walkthrough: atpB ---\n")

# Define the complex
atp_def <- DEFAULT_COMPLEXES$atp_synthase

# atpB has three chains (β subunits): B, D, F
atpB_chains <- atp_def$chain_map$atpB
message("atpB chains: ", paste(atpB_chains, collapse = ", "))

# Map all three chains to the alignment
atpB_maps <- map_multiple_chains(
  atp_structure,
  file.path(ALIGNMENT_DIR, "atpB_AA_aligned.fasta"),
  REFERENCE_PATTERN,
  "atpB",
  atpB_chains
)

# Check: how many residues mapped per chain?
message("\nResidues mapped per chain:")
print(atpB_maps[, .N, by = chain])

# Check: are the chains in different conformations?
# Look at ligand distances per chain
nuc_coords <- atp_structure$ligands[ligand_type %in% c("ATP", "ADP"), 
                                     .(x = lig_x, y = lig_y, z = lig_z)]

atpB_maps[, dist_to_nucleotide := sapply(1:.N, function(i) {
  min(sqrt((x[i] - nuc_coords$x)^2 + (y[i] - nuc_coords$y)^2 + (z[i] - nuc_coords$z)^2))
})]

# Compare conformations
message("\nConformational states (distance to nearest nucleotide):")
print(atpB_maps[, .(
  mean_dist_nuc = mean(dist_to_nucleotide),
  min_dist_nuc = min(dist_to_nucleotide),
  n_near_nuc = sum(dist_to_nucleotide < 15)
), by = chain])

# This shows the three β subunits are in different states:
# F - βTP (tight, ATP-bound): closest to nucleotide
# D - βDP (loose, ADP-bound): intermediate
# B - βE (empty): farthest from nucleotide


# --- Step 3: Compute aggregated features ---

# Get other chains for interface calculation
other_chains <- setdiff(unique(atp_structure$ca_df$chain), atpB_chains)

atpB_features <- compute_residue_features(
  atpB_maps, 
  atp_structure, 
  atp_def,
  other_chains
)

# Add axis distance (c-ring defines the rotation axis)
atpB_features <- add_axis_distance(
  atpB_features, 
  atp_structure, 
  atp_def$chain_map$atpH
)

message("\nFeature summary:")
print(summary(atpB_features[, .(dist_ligand_min, dist_ligand_mean, 
                                 dist_interface_min, conformational_range,
                                 z_mean, dist_to_axis)]))

# Key insight: conformational_range shows how much each residue moves
# between the three β states
message("\nTop 10 most flexible residues (highest conformational range):")
print(atpB_features[order(-conformational_range)][1:10, 
                    .(aln_pos, conformational_range, z_mean, at_interface, near_ligand)])


# --- Step 4: Merge with GWAS and run tests ---

atpB_merged <- merge_gwas_with_structure(
  sites_df, 
  atpB_features,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

plot(atpB_merged$x_mean, atpB_merged$y_mean)

message("\nSite classification:")
print(table(atpB_merged$site_class))

# Run all enrichment tests
atpB_tests <- run_all_enrichment_tests(atpB_merged)
print_test_summary(atpB_tests)


# --- Step 5: Visualize ---

# Structure projections
p_atpB_xy <- plot_gwas_on_structure(atpB_merged, atp_structure, atpB_maps, "xy")
p_atpB_xz <- plot_gwas_on_structure(atpB_merged, atp_structure, atpB_maps, "xz")

# Feature boxplots
p_box_lig <- plot_feature_boxplot(atpB_merged, "dist_ligand_min", 
                                   "Distance to nucleotide (Å)", "atpB: Ligand proximity")
p_box_int <- plot_feature_boxplot(atpB_merged, "dist_interface_min",
                                   "Distance to interface (Å)", "atpB: Interface proximity")
p_box_axis <- plot_feature_boxplot(atpB_merged, "dist_to_axis",
                                    "Distance to axis (Å)", "atpB: Radial position")
p_box_z <- plot_feature_boxplot(atpB_merged, "z_mean",
                                 "Z coordinate (Å)", "atpB: Membrane position")

# Display
print(p_atpB_xy + p_atpB_xz)
print((p_box_lig + p_box_int) / (p_box_axis + p_box_z))


# --- Step 6: Analyze ALL ATP synthase genes at once ---

message("\n\n--- Full ATP synthase analysis ---\n")

# Use the high-level function
atp_results <- analyze_complex(
  "atp_synthase",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Combine test results across all genes
atp_combined_tests <- combine_test_results(atp_results)

message("\n=== COMBINED ATP SYNTHASE RESULTS ===\n")

message("Continuous features (significant results only):")
sig_continuous <- atp_combined_tests$continuous[wilcox_p < 0.05]
if (nrow(sig_continuous) > 0) {
  print(sig_continuous[order(wilcox_p), .(gene, feature, site_class, 
                                           median_sig, median_bg, direction, wilcox_p)])
} else {
  message("  No significant results")
}

message("\nBinary features (significant results only):")
sig_binary <- atp_combined_tests$binary[fisher_p < 0.05]
if (nrow(sig_binary) > 0) {
  print(sig_binary[order(fisher_p), .(gene, feature, site_class,
                                       pct_sig, pct_bg, odds_ratio, fisher_p)])
} else {
  message("  No significant results")
}


# --- Hypothesis testing on combined data ---

message("\n=== HYPOTHESIS TESTS ON COMBINED ATP SYNTHASE DATA ===\n")

combined_atp <- atp_results$combined_merged

# H1: Are sig sites more stromal (higher Z)?
message("H1: Stromal enrichment")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- combined_atp[site_class == sc & !is.na(z_mean), z_mean]
  bg_z <- combined_atp[site_class == "not_sig" & !is.na(z_mean), z_mean]
  if (length(test_z) >= 3) {
    wt <- wilcox.test(test_z, bg_z)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f", 
                    sc, median(test_z), median(bg_z), wt$p.value))
  }
}

# H2: Are sig sites more peripheral (far from axis)?
message("\nH2: Peripheral enrichment")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_atp[site_class == sc & !is.na(dist_to_axis), dist_to_axis]
  bg_d <- combined_atp[site_class == "not_sig" & !is.na(dist_to_axis), dist_to_axis]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f",
                    sc, median(test_d), median(bg_d), wt$p.value))
  }
}

# H3: Do sig sites avoid ligand binding sites?
message("\nH3: Ligand avoidance")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_atp[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- combined_atp[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f",
                    sc, median(test_d), median(bg_d), wt$p.value))
  }
}

# Combined visualization
p_atp_all_xz <- plot_complex_gwas(atp_results, "xz")
p_atp_all_xy <- plot_complex_gwas(atp_results, "xy")

p_atp_all_xz + p_atp_all_xy


# Plot atpF only
p_atpF_xy <- plot_gwas_on_structure(
  atp_results$gene_results$atpF$merged, 
  atp_results$structure, 
  atp_results$gene_results$atpF$position_maps, 
  "xy", 
  "atpF - XY (top view)"
)

p_atpF_xz <- plot_gwas_on_structure(
  atp_results$gene_results$atpF$merged, 
  atp_results$structure, 
  atp_results$gene_results$atpF$position_maps, 
  "xz", 
  "atpF - XZ (side view)"
)

print(p_atpF_xy + p_atpF_xz)

# Check the data
print(atp_results$gene_results$atpF$merged[site_class != "not_sig", 
                                           .(Position, P_aa_only, P_aa_with_pcs, site_class, z_mean, compartment)])


# Check what alignment type we're using and how atpF maps
atpF_maps <- atp_results$gene_results$atpF$position_maps

# Look at the mapping
print(head(atpF_maps, 20))
print(tail(atpF_maps, 20))

# Compare Z coordinates - are they inverted?
print(atpF_maps[, .(aln_pos, z, chain)][order(aln_pos)])

# What's the Z range?
print(atpF_maps[, .(min_z = min(z), max_z = max(z))])

# Re-source and re-run

# Re-analyze ATP synthase
atp_results <- analyze_complex(
  "atp_synthase",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Check atpF specifically
p_atpF_xz <- plot_gwas_on_structure(
  atp_results$gene_results$atpF$merged, 
  atp_results$structure, 
  atp_results$gene_results$atpF$position_maps, 
  "xz"
)
print(p_atpF_xz)

# Re-run full ATP synthase analysis with fixed mapping

# Combine test results across all genes
atp_combined_tests <- combine_test_results(atp_results)

# Print summary of significant results
message("\n=== SIGNIFICANT CONTINUOUS FEATURES ===")
print(atp_combined_tests$continuous[wilcox_p < 0.05][order(wilcox_p), 
                                                     .(gene, feature, site_class, n_sig, median_sig, median_bg, direction, wilcox_p)])

message("\n=== SIGNIFICANT BINARY FEATURES ===")
print(atp_combined_tests$binary[fisher_p < 0.05][order(fisher_p),
                                                 .(gene, feature, site_class, n_sig_in, pct_sig, pct_bg, odds_ratio, fisher_p)])

# ---- HYPOTHESIS TESTS ON COMBINED DATA ----

combined_atp <- atp_results$combined_merged

message("\n", paste(rep("=", 60), collapse = ""))
message("COMBINED ATP SYNTHASE HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

# H1: Stromal enrichment (higher Z)?
message("\nH1: Stromal enrichment (higher Z)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- combined_atp[site_class == sc & !is.na(z_mean), z_mean]
  bg_z <- combined_atp[site_class == "not_sig" & !is.na(z_mean), z_mean]
  if (length(test_z) >= 3) {
    wt <- wilcox.test(test_z, bg_z)
    dir <- ifelse(median(test_z) > median(bg_z), "↑", "↓")
    message(sprintf("  %s (n=%d): median=%.1f %s vs %.1f, p=%.4f", 
                    sc, length(test_z), median(test_z), dir, median(bg_z), wt$p.value))
  }
}

# H2: Peripheral (far from axis)?
message("\nH2: Peripheral (far from rotation axis)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_atp[site_class == sc & !is.na(dist_to_axis), dist_to_axis]
  bg_d <- combined_atp[site_class == "not_sig" & !is.na(dist_to_axis), dist_to_axis]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ peripheral", "↓ central")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: Avoid ligand binding sites?
message("\nH3: Avoid ligand binding sites")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_atp[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- combined_atp[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ avoiding", "↓ closer")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H4: Interface enrichment?
message("\nH4: Interface enrichment (<8Å from other subunit)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_atp[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_atp[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ in bulk")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# ---- SUMMARY PLOTS ----

# 1. Structure overview with all GWAS hits
p_struct_xy <- plot_complex_gwas(atp_results, "xy", "ATP synthase - top view")
p_struct_xz <- plot_complex_gwas(atp_results, "xz", "ATP synthase - side view")

# 2. Hypothesis test boxplots
combined_atp[, site_class := factor(site_class, 
                                    levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

p_h1 <- ggplot(combined_atp[!is.na(z_mean)], aes(x = site_class, y = z_mean, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = c(105, 135), linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H1: Membrane position", y = "Z coordinate (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h2 <- ggplot(combined_atp[!is.na(dist_to_axis)], aes(x = site_class, y = dist_to_axis, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H2: Radial position", y = "Distance to axis (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h3 <- ggplot(combined_atp[!is.na(dist_ligand_min)], aes(x = site_class, y = dist_ligand_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H3: Ligand proximity", y = "Distance to ligand (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h4 <- ggplot(combined_atp[!is.na(dist_interface_min)], aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = 8, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H4: Interface proximity", y = "Distance to interface (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Display
print(p_struct_xy + p_struct_xz)
print((p_h1 + p_h2) / (p_h3 + p_h4))

# 3. Per-gene structure plots (3 projections each)
for (gene_name in names(atp_results$gene_results)) {
  gr <- atp_results$gene_results[[gene_name]]
  if (is.null(gr)) next
  
  p_xy <- plot_gwas_on_structure(gr$merged, atp_results$structure, gr$position_maps, 
                                 "xy", paste(gene_name, "- XY"))
  p_xz <- plot_gwas_on_structure(gr$merged, atp_results$structure, gr$position_maps, 
                                 "xz", paste(gene_name, "- XZ"))
  p_yz <- plot_gwas_on_structure(gr$merged, atp_results$structure, gr$position_maps, 
                                 "yz", paste(gene_name, "- YZ"))
  
  print(p_xy + p_xz + p_yz + plot_layout(ncol = 3))
}

# 4. Summary table
message("\n=== SITE COUNTS BY GENE ===")
site_summary <- combined_atp[, .(
  total = .N,
  with_structure = sum(!is.na(x_mean)),
  sig_no_ctrl = sum(site_class == "sig_no_control"),
  sig_with_ctrl = sum(site_class == "sig_with_control"),
  sig_both = sum(site_class == "sig_both")
), by = Gene]
print(site_summary)

# ---- SURFACE ACCESSIBILITY APPROXIMATION ----

#' Compute neighbor-based surface accessibility proxy
#' 
#' Buried residues have many neighbors, surface residues have few.
#' This is a rough approximation - true RSA requires DSSP or FreeSASA.
#' 
#' @param ca_df CA atom coordinates (data.table with x, y, z)
#' @param radius Shell radius in Angstroms (default 10)
#' @return Vector of neighbor counts (lower = more surface exposed)
compute_neighbor_count <- function(x, y, z, radius = 10) {
  n <- length(x)
  sapply(1:n, function(i) {
    dists <- sqrt((x - x[i])^2 + (y - y[i])^2 + (z - z[i])^2)
    sum(dists > 0 & dists <= radius)
  })
}

atp_structure$ca_df$neighbor_count <- with(atp_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

# Add neighbor count (surface proxy) to combined data
# Need to compute per-chain, then aggregate

# Compute for all CA atoms in structure
# Check distribution
message("Neighbor count distribution (all CA atoms):")
print(summary(atp_structure$ca_df$neighbor_count))
hist(atp_structure$ca_df$neighbor_count)

# Now add to position_maps for each gene and aggregate
for (gene_name in names(atp_results$gene_results)) {
  gr <- atp_results$gene_results[[gene_name]]
  if (is.null(gr)) next
  
  # Get neighbor counts for residues in this gene's chains
  gr$position_maps <- merge(
    gr$position_maps,
    atp_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
}

# Aggregate across chains and add to combined_merged
# Aggregate directly without modifying atp_results
surface_features <- rbindlist(lapply(names(atp_results$gene_results), function(gene_name) {
  gr <- atp_results$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  # Merge here instead
  pm <- merge(
    gr$position_maps,
    atp_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
  
  pm[, .(
    neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
    neighbor_count_min = min(neighbor_count, na.rm = TRUE),
    neighbor_count_max = max(neighbor_count, na.rm = TRUE)
  ), by = .(gene, aln_pos)]
}))

# Merge with combined data
combined_atp <- merge(combined_atp, surface_features,
                      by.x = c("Gene", "Position"),
                      by.y = c("gene", "aln_pos"),
                      all.x = TRUE)

# Define "surface exposed" as below median neighbor count
median_neighbors <- median(combined_atp$neighbor_count_min, na.rm = TRUE)
combined_atp[, surface_exposed := neighbor_count_min < median_neighbors]

message("\nNeighbor count by site class (lower = more surface exposed):")
print(combined_atp[!is.na(neighbor_count_min), .(
  n = .N,
  median_neighbors = as.numeric(median(neighbor_count_min)),
  pct_surface = 100 * mean(surface_exposed, na.rm = TRUE)
), by = site_class])

# ---- H5: Surface exposure test ----

message("\n", paste(rep("=", 60), collapse = ""))
message("H5: Are significant sites more SURFACE EXPOSED?")
message("(Lower neighbor count = more exposed)")
message(paste(rep("=", 60), collapse = ""))

for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_atp[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_atp[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# Fisher's test for surface enrichment
message("\nSurface exposure enrichment (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_surf <- sum(combined_atp$site_class == sc & combined_atp$surface_exposed == TRUE, na.rm = TRUE)
  sig_buried <- sum(combined_atp$site_class == sc & combined_atp$surface_exposed == FALSE, na.rm = TRUE)
  bg_surf <- sum(combined_atp$site_class == "not_sig" & combined_atp$surface_exposed == TRUE, na.rm = TRUE)
  bg_buried <- sum(combined_atp$site_class == "not_sig" & combined_atp$surface_exposed == FALSE, na.rm = TRUE)
  
  if (sig_surf + sig_buried == 0) next
  mat <- matrix(c(sig_surf, sig_buried, bg_surf, bg_buried), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("  %s: %d/%d surface (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_surf, sig_surf + sig_buried,
                  100 * sig_surf / (sig_surf + sig_buried),
                  100 * bg_surf / (bg_surf + bg_buried),
                  ft$estimate, ft$p.value))
}

# ---- Add to summary plot ----

p_h5 <- ggplot(combined_atp[!is.na(neighbor_count_min)], 
               aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = median_neighbors, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H5: Surface exposure", 
       subtitle = "Lower = more exposed",
       y = "Neighbor count (10Å)", x = "") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Updated combined plot
print((p_h1 + p_h2 + p_h3) / (p_h4 + p_h5))

# Correlation with p-value
cor_surface <- cor.test(combined_atp$neighbor_count_min, combined_atp$neglog_p_only, 
                        method = "spearman", use = "complete")
message(sprintf("\nCorrelation (neighbor count vs -log10 P_only): rho = %.3f, p = %.2e",
                cor_surface$estimate, cor_surface$p.value))
message("(Negative rho = more exposed residues have higher significance)")

plot(combined_atp$neighbor_count_min, combined_atp$neglog_p_pcs)

# Visualize surface vs buried on structure
# Color by neighbor count (low = surface/exposed, high = buried)

# Visualize surface vs buried on structure
p_surface_xz <- ggplot() +
  geom_point(data = atp_structure$ca_df, 
             aes(x = x, y = z, color = neighbor_count),
             size = 0.8, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = "Neighbors\n(low=exposed)") +
  coord_fixed() +
  labs(title = "ATP synthase - neighbor count (surface proxy)",
       subtitle = "Yellow = exposed, Purple = buried",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

print(p_surface_xz)

# Get sig_both positions with per-chain coordinates
sig_both_coords <- rbindlist(lapply(names(atp_results$gene_results), function(gene_name) {
  gr <- atp_results$gene_results[[gene_name]]
  if (is.null(gr)) return(NULL)
  
  # Merge position_maps with site_class
  merge(
    gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
    gr$merged[site_class == "sig_both", .(Gene, Position)],
    by.x = c("gene", "aln_pos"),
    by.y = c("Gene", "Position")
  )
}))

# Overlay sig_both sites 
p_surface_sig_xz <- ggplot() +
  geom_point(data = atp_structure$ca_df, 
             aes(x = x, y = z, color = neighbor_count),
             size = 0.5, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = "Neighbors") +
  geom_point(data = sig_both_coords,
             aes(x = x, y = z),
             color = "red", size = 3, shape = 1, stroke = 1.5) +
  coord_fixed() +
  labs(title = "sig_both sites (red circles) on surface map",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_surface_sig_yz <- ggplot() +
  geom_point(data = atp_structure$ca_df, 
             aes(x = y, y = z, color = neighbor_count),
             size = 0.5, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = "Neighbors") +
  geom_point(data = sig_both_coords,
             aes(x = x, y = z),
             color = "red", size = 3, shape = 1, stroke = 1.5) +
  coord_fixed() +
  labs(title = "sig_both sites (red circles) on surface map",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()


print(p_surface_sig_xz + p_surface_sig_yz)

# Histogram comparison
p_hist <- ggplot(combined_atp[!is.na(neighbor_count_min)], 
                 aes(x = neighbor_count_min, fill = site_class == "not_sig")) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("TRUE" = "gray50", "FALSE" = "red"),
                    labels = c("TRUE" = "not_sig", "FALSE" = "significant"),
                    name = "") +
  geom_vline(xintercept = median_neighbors, linetype = "dashed") +
  labs(title = "Neighbor count distribution",
       subtitle = "Dashed line = median (surface threshold)",
       x = "Neighbor count (10Å)", y = "Count") +
  theme_classic()

print(p_hist)


# Color structure by significance level!

# First, add neighbor_count to position_maps and merge with site_class
all_coords_with_sig <- rbindlist(lapply(names(atp_results$gene_results), function(gene_name) {
  gr <- atp_results$gene_results[[gene_name]]
  if (is.null(gr)) return(NULL)
  
  # Merge position_maps with site_class from merged data
  merge(
    gr$position_maps[, .(gene, aln_pos, chain, pdb_resno, x, y, z)],
    gr$merged[, .(Gene, Position, site_class, P_aa_only, P_aa_with_pcs)],
    by.x = c("gene", "aln_pos"),
    by.y = c("Gene", "Position"),
    all.x = TRUE
  )
}))

all_coords_with_sig[is.na(site_class), site_class := "no_data"]
all_coords_with_sig[, neglog_p := -log10(P_aa_only)]

# Plot 1: Color by site_class (categorical)
p_sig_class_xz <- ggplot() +
  # Background atoms not in our genes (gray)
  geom_point(data = atp_structure$ca_df[!chain %in% unique(all_coords_with_sig$chain)],
             aes(x = x, y = z), color = "gray90", size = 0.5, alpha = 0.3) +
  # Our genes colored by significance
  geom_point(data = all_coords_with_sig[site_class == "no_data"],
             aes(x = x, y = z), color = "gray70", size = 1) +
  geom_point(data = all_coords_with_sig[site_class == "not_sig"],
             aes(x = x, y = z), color = "gray50", size = 1.2) +
  geom_point(data = all_coords_with_sig[site_class == "sig_no_control"],
             aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_coords_with_sig[site_class == "sig_with_control"],
             aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_coords_with_sig[site_class == "sig_both"],
             aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() +
  labs(title = "ATP synthase colored by GWAS significance",
       subtitle = "Gray=not_sig, Gold=sig_no_ctrl, Blue=sig_with_ctrl, Red=sig_both",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_sig_class_xy <- p_sig_class_xz + aes(x = x, y = y) + 
  labs(title = "ATP synthase - top view", x = "X (Å)", y = "Y (Å)")

# Plot 2: Continuous color by -log10(P)
p_pval_xz <- ggplot() +
  geom_point(data = atp_structure$ca_df[!chain %in% unique(all_coords_with_sig$chain)],
             aes(x = x, y = z), color = "gray90", size = 0.5, alpha = 0.3) +
  geom_point(data = all_coords_with_sig[!is.na(neglog_p)],
             aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() +
  labs(title = "ATP synthase colored by -log10(P_aa_only)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_pval_xy <- ggplot() +
  geom_point(data = atp_structure$ca_df[!chain %in% unique(all_coords_with_sig$chain)],
             aes(x = x, y = y), color = "gray90", size = 0.5, alpha = 0.3) +
  geom_point(data = all_coords_with_sig[!is.na(neglog_p)],
             aes(x = x, y = y, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() +
  labs(title = "ATP synthase - top view",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

print(p_sig_class_xz + p_sig_class_xy)
print(p_pval_xz + p_pval_xy)

# Plot 3: Combine - surface exposure + significance
# Merge neighbor_count onto all_coords
all_coords_with_sig <- merge(
  all_coords_with_sig,
  atp_structure$ca_df[, .(chain, resno, neighbor_count)],
  by.x = c("chain", "pdb_resno"),
  by.y = c("chain", "resno"),
  all.x = TRUE
)

# Scatterplot: neighbor count vs p-value
p_scatter <- ggplot(all_coords_with_sig[!is.na(neglog_p) & !is.na(neighbor_count)],
                    aes(x = neighbor_count, y = neglog_p, color = site_class)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("not_sig" = "gray50", "no_data" = "gray80",
                                "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", 
                                "sig_both" = "darkred")) +
  geom_smooth(data = all_coords_with_sig[!is.na(neglog_p) & !is.na(neighbor_count)],
              aes(x = neighbor_count, y = neglog_p), 
              inherit.aes = FALSE, method = "loess", color = "black", se = TRUE) +
  labs(title = "Surface exposure vs significance",
       subtitle = "Lower neighbor count = more surface exposed",
       x = "Neighbor count (10Å)", y = "-log10(P_aa_only)") +
  theme_classic()

print(p_scatter)

# ---- rbcL ----
source("src/analyses/structureAnalysis.r")
message("\n\n")
message(paste(rep("*", 70), collapse = ""))
message("*  RUBISCO ANALYSIS")
message(paste(rep("*", 70), collapse = ""))

# Load structure
rbc_structure <- load_pdb_structure("1RCX")

message("\nStructure summary:")
print(rbc_structure$chain_summary)

# Rubisco has L8S8 architecture:
# Large subunits: L, M, J, K, B, C, H, I (plastid-encoded)
# Small subunits: S, T, A, D, E, F, G, N (nuclear-encoded)

message("\nLigands:")
print(rbc_structure$ligands[, .N, by = ligand_type])

rbc_def <- DEFAULT_COMPLEXES$rubisco

rbcL_result <- analyze_gene(
  
  "rbcL",
  rbc_structure,
  rbc_def,
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Print test results
print_test_summary(rbcL_result$tests)

# ---- Visualize with per-chain coordinates ----
p_rbcL_xy <- plot_gwas_on_structure(rbcL_result$merged, rbc_structure, 
                                    rbcL_result$position_maps, "xy", "RbcL - XY (top view)")
p_rbcL_xz <- plot_gwas_on_structure(rbcL_result$merged, rbc_structure, 
                                    rbcL_result$position_maps, "xz", "RbcL - XZ (side view)")
p_rbcL_yz <- plot_gwas_on_structure(rbcL_result$merged, rbc_structure, 
                                    rbcL_result$position_maps, "yz", "RbcL - YZ (side view)")

print(p_rbcL_xy + p_rbcL_xz)

# ---- Feature boxplots ----
p_rbcL_lig <- plot_feature_boxplot(rbcL_result$merged, "dist_ligand_min",
                                   "Distance to RuBP (Å)", "RbcL: Ligand proximity")
p_rbcL_int <- plot_feature_boxplot(rbcL_result$merged, "dist_interface_min",
                                   "Distance to interface (Å)", "RbcL: Interface proximity")

print(p_rbcL_lig + p_rbcL_int)

# ---- Surface accessibility ----
rbc_structure$ca_df$neighbor_count <- with(rbc_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

message("\nNeighbor count distribution (Rubisco):")
print(summary(rbc_structure$ca_df$neighbor_count))

# Aggregate surface features for rbcL
rbcL_surface <- merge(
  rbcL_result$position_maps,
  rbc_structure$ca_df[, .(chain, resno, neighbor_count)],
  by.x = c("chain", "pdb_resno"),
  by.y = c("chain", "resno"),
  all.x = TRUE
)[, .(
  neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
  neighbor_count_min = min(neighbor_count, na.rm = TRUE),
  neighbor_count_max = max(neighbor_count, na.rm = TRUE)
), by = .(gene, aln_pos)]

# Merge with rbcL results
rbcL_merged <- merge(rbcL_result$merged, rbcL_surface,
                     by.x = c("Gene", "Position"),
                     by.y = c("gene", "aln_pos"),
                     all.x = TRUE)

median_neighbors_rbc <- median(rbcL_merged$neighbor_count_min, na.rm = TRUE)
rbcL_merged[, surface_exposed := neighbor_count_min < median_neighbors_rbc]

# ---- L-L dimer interface (Rubisco-specific) ----
message("\n--- Rubisco L-L interface analysis ---\n")

rbcL_chains <- rbc_def$chain_map$rbcL
L_coords <- rbc_structure$ca_df[chain %in% rbcL_chains]

# Distance to other L subunits
rbcL_result$position_maps[, dist_to_L_interface := sapply(1:.N, function(i) {
  my_chain <- chain[i]
  other_L <- L_coords[chain != my_chain]
  min(sqrt((x[i] - other_L$x)^2 + (y[i] - other_L$y)^2 + (z[i] - other_L$z)^2))
})]

# Distance to S subunits (L-S interface)
S_chains <- setdiff(unique(rbc_structure$ca_df$chain), rbcL_chains)
S_coords <- rbc_structure$ca_df[chain %in% S_chains]

rbcL_result$position_maps[, dist_to_S_interface := sapply(1:.N, function(i) {
  min(sqrt((x[i] - S_coords$x)^2 + (y[i] - S_coords$y)^2 + (z[i] - S_coords$z)^2))
})]

# Aggregate
rbcL_interface_features <- rbcL_result$position_maps[, .(
  dist_LL_min = min(dist_to_L_interface),
  dist_LL_mean = mean(dist_to_L_interface),
  dist_LS_min = min(dist_to_S_interface),
  dist_LS_mean = mean(dist_to_S_interface)
), by = .(gene, aln_pos)]

rbcL_merged <- merge(rbcL_merged, rbcL_interface_features,
                     by.x = c("Gene", "Position"),
                     by.y = c("gene", "aln_pos"),
                     all.x = TRUE)

rbcL_merged[, at_LL_interface := dist_LL_min < 8]
rbcL_merged[, at_LS_interface := dist_LS_min < 8]

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("RUBISCO HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

# H1: Ligand proximity
message("\nH1: Distance to ligand (RuBP)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- rbcL_merged[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- rbcL_merged[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ avoiding", "↓ closer")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H2: L-L interface
message("\nH2: Distance to L-L interface")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- rbcL_merged[site_class == sc & !is.na(dist_LL_min), dist_LL_min]
  bg_d <- rbcL_merged[site_class == "not_sig" & !is.na(dist_LL_min), dist_LL_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: L-S interface
message("\nH3: Distance to L-S interface")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- rbcL_merged[site_class == sc & !is.na(dist_LS_min), dist_LS_min]
  bg_d <- rbcL_merged[site_class == "not_sig" & !is.na(dist_LS_min), dist_LS_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H4: Surface exposure
message("\nH4: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- rbcL_merged[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- rbcL_merged[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# Fisher's tests
message("\nSurface enrichment (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_surf <- sum(rbcL_merged$site_class == sc & rbcL_merged$surface_exposed == TRUE, na.rm = TRUE)
  sig_buried <- sum(rbcL_merged$site_class == sc & rbcL_merged$surface_exposed == FALSE, na.rm = TRUE)
  bg_surf <- sum(rbcL_merged$site_class == "not_sig" & rbcL_merged$surface_exposed == TRUE, na.rm = TRUE)
  bg_buried <- sum(rbcL_merged$site_class == "not_sig" & rbcL_merged$surface_exposed == FALSE, na.rm = TRUE)
  
  if (sig_surf + sig_buried == 0) next
  mat <- matrix(c(sig_surf, sig_buried, bg_surf, bg_buried), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("  %s: %d/%d surface (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_surf, sig_surf + sig_buried,
                  100 * sig_surf / (sig_surf + sig_buried),
                  100 * bg_surf / (bg_surf + bg_buried),
                  ft$estimate, ft$p.value))
}

# ---- VISUALIZATION ----

# Color structure by -log10(P)
rbcL_coords_with_sig <- merge(
  rbcL_result$position_maps[, .(gene, aln_pos, chain, pdb_resno, x, y, z)],
  rbcL_result$merged[, .(Gene, Position, site_class, P_aa_only, P_aa_with_pcs)],
  by.x = c("gene", "aln_pos"),
  by.y = c("Gene", "Position"),
  all.x = TRUE
)
rbcL_coords_with_sig[, neglog_p := -log10(P_aa_only)]

p_rbc_pval_xy <- ggplot() +
  geom_point(data = rbc_structure$ca_df[!chain %in% rbcL_chains],
             aes(x = x, y = y), color = "gray90", size = 0.5, alpha = 0.3) +
  geom_point(data = rbcL_coords_with_sig[!is.na(neglog_p)],
             aes(x = x, y = y, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() +
  labs(title = "Rubisco L8S8 - colored by -log10(P_aa_only)",
       subtitle = "Gray = small subunits (nuclear)",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_rbc_pval_xz <- ggplot() +
  geom_point(data = rbc_structure$ca_df[!chain %in% rbcL_chains],
             aes(x = x, y = z), color = "gray90", size = 0.5, alpha = 0.3) +
  geom_point(data = rbcL_coords_with_sig[!is.na(neglog_p)],
             aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() +
  labs(title = "Rubisco - XZ view",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

print(p_rbc_pval_xy + p_rbc_pval_xz)

# Surface map with sig sites
p_rbc_surface <- ggplot() +
  geom_point(data = rbc_structure$ca_df,
             aes(x = x, y = y, color = neighbor_count),
             size = 0.8, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = "Neighbors\n(low=exposed)") +
  coord_fixed() +
  labs(title = "Rubisco surface map",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

print(p_rbc_surface)

# Boxplots
rbcL_merged[, site_class := factor(site_class, 
                                   levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

p_box_lig <- plot_feature_boxplot(rbcL_merged, "dist_ligand_min", 
                                  "Distance to RuBP (Å)", "Ligand proximity")
p_box_LL <- plot_feature_boxplot(rbcL_merged, "dist_LL_min",
                                 "Distance to L-L interface (Å)", "L-L interface")
p_box_LS <- plot_feature_boxplot(rbcL_merged, "dist_LS_min",
                                 "Distance to L-S interface (Å)", "L-S interface")
p_box_surf <- ggplot(rbcL_merged[!is.na(neighbor_count_min)], 
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = median_neighbors_rbc, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Surface exposure", y = "Neighbor count", x = "") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print((p_box_lig + p_box_LL) / (p_box_LS + p_box_surf))

# ---- Summary ----
message("\n=== SITE COUNTS ===")
print(rbcL_merged[, .(
  total = .N,
  with_structure = sum(!is.na(dist_ligand_min)),
  sig_no_ctrl = sum(site_class == "sig_no_control"),
  sig_with_ctrl = sum(site_class == "sig_with_control"),
  sig_both = sum(site_class == "sig_both")
)])

# Save
saveRDS(rbcL_result, "results/rbcL_3d_analysis.rds")
saveRDS(rbcL_merged, "results/rbcL_merged_with_features.rds")

# ---- photosynthetic ----
# ANALYSIS OF PHOTOSYNTHETIC COMPLEXES: Cyt b6f, PSII, PSI
# 
# Run AFTER sourcing structureAnalysis.r and loading GWAS data
#

# STEP 0: VALIDATE CHAIN MAPPINGS BEFORE ANALYSIS
# The DEFAULT_COMPLEXES chain mappings may be wrong (like rubisco was).
# Let's check each structure first.

message("\n", paste(rep("=", 70), collapse = ""))
message("VALIDATING PDB STRUCTURES AND CHAIN MAPPINGS")
message(paste(rep("=", 70), collapse = ""))

# --- Cytochrome b6f (6RQF) ---
message("\n--- Cytochrome b6f (6RQF) ---")
b6f_structure <- load_pdb_structure("6RQF")
print(b6f_structure$chain_summary)
message("\nLigands:")
print(b6f_structure$ligands[, .N, by = ligand_type])

b6f_def <- DEFAULT_COMPLEXES$cytochrome_b6f

b6f_result <- analyze_complex(
  "cytochrome_b6f",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Combine test results across all genes
b6f_combined_tests <- combine_test_results(b6f_result)

message("\n=== COMBINED b6f SYNTHASE RESULTS ===\n")

message("Continuous features (significant results only):")
sig_continuous <- b6f_combined_tests$continuous[wilcox_p < 0.05]
if (nrow(sig_continuous) > 0) {
  print(sig_continuous[order(wilcox_p), .(gene, feature, site_class, 
                                          median_sig, median_bg, direction, wilcox_p)])
} else {
  message("  No significant results")
}

message("\nBinary features (significant results only):")
sig_binary <- b6f_combined_tests$binary[fisher_p < 0.05]
if (nrow(sig_binary) > 0) {
  print(sig_binary[order(fisher_p), .(gene, feature, site_class,
                                      pct_sig, pct_bg, odds_ratio, fisher_p)])
} else {
  message("  No significant results")
}


# --- Hypothesis testing on combined data ---

message("\n=== HYPOTHESIS TESTS ON COMBINED b6f DATA ===\n")

combined_b6f <- b6f_result$combined_merged

# H1: Are sig sites more stromal (higher Z)?
message("H1: Stromal enrichment")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- combined_b6f[site_class == sc & !is.na(z_mean), z_mean]
  bg_z <- combined_b6f[site_class == "not_sig" & !is.na(z_mean), z_mean]
  if (length(test_z) >= 3) {
    wt <- wilcox.test(test_z, bg_z)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f", 
                    sc, median(test_z), median(bg_z), wt$p.value))
  }
}



# H3: Do sig sites avoid ligand binding sites?
message("\nH3: Ligand avoidance")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_b6f[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- combined_b6f[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f",
                    sc, median(test_d), median(bg_d), wt$p.value))
  }
}

# Combined visualization
p_b6f_all_xz <- plot_complex_gwas(b6f_result, "xz")
p_b6f_all_xy <- plot_complex_gwas(b6f_result, "xy")

p_b6f_all_xz + p_b6f_all_xy

b6f_structure$ca_df$neighbor_count <- with(b6f_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

# Add neighbor count (surface proxy) to combined data
# Need to compute per-chain, then aggregate

# Compute for all CA atoms in structure
# Check distribution
message("Neighbor count distribution (all CA atoms):")
print(summary(b6f_structure$ca_df$neighbor_count))
hist(b6f_structure$ca_df$neighbor_count)

# Now add to position_maps for each gene and aggregate
for (gene_name in names(b6f_result$gene_results)) {
  gr <- b6f_result$gene_results[[gene_name]]
  if (is.null(gr)) next
  
  # Get neighbor counts for residues in this gene's chains
  gr$position_maps <- merge(
    gr$position_maps,
    b6f_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
}

# Aggregate across chains and add to combined_merged
# Aggregate directly without modifying atp_results
surface_features <- rbindlist(lapply(names(b6f_result$gene_results), function(gene_name) {
  gr <- b6f_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  # Merge here instead
  pm <- merge(
    gr$position_maps,
    b6f_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
  
  pm[, .(
    neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
    neighbor_count_min = min(neighbor_count, na.rm = TRUE),
    neighbor_count_max = max(neighbor_count, na.rm = TRUE)
  ), by = .(gene, aln_pos)]
}))

# Merge with combined data
combined_b6f <- merge(combined_b6f, surface_features,
                      by.x = c("Gene", "Position"),
                      by.y = c("gene", "aln_pos"),
                      all.x = TRUE)

# Define "surface exposed" as below median neighbor count
median_neighbors <- median(combined_b6f$neighbor_count_min, na.rm = TRUE)
combined_b6f[, surface_exposed := neighbor_count_min < median_neighbors]

message("\nNeighbor count by site class (lower = more surface exposed):")
print(combined_b6f[!is.na(neighbor_count_min), .(
  n = .N,
  median_neighbors = as.numeric(median(neighbor_count_min)),
  pct_surface = 100 * mean(surface_exposed, na.rm = TRUE)
), by = site_class])


boxplot(neighbor_count_min ~ site_class, combined_b6f)


# ---- PSII (5XNL) ----
message("\n--- Photosystem II (5XNL) ---")
psII_structure <- load_pdb_structure("5XNL")
print(psII_structure$chain_summary)
message("\nLigands:")
print(psII_structure$ligands[, .N, by = ligand_type])

psII_def <- DEFAULT_COMPLEXES$photosystem_ii

psII_result <- analyze_complex(
  "photosystem_ii",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Combine test results across all genes
psII_combined_tests <- combine_test_results(psII_result)

message("\n=== COMBINED psII SYNTHASE RESULTS ===\n")

message("Continuous features (significant results only):")
sig_continuous <- psII_combined_tests$continuous[wilcox_p < 0.05]
if (nrow(sig_continuous) > 0) {
  print(sig_continuous[order(wilcox_p), .(gene, feature, site_class, 
                                          median_sig, median_bg, direction, wilcox_p)])
} else {
  message("  No significant results")
}

message("\nBinary features (significant results only):")
sig_binary <- psII_combined_tests$binary[fisher_p < 0.05]
if (nrow(sig_binary) > 0) {
  print(sig_binary[order(fisher_p), .(gene, feature, site_class,
                                      pct_sig, pct_bg, odds_ratio, fisher_p)])
} else {
  message("  No significant results")
}


# --- Hypothesis testing on combined data ---

message("\n=== HYPOTHESIS TESTS ON COMBINED psII DATA ===\n")

combined_psII <- psII_result$combined_merged

# H1: Are sig sites more stromal (higher Z)?
message("H1: Stromal enrichment")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- combined_psII[site_class == sc & !is.na(z_mean), z_mean]
  bg_z <- combined_psII[site_class == "not_sig" & !is.na(z_mean), z_mean]
  if (length(test_z) >= 3) {
    wt <- wilcox.test(test_z, bg_z)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f", 
                    sc, median(test_z), median(bg_z), wt$p.value))
  }
}



# H3: Do sig sites avoid ligand binding sites?
message("\nH3: Ligand avoidance")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_psII[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- combined_psII[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f",
                    sc, median(test_d), median(bg_d), wt$p.value))
  }
}

# Combined visualization
p_psII_all_xz <- plot_complex_gwas(psII_result, "xz")
p_psII_all_xy <- plot_complex_gwas(psII_result, "xy")

p_psII_all_xz + p_psII_all_xy

psII_structure$ca_df$neighbor_count <- with(psII_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

# Add neighbor count (surface proxy) to combined data
# Need to compute per-chain, then aggregate

# Compute for all CA atoms in structure
# Check distribution
message("Neighbor count distribution (all CA atoms):")
print(summary(psII_structure$ca_df$neighbor_count))
hist(psII_structure$ca_df$neighbor_count)

# Now add to position_maps for each gene and aggregate
for (gene_name in names(psII_result$gene_results)) {
  gr <- psII_result$gene_results[[gene_name]]
  if (is.null(gr)) next
  
  # Get neighbor counts for residues in this gene's chains
  gr$position_maps <- merge(
    gr$position_maps,
    psII_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
}

# Aggregate across chains and add to combined_merged
# Aggregate directly without modifying atp_results
surface_features <- rbindlist(lapply(names(psII_result$gene_results), function(gene_name) {
  gr <- psII_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  # Merge here instead
  pm <- merge(
    gr$position_maps,
    psII_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
  
  pm[, .(
    neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
    neighbor_count_min = min(neighbor_count, na.rm = TRUE),
    neighbor_count_max = max(neighbor_count, na.rm = TRUE)
  ), by = .(gene, aln_pos)]
}))

# Merge with combined data
combined_psII <- merge(combined_psII, surface_features,
                      by.x = c("Gene", "Position"),
                      by.y = c("gene", "aln_pos"),
                      all.x = TRUE)

# Define "surface exposed" as below median neighbor count
median_neighbors <- median(combined_psII$neighbor_count_min, na.rm = TRUE)
combined_psII[, surface_exposed := neighbor_count_min < median_neighbors]

message("\nNeighbor count by site class (lower = more surface exposed):")
print(combined_psII[!is.na(neighbor_count_min), .(
  n = .N,
  median_neighbors = as.numeric(median(neighbor_count_min)),
  pct_surface = 100 * mean(surface_exposed, na.rm = TRUE)
), by = site_class])


boxplot(neighbor_count_min ~ site_class, combined_psII)

#psII pigments
# PSII PIGMENT PROXIMITY ANALYSIS

# Ligand key for PSII:
# CHL = chlorophyll (generic)
# CLA = chlorophyll a
# LUT = lutein (carotenoid)
# XAT, NEX = xanthophylls (carotenoids)
# BCR = beta-carotene
# PHO = pheophytin
# HEM = heme

# Get coordinates for different pigment classes
chl_coords <- psII_structure$ligands[ligand_type %in% c("CHL", "CLA"), 
                                     .(x = lig_x, y = lig_y, z = lig_z)]
carotenoid_coords <- psII_structure$ligands[ligand_type %in% c("LUT", "XAT", "NEX", "BCR"),
                                            .(x = lig_x, y = lig_y, z = lig_z)]
pheophytin_coords <- psII_structure$ligands[ligand_type == "PHO",
                                            .(x = lig_x, y = lig_y, z = lig_z)]

message("Pigment counts:")
message("  Chlorophylls (CHL+CLA): ", nrow(chl_coords))
message("  Carotenoids (LUT+XAT+NEX+BCR): ", nrow(carotenoid_coords))
message("  Pheophytins: ", nrow(pheophytin_coords))

# Compute distances for each position in combined_psII
# Using the per-chain coordinates from position_maps

pigment_distances <- rbindlist(lapply(names(psII_result$gene_results), function(gene_name) {
  gr <- psII_result$gene_results[[gene_name]]
  if (is.null(gr$position_maps)) return(NULL)
  
  pm <- gr$position_maps
  
  pm[, `:=`(
    dist_chlorophyll = sapply(1:.N, function(i) {
      min(sqrt((x[i] - chl_coords$x)^2 + (y[i] - chl_coords$y)^2 + (z[i] - chl_coords$z)^2))
    }),
    dist_carotenoid = sapply(1:.N, function(i) {
      min(sqrt((x[i] - carotenoid_coords$x)^2 + (y[i] - carotenoid_coords$y)^2 + (z[i] - carotenoid_coords$z)^2))
    }),
    dist_pheophytin = if (nrow(pheophytin_coords) > 0) {
      sapply(1:.N, function(i) {
        min(sqrt((x[i] - pheophytin_coords$x)^2 + (y[i] - pheophytin_coords$y)^2 + (z[i] - pheophytin_coords$z)^2))
      })
    } else NA_real_
  )]
  
  # Aggregate across chains
  pm[, .(
    dist_chl_min = min(dist_chlorophyll),
    dist_chl_mean = mean(dist_chlorophyll),
    dist_car_min = min(dist_carotenoid),
    dist_car_mean = mean(dist_carotenoid),
    dist_pheo_min = min(dist_pheophytin),
    dist_pheo_mean = mean(dist_pheophytin)
  ), by = .(gene, aln_pos)]
}))

# Merge with combined data
combined_psII <- merge(combined_psII, pigment_distances,
                       by.x = c("Gene", "Position"),
                       by.y = c("gene", "aln_pos"),
                       all.x = TRUE)

# ---- HYPOTHESIS TESTS: PIGMENT PROXIMITY ----

message("\n", paste(rep("=", 60), collapse = ""))
message("PSII PIGMENT PROXIMITY TESTS")
message(paste(rep("=", 60), collapse = ""))

# H: Distance to chlorophyll
message("\nDistance to nearest chlorophyll:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_psII[site_class == sc & !is.na(dist_chl_min), dist_chl_min]
  bg_d <- combined_psII[site_class == "not_sig" & !is.na(dist_chl_min), dist_chl_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H: Distance to carotenoid
message("\nDistance to nearest carotenoid:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_psII[site_class == sc & !is.na(dist_car_min), dist_car_min]
  bg_d <- combined_psII[site_class == "not_sig" & !is.na(dist_car_min), dist_car_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H: Distance to pheophytin (special pair electron acceptor)
message("\nDistance to nearest pheophytin:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_psII[site_class == sc & !is.na(dist_pheo_min), dist_pheo_min]
  bg_d <- combined_psII[site_class == "not_sig" & !is.na(dist_pheo_min), dist_pheo_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# Define "near pigment" thresholds (within 10Å is typical for direct contact)
combined_psII[, near_chlorophyll := dist_chl_min < 10]
combined_psII[, near_carotenoid := dist_car_min < 10]

# Fisher's tests
message("\nEnrichment near chlorophyll (<10Å):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_near <- sum(combined_psII$site_class == sc & combined_psII$near_chlorophyll == TRUE, na.rm = TRUE)
  sig_far <- sum(combined_psII$site_class == sc & combined_psII$near_chlorophyll == FALSE, na.rm = TRUE)
  bg_near <- sum(combined_psII$site_class == "not_sig" & combined_psII$near_chlorophyll == TRUE, na.rm = TRUE)
  bg_far <- sum(combined_psII$site_class == "not_sig" & combined_psII$near_chlorophyll == FALSE, na.rm = TRUE)
  
  if (sig_near + sig_far == 0) next
  mat <- matrix(c(sig_near, sig_far, bg_near, bg_far), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("  %s: %d/%d near (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_near, sig_near + sig_far,
                  100 * sig_near / (sig_near + sig_far),
                  100 * bg_near / (bg_near + bg_far),
                  ft$estimate, ft$p.value))
}

# ---- VISUALIZATION ----

# Boxplots
p_chl <- ggplot(combined_psII[!is.na(dist_chl_min)],
                aes(x = site_class, y = dist_chl_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to chlorophyll", y = "Distance (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_car <- ggplot(combined_psII[!is.na(dist_car_min)],
                aes(x = site_class, y = dist_car_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to carotenoid", y = "Distance (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_chl + p_car)

# Structure colored by chlorophyll distance
all_psII_coords <- rbindlist(lapply(psII_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(
    gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
    gr$merged[, .(Gene, Position, site_class)],
    by.x = c("gene", "aln_pos"),
    by.y = c("Gene", "Position"),
    all.x = TRUE
  )
}), fill = TRUE)

# Add chlorophyll distance to coords
all_psII_coords[, dist_chl := sapply(1:.N, function(i) {
  min(sqrt((x[i] - chl_coords$x)^2 + (y[i] - chl_coords$y)^2 + (z[i] - chl_coords$z)^2))
})]

p_psII_chl_xz <- ggplot() +
  geom_point(data = psII_structure$ca_df, aes(x = x, y = z), 
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = psII_structure$ligands[ligand_type %in% c("CHL", "CLA")],
             aes(x = lig_x, y = lig_z), color = "green3", size = 2, shape = 18) +
  geom_point(data = all_psII_coords[site_class == "sig_both"],
             aes(x = x, y = z), color = "darkred", size = 2.5) +
  geom_point(data = all_psII_coords[site_class == "sig_with_control"],
             aes(x = x, y = z), color = "steelblue", size = 2) +
  coord_fixed() +
  labs(title = "PSII: sig sites vs chlorophylls (green)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

print(p_psII_chl_xz)

# Distance to plastoquinone (PQN)
pqn_coords <- psII_structure$ligands[ligand_type == "DGD", .(x = lig_x, y = lig_y, z = lig_z)]
stopifnot(nrow(pqn_coords) > 0)

combined_psII[, dist_pqn := sapply(1:.N, function(i) {
  min(sqrt((x_mean[i] - pqn_coords$x)^2 + (y_mean[i] - pqn_coords$y)^2 + (z_mean[i] - pqn_coords$z)^2))
})]

# Test
message("Distance to DGD:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  wt <- wilcox.test(combined_psII[site_class == sc, dist_pqn], 
                    combined_psII[site_class == "not_sig", dist_pqn])
  message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f", sc,
                  median(combined_psII[site_class == sc, dist_pqn], na.rm = TRUE),
                  median(combined_psII[site_class == "not_sig", dist_pqn], na.rm = TRUE),
                  wt$p.value))
}

# Plots
p_pqn_box <- ggplot(combined_psII, aes(x = site_class, y = dist_pqn, fill = site_class)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig"="gray60", "sig_no_control"="gold", 
                               "sig_with_control"="steelblue", "sig_both"="darkred")) +
  labs(title = "Distance to PQN", y = "Å") + theme_classic() + theme(legend.position = "none")

p_pqn_scatter <- ggplot(combined_psII, aes(x = dist_pqn, y = -log10(P_aa_only))) +
  geom_point(aes(color = site_class), alpha = 0.6) +
  scale_color_manual(values = c("not_sig"="gray60", "sig_no_control"="gold",
                                "sig_with_control"="steelblue", "sig_both"="darkred")) +
  geom_smooth(method = "loess", color = "black") +
  labs(title = "PQN distance vs significance", x = "Distance to PQN (Å)", y = "-log10(P)") +
  theme_classic()

print(p_pqn_box + p_pqn_scatter)

p_psII_pqn <- ggplot() +
  geom_point(data = psII_structure$ca_df, aes(x = x, y = z), 
             color = "gray80", size = 0.3, alpha = 0.3) +
  geom_point(data = psII_structure$ligands[ligand_type == "DGD"],
             aes(x = lig_x, y = lig_z), color = "purple", size = 4, shape = 18) +
  geom_point(data = all_psII_coords[site_class == "sig_both"],
             aes(x = x, y = z), color = "darkred", size = 2.5) +
  geom_point(data = all_psII_coords[site_class == "sig_with_control"],
             aes(x = x, y = z), color = "steelblue", size = 2) +
  coord_fixed() +
  labs(title = "PSII: sig sites vs plastoquinone (purple)", x = "X (Å)", y = "Z (Å)") +
  theme_classic()

print(p_psII_pqn)


# Check what we have for orientation markers
message("PSII orientation markers:")
print(psII_structure$ligands[ligand_type %in% c("OEX", "FE2", "HEM"), 
                             .(ligand_type, chain, lig_x, lig_y, lig_z)])

# OEX = Mn4Ca cluster (lumenal)
# For stromal side, we can use the terminal acceptors or just the opposite Z extreme

oec_coords <- psII_structure$ligands[ligand_type == "OEX", 
                                     .(x = mean(lig_x), y = mean(lig_y), z = mean(lig_z))]

# Find stromal extreme (opposite Z from OEC)
z_range <- range(psII_structure$ca_df$z)
if (oec_coords$z < mean(z_range)) {
  # OEC is at low Z, so stroma is high Z
  stroma_z <- z_range[2]
  lumen_z <- z_range[1]
} else {
  stroma_z <- z_range[1]
  lumen_z <- z_range[2]
}

message(sprintf("PSII orientation: OEC z=%.1f, lumen~%.1f, stroma~%.1f",
                oec_coords$z, lumen_z, stroma_z))

# Simple approach: just use Z directly if OEC is at one extreme
# More sophisticated: use compute_membrane_depth with proper references

# For PSII, let's just flip Z if needed so that higher Z = more stromal
psII_flip_z <- oec_coords$z > mean(z_range)  # TRUE if we need to flip

if (psII_flip_z) {
  combined_psII[, z_membrane := -z_mean + max(z_mean, na.rm = TRUE)]
  message("Flipping Z for PSII (OEC was at high Z)")
} else {
  combined_psII[, z_membrane := z_mean]
  message("Keeping Z orientation for PSII (OEC at low Z = lumen)")
}

# Define compartments based on membrane depth
# Typical thylakoid membrane is ~4nm thick, so use ±20Å from center
z_center <- mean(combined_psII$z_membrane, na.rm = TRUE)
combined_psII[, compartment := fcase(
  z_mean < -45, "lumen",
  z_mean > -20, "stroma",
  default = "membrane"
)]

message("PSII compartment distribution:")
print(combined_psII[!is.na(compartment), .N, by = compartment])
print(combined_psII[!is.na(compartment), .N, by = .(site_class, compartment)])

# Test compartment enrichment
message("\n=== PSII COMPARTMENT ENRICHMENT ===")
for (comp in c("lumen", "membrane", "stroma")) {
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    sig_in <- sum(combined_psII$site_class == sc & combined_psII$compartment == comp, na.rm = TRUE)
    sig_out <- sum(combined_psII$site_class == sc & combined_psII$compartment != comp, na.rm = TRUE)
    bg_in <- sum(combined_psII$site_class == "not_sig" & combined_psII$compartment == comp, na.rm = TRUE)
    bg_out <- sum(combined_psII$site_class == "not_sig" & combined_psII$compartment != comp, na.rm = TRUE)
    
    if (sig_in + sig_out < 3) next
    mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
    ft <- fisher.test(mat)
    message(sprintf("  %s in %s: %d/%d (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                    sc, comp, sig_in, sig_in + sig_out,
                    100 * sig_in / (sig_in + sig_out),
                    100 * bg_in / (bg_in + bg_out),
                    ft$estimate, ft$p.value))
  }
}

# Visualize with compartment boundaries
p_psII_compartments <- ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -45, ymax = -20), 
            fill = "lightblue", alpha = 0.3) +
  annotate("text", x = min(psII_structure$ca_df$x), y = -55, label = "LUMEN", hjust = 0) +
  annotate("text", x = min(psII_structure$ca_df$x), y = -32, label = "MEMBRANE", hjust = 0) +
  annotate("text", x = min(psII_structure$ca_df$x), y = -10, label = "STROMA", hjust = 0) +
  geom_point(data = psII_structure$ca_df, aes(x = x, y = z), 
             color = "gray80", size = 0.3, alpha = 0.3) +
  geom_point(data = all_psII_coords[site_class == "not_sig"],
             aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_psII_coords[site_class == "sig_no_control"],
             aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_psII_coords[site_class == "sig_with_control"],
             aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_psII_coords[site_class == "sig_both"],
             aes(x = x, y = z), color = "darkred", size = 2.5) +
  geom_hline(yintercept = c(-45, -20), linetype = "dashed", alpha = 0.5) +
  coord_fixed() +
  labs(title = "PSII with compartment boundaries",
       subtitle = "Blue shading = membrane region",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

print(p_psII_compartments)

# Boxplot by compartment
p_psII_comp_box <- ggplot(combined_psII[!is.na(compartment)],
                          aes(x = compartment, fill = site_class)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "PSII: Site class by compartment", y = "Proportion", x = "") +
  theme_classic()

print(p_psII_comp_box)

# ---- PSI (4XK8) ----
message("\n--- Photosystem I (4XK8) ---")
psI_structure <- load_pdb_structure("5ZJI")
print(psI_structure$chain_summary)



#this has another oligomer attatched to it, we should drop that 

#psI_structure$ca_df <- psI_structure$ca_df[chain == toupper(chain)]

message("\nLigands:")
print(psI_structure$ligands[, .N, by = ligand_type])

psI_def <- DEFAULT_COMPLEXES$photosystem_i

psI_result <- analyze_complex(
  "photosystem_i",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS
)

# Combine test results across all genes
psI_combined_tests <- combine_test_results(psI_result)

message("\n=== COMBINED psI SYNTHASE RESULTS ===\n")

message("Continuous features (significant results only):")
sig_continuous <- psI_combined_tests$continuous[wilcox_p < 0.05]
if (nrow(sig_continuous) > 0) {
  print(sig_continuous[order(wilcox_p), .(gene, feature, site_class, 
                                          median_sig, median_bg, direction, wilcox_p)])
} else {
  message("  No significant results")
}

message("\nBinary features (significant results only):")
sig_binary <- psI_combined_tests$binary[fisher_p < 0.05]
if (nrow(sig_binary) > 0) {
  print(sig_binary[order(fisher_p), .(gene, feature, site_class,
                                      pct_sig, pct_bg, odds_ratio, fisher_p)])
} else {
  message("  No significant results")
}


# --- Hypothesis testing on combined data ---

message("\n=== HYPOTHESIS TESTS ON COMBINED psI DATA ===\n")

combined_psI <- psI_result$combined_merged

# H1: Are sig sites more stromal (higher Z)?
message("H1: Stromal enrichment")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- combined_psI[site_class == sc & !is.na(z_mean), z_mean]
  bg_z <- combined_psI[site_class == "not_sig" & !is.na(z_mean), z_mean]
  if (length(test_z) >= 3) {
    wt <- wilcox.test(test_z, bg_z)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f", 
                    sc, median(test_z), median(bg_z), wt$p.value))
  }
}



# H3: Do sig sites avoid ligand binding sites?
message("\nH3: Ligand avoidance")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_psI[site_class == sc & !is.na(dist_ligand_min), dist_ligand_min]
  bg_d <- combined_psI[site_class == "not_sig" & !is.na(dist_ligand_min), dist_ligand_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    message(sprintf("  %s: median=%.1f vs %.1f, p=%.4f",
                    sc, median(test_d), median(bg_d), wt$p.value))
  }
}

# Combined visualization
p_psI_all_xz <- plot_complex_gwas(psI_result, "xz")
p_psI_all_xy <- plot_complex_gwas(psI_result, "xy")

p_psI_all_xz + p_psI_all_xy

# ---- psI mem orientation ----
# Chain identification by size:
# A (742 aa) = psaA
# B (733 aa) = psaB  
# C (81 aa) = psaC (stromal, ferredoxin docking)
# Small PSI subunits: D, E, F, G, H, I, J, K, L, N, O
# LHCI subunits: 1, 2, 3, 4, X, Y, Z (numbered/lettered ~200-230 aa)

# Let's check the orientation - plot XY colored by chain type
psI_structure$ca_df[, subcomplex := fifelse(chain %in% c("A", "B"), "PSI_core",
                                            fifelse(chain %in% c("1","2","3","4","X","Y","Z"), "LHCI", "PSI_small"))]

p_psI_xy <- ggplot(psI_structure$ca_df, aes(x = x, y = y, color = subcomplex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() + labs(title = "PSI-LHCI - XY view") + theme_classic()

p_psI_xz <- ggplot(psI_structure$ca_df, aes(x = x, y = z, color = subcomplex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() + labs(title = "PSI-LHCI - XZ view") + theme_classic()

p_psI_yz <- ggplot(psI_structure$ca_df, aes(x = y, y = z, color = subcomplex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() + labs(title = "PSI-LHCI - YZ view") + theme_classic()

print(p_psI_xy + p_psI_xz + p_psI_yz)

# Check where psaC is (stromal marker - ferredoxin binding)
psaC_center <- psI_structure$ca_df[chain == "C", .(x = mean(x), y = mean(y), z = mean(z))]
message("psaC (stromal) center: x=", round(psaC_center$x, 1), 
        " y=", round(psaC_center$y, 1), " z=", round(psaC_center$z, 1))

# Check SF4 clusters (iron-sulfur centers in psaC - stromal side)
print(psI_structure$ligands[ligand_type == "SF4", .(ligand_type, lig_x, lig_y, lig_z)])

# P700 special pair chlorophylls should be in membrane
# Plastocyanin binds on lumen side

# Rotate around Z axis to align membrane normal with Y
# Clockwise rotation by theta means: x' = x*cos(theta) + y*sin(theta)
#                                    y' = -x*sin(theta) + y*cos(theta)

theta <- 30 * pi / 180  # -35 degrees (clockwise)

psI_structure$ca_df[, `:=`(
  x_rot = x * cos(theta) + y * sin(theta),
  y_rot = -x * sin(theta) + y * cos(theta)
)]

# Also rotate ligands
psI_structure$ligands[, `:=`(
  lig_x_rot = lig_x * cos(theta) + lig_y * sin(theta),
  lig_y_rot = -lig_x * sin(theta) + lig_y * cos(theta)
)]

# Check the rotation - SF4 should now have similar y_rot (stromal side)
message("SF4 after rotation:")
print(psI_structure$ligands[ligand_type == "SF4", .(lig_x_rot, lig_y_rot, lig_z)])

# Visualize rotated structure
p_rot_xy <- ggplot(psI_structure$ca_df, aes(x = x_rot, y = y_rot, color = subcomplex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() + labs(title = "PSI - XY rotated") + theme_classic()

p_rot_xz <- ggplot(psI_structure$ca_df, aes(x = x_rot, y = z, color = subcomplex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  coord_fixed() + labs(title = "PSI - XZ rotated (side view)") + theme_classic()

print(p_rot_xy + p_rot_xz)

# Now Y_rot should be membrane normal - check range
message("\nY_rot range: ", round(min(psI_structure$ca_df$y_rot), 1), " to ", 
        round(max(psI_structure$ca_df$y_rot), 1))
message("SF4 y_rot (should be stromal extreme): ", 
        round(mean(psI_structure$ligands[ligand_type == "SF4", lig_y_rot]), 1))

# psaF (chain F) has lumenal domain for plastocyanin binding
# psaC (chain C) is stromal (we know SF4 is there)

psaF_center <- psI_structure$ca_df[chain == "F", .(x_rot = mean(x_rot), y_rot = mean(y_rot), z = mean(z))]
psaC_center <- psI_structure$ca_df[chain == "C", .(x_rot = mean(x_rot), y_rot = mean(y_rot), z = mean(z))]

message("psaC (stromal): y_rot = ", round(psaC_center$y_rot, 1))
message("psaF (lumenal): y_rot = ", round(psaF_center$y_rot, 1))

# Plot with psaC and psaF highlighted
p_check <- ggplot() +
  geom_point(data = psI_structure$ca_df, aes(x = x_rot, y = y_rot), 
             color = "gray70", size = 0.3, alpha = 0.3) +
  geom_point(data = psI_structure$ca_df[chain == "C"], aes(x = x_rot, y = y_rot), 
             color = "red", size = 1.5) +
  geom_point(data = psI_structure$ca_df[chain == "F"], aes(x = x_rot, y = y_rot), 
             color = "blue", size = 1.5) +
  geom_point(data = psI_structure$ligands[ligand_type == "SF4"], 
             aes(x = lig_x_rot, y = lig_y_rot), color = "red", size = 4, shape = 18) +
  annotate("text", x = psaC_center$x_rot, y = psaC_center$y_rot + 5, label = "psaC (stroma)", color = "red") +
  annotate("text", x = psaF_center$x_rot, y = psaF_center$y_rot - 5, label = "psaF (lumen)", color = "blue") +
  coord_fixed() +
  labs(title = "PSI orientation check", subtitle = "Red=stromal (psaC/SF4), Blue=lumenal (psaF)") +
  theme_classic()

print(p_check)

# Add rotated coordinates to position_maps for each gene
for (gene_name in names(psI_result$gene_results)) {
  gr <- psI_result$gene_results[[gene_name]]
  if (is.null(gr$position_maps)) next
  
  # Apply same rotation to position_maps
  theta <- 30 * pi / 180
  gr$position_maps[, `:=`(
    x_rot = x * cos(theta) + y * sin(theta),
    y_rot = -x * sin(theta) + y * cos(theta),
    y_membrane = -(-x * sin(theta) + y * cos(theta))  # flipped
  )]
}

# Aggregate y_membrane into combined_psI
membrane_coords <- rbindlist(lapply(psI_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  gr$position_maps[, .(y_membrane_mean = mean(y_membrane)), by = .(gene, aln_pos)]
}))

combined_psI <- merge(combined_psI, membrane_coords,
                      by.x = c("Gene", "Position"),
                      by.y = c("gene", "aln_pos"),
                      all.x = TRUE)

# Set compartment boundaries - eyeball from your plot
# Membrane core is roughly y_membrane = -20 to +20?
message("y_membrane range in GWAS data:")
print(summary(combined_psI$y_membrane_mean))

# Set boundaries (adjust after seeing summary)
combined_psI[, compartment := fcase(
  y_membrane_mean < -20, "lumen",
  y_membrane_mean > 20, "stroma",
  default = "membrane"
)]

message("\nPSI compartment distribution:")
print(combined_psI[!is.na(compartment), .N, by = .(site_class, compartment)])

# Compartment enrichment tests
message("\n=== PSI COMPARTMENT ENRICHMENT ===")
for (comp in c("lumen", "membrane", "stroma")) {
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    sig_in <- sum(combined_psI$site_class == sc & combined_psI$compartment == comp, na.rm = TRUE)
    sig_out <- sum(combined_psI$site_class == sc & combined_psI$compartment != comp, na.rm = TRUE)
    bg_in <- sum(combined_psI$site_class == "not_sig" & combined_psI$compartment == comp, na.rm = TRUE)
    bg_out <- sum(combined_psI$site_class == "not_sig" & combined_psI$compartment != comp, na.rm = TRUE)
    
    if (sig_in + sig_out < 3) next
    mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
    ft <- fisher.test(mat)
    message(sprintf("  %s in %s: %d/%d (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                    sc, comp, sig_in, sig_in + sig_out,
                    100 * sig_in / (sig_in + sig_out),
                    100 * bg_in / (bg_in + bg_out),
                    ft$estimate, ft$p.value))
  }
}

# Visualize with compartment boundaries and GWAS hits
all_psI_coords <- rbindlist(lapply(psI_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x_rot, y_membrane)],
        gr$merged[, .(Gene, Position, site_class)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)

all_psI_coords <- rbindlist(lapply(psI_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x_rot, y_membrane)],
        gr$merged[, .(Gene, Position, site_class)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)

# Check it worked
stopifnot("y_membrane" %in% names(all_psI_coords))

# Now plot
p_psI_compartments <- ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -20, ymax = 20), fill = "lightblue", alpha = 0.3) +
  geom_point(data = psI_structure$ca_df, aes(x = x_rot, y = y_membrane), 
             color = "gray80", size = 0.3, alpha = 0.3) +
  geom_point(data = all_psI_coords[site_class == "not_sig"], aes(x = x_rot, y = y_membrane), 
             color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_psI_coords[site_class == "sig_no_control"], aes(x = x_rot, y = y_membrane), 
             color = "gold", size = 2) +
  geom_point(data = all_psI_coords[site_class == "sig_with_control"], aes(x = x_rot, y = y_membrane), 
             color = "steelblue", size = 2) +
  geom_point(data = all_psI_coords[site_class == "sig_both"], aes(x = x_rot, y = y_membrane), 
             color = "darkred", size = 2.5) +
  geom_hline(yintercept = c(-30, 10), linetype = "dashed") +
  annotate("text", x = -120, y = -35, label = "LUMEN") +
  annotate("text", x = -120, y = 0, label = "MEMBRANE") +
  annotate("text", x = -120, y = 35, label = "STROMA") +
  coord_fixed() +
  labs(title = "PSI compartment enrichment", x = "X (Å)", y = "Y membrane (Å)") +
  theme_classic()

print(p_psI_compartments)

combined_psI[, compartment := fcase(
  y_membrane_mean < -20, "lumen",
  y_membrane_mean > 20, "stroma",
  default = "membrane"
)]

download_opm <- function(pdb_id, outdir = "data/pdb_cache") {
  url <- sprintf("https://opm-assets.storage.googleapis.com/pdb/%s.pdb", tolower(pdb_id))
  outfile <- file.path(outdir, paste0(pdb_id, "_opm.pdb"))
  download.file(url, outfile)
  message("Downloaded OPM structure: ", outfile)
  outfile
}

download_opm("5ZJI")
psI_opm <- load_pdb_structure("data/pdb_cache/5ZJI_opm.pdb")


# Formal test
message("\n=== PSI COMPARTMENT ENRICHMENT ===")
for (comp in c("lumen", "membrane", "stroma")) {
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    sig_in <- sum(combined_psI$site_class == sc & combined_psI$compartment == comp, na.rm = TRUE)
    sig_out <- sum(combined_psI$site_class == sc & combined_psI$compartment != comp, na.rm = TRUE)
    bg_in <- sum(combined_psI$site_class == "not_sig" & combined_psI$compartment == comp, na.rm = TRUE)
    bg_out <- sum(combined_psI$site_class == "not_sig" & combined_psI$compartment != comp, na.rm = TRUE)
    
    if (sig_in + sig_out < 3) next
    mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
    ft <- fisher.test(mat)
    message(sprintf("  %s in %s: %d/%d (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                    sc, comp, sig_in, sig_in + sig_out,
                    100 * sig_in / (sig_in + sig_out),
                    100 * bg_in / (bg_in + bg_out),
                    ft$estimate, ft$p.value))
  }
}

psI_structure$ca_df$neighbor_count <- with(psI_structure$ca_df, 
                                            compute_neighbor_count(x, y, z))

# Add neighbor count (surface proxy) to combined data
# Need to compute per-chain, then aggregate

# Compute for all CA atoms in structure
# Check distribution
message("Neighbor count distribution (all CA atoms):")
print(summary(psII_structure$ca_df$neighbor_count))
hist(psI_structure$ca_df$neighbor_count)

# Now add to position_maps for each gene and aggregate
for (gene_name in names(psI_result$gene_results)) {
  gr <- psI_result$gene_results[[gene_name]]
  if (is.null(gr)) next
  
  # Get neighbor counts for residues in this gene's chains
  gr$position_maps <- merge(
    gr$position_maps,
    psI_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
}

# Aggregate across chains and add to combined_merged
# Aggregate directly without modifying atp_results
surface_features <- rbindlist(lapply(names(psI_result$gene_results), function(gene_name) {
  gr <- psI_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  # Merge here instead
  pm <- merge(
    gr$position_maps,
    psI_structure$ca_df[, .(chain, resno, neighbor_count)],
    by.x = c("chain", "pdb_resno"),
    by.y = c("chain", "resno"),
    all.x = TRUE
  )
  
  pm[, .(
    neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
    neighbor_count_min = min(neighbor_count, na.rm = TRUE),
    neighbor_count_max = max(neighbor_count, na.rm = TRUE)
  ), by = .(gene, aln_pos)]
}))

# Merge with combined data
combined_psI <- merge(combined_psI, surface_features,
                       by.x = c("Gene", "Position"),
                       by.y = c("gene", "aln_pos"),
                       all.x = TRUE)

# Define "surface exposed" as below median neighbor count
median_neighbors <- median(combined_psI$neighbor_count_min, na.rm = TRUE)
combined_psI[, surface_exposed := neighbor_count_min < median_neighbors]

message("\nNeighbor count by site class (lower = more surface exposed):")
print(combined_psI[!is.na(neighbor_count_min), .(
  n = .N,
  median_neighbors = as.numeric(median(neighbor_count_min)),
  pct_surface = 100 * mean(surface_exposed, na.rm = TRUE)
), by = site_class])


boxplot(neighbor_count_min ~ site_class, combined_psI)

# ---- ndh ----

# NDH COMPLEX ANALYSIS - 7EU3 (Hordeum vulgare)

message("\n", paste(rep("#", 70), collapse = ""))
message("NDH COMPLEX ANALYSIS")
message(paste(rep("#", 70), collapse = ""))

# Load structure
ndh_structure <- load_pdb_structure("7EU3")
print(ndh_structure$chain_summary)

message("\nLigands:")
print(ndh_structure$ligands[, .N, by = ligand_type])

# Define complex
ndh_def <- list(
  pdb = "7EU3",
  organism = "Hordeum vulgare",
  description = "NDH-1 complex",
  chain_map = list(
    ndhA = "A",
    ndhB = "B",
    ndhC = "C",
    ndhD = "D",
    ndhG = "G",
    ndhH = "H",
    ndhI = "I",
    ndhJ = "J",
    ndhK = "K"
  ),
  compartment_z = NULL,
  compartment_labels = NULL
)


# Quick structure visualization
p_ndh_xy <- ggplot(ndh_structure$ca_df, aes(x = x, y = y, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XY view") + theme_classic()

p_ndh_xz <- ggplot(ndh_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XZ view") + theme_classic()

p_ndh_yz <- ggplot(ndh_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - YZ view") + theme_classic()

print(p_ndh_xy + p_ndh_xz + p_ndh_yz)

# ---- RUN MAIN ANALYSIS ----
ndh_result <- analyze_complex(
  "ndh_complex",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS,
  genes_to_analyze = names(ndh_def$chain_map)
)

# ---- SURFACE ACCESSIBILITY ----
ndh_structure$ca_df$neighbor_count <- with(ndh_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

surface_features <- rbindlist(lapply(names(ndh_result$gene_results), function(gene_name) {
  gr <- ndh_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  pm <- merge(gr$position_maps,
              ndh_structure$ca_df[, .(chain, resno, neighbor_count)],
              by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
  
  pm[, .(neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
         neighbor_count_min = min(neighbor_count, na.rm = TRUE)), by = .(gene, aln_pos)]
}))

combined_ndh <- ndh_result$combined_merged
combined_ndh <- merge(combined_ndh, surface_features,
                      by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)

median_neighbors <- median(combined_ndh$neighbor_count_min, na.rm = TRUE)
combined_ndh[, surface_exposed := neighbor_count_min < median_neighbors]

# ---- LIGAND DISTANCES ----
# NDH has Fe-S clusters, quinone binding sites
message("\nLigand types present:")
print(ndh_structure$ligands[, .N, by = ligand_type])

# Iron-sulfur clusters (electron transfer chain)
fes_types <- c("SF4", "FES", "F3S")  # different Fe-S nomenclature
fes_coords <- ndh_structure$ligands[ligand_type %in% fes_types, .(x = lig_x, y = lig_y, z = lig_z)]

# Plastoquinone
pq_coords <- ndh_structure$ligands[ligand_type %in% c("PL9", "PLQ", "PQN"), .(x = lig_x, y = lig_y, z = lig_z)]

if (nrow(fes_coords) > 0) {
  combined_ndh[, dist_fes := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - fes_coords$x)^2 + (y_mean[i] - fes_coords$y)^2 + (z_mean[i] - fes_coords$z)^2))
  })]
  message("Added Fe-S distance (", nrow(fes_coords), " clusters)")
}

if (nrow(pq_coords) > 0) {
  combined_ndh[, dist_pq := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - pq_coords$x)^2 + (y_mean[i] - pq_coords$y)^2 + (z_mean[i] - pq_coords$z)^2))
  })]
  message("Added plastoquinone distance (", nrow(pq_coords), " sites)")
}

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("NDH HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

message("\nSite class distribution:")
print(combined_ndh[, .N, by = site_class])

# H1: Surface exposure
message("\nH1: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_ndh[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_ndh[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# H2: Interface proximity
message("\nH2: Interface proximity")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_ndh[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_ndh[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: Fe-S cluster proximity (if present)
if ("dist_fes" %in% names(combined_ndh)) {
  message("\nH3: Distance to Fe-S clusters")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_ndh[site_class == sc & !is.na(dist_fes), dist_fes]
    bg_d <- combined_ndh[site_class == "not_sig" & !is.na(dist_fes), dist_fes]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# H4: Plastoquinone proximity (if present)
if ("dist_pq" %in% names(combined_ndh)) {
  message("\nH4: Distance to plastoquinone")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_ndh[site_class == sc & !is.na(dist_pq), dist_pq]
    bg_d <- combined_ndh[site_class == "not_sig" & !is.na(dist_pq), dist_pq]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# ---- VISUALIZATIONS ----

# Build coordinate table for plotting
all_ndh_coords <- rbindlist(lapply(ndh_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
        gr$merged[, .(Gene, Position, site_class, P_aa_only)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)
all_ndh_coords[, neglog_p := -log10(P_aa_only)]

# Structure with GWAS hits - all 3 projections
p_ndh_gwas_xy <- ggplot() +
  geom_point(data = ndh_structure$ca_df, aes(x = x, y = y), color = "gray80", size = 0.3, alpha = 0.3) +
  geom_point(data = all_ndh_coords[site_class == "not_sig"], aes(x = x, y = y), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_ndh_coords[site_class == "sig_no_control"], aes(x = x, y = y), color = "gold", size = 2) +
  geom_point(data = all_ndh_coords[site_class == "sig_with_control"], aes(x = x, y = y), color = "steelblue", size = 2) +
  geom_point(data = all_ndh_coords[site_class == "sig_both"], aes(x = x, y = y), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "NDH - XY (top)") + theme_classic()

p_ndh_gwas_xz <- ggplot() +
  geom_point(data = ndh_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
  geom_point(data = all_ndh_coords[site_class == "not_sig"], aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_ndh_coords[site_class == "sig_no_control"], aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_ndh_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_ndh_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "NDH - XZ (side)") + theme_classic()

print(p_ndh_gwas_xy + p_ndh_gwas_xz)

# Color by -log10(P)
p_ndh_pval <- ggplot() +
  geom_point(data = ndh_structure$ca_df, aes(x = x, y = z), color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = all_ndh_coords[!is.na(neglog_p)], aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() + labs(title = "NDH colored by significance") + theme_classic()

print(p_ndh_pval)

# Boxplots
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_box_surf <- ggplot(combined_ndh[!is.na(neighbor_count_min)],
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure", y = "Neighbor count (lower=exposed)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_box_int <- ggplot(combined_ndh[!is.na(dist_interface_min)],
                    aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Interface proximity", y = "Distance to interface (Å)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_surf + p_box_int)

# Ligand-specific plots if available
if ("dist_fes" %in% names(combined_ndh)) {
  p_box_fes <- ggplot(combined_ndh[!is.na(dist_fes)],
                      aes(x = site_class, y = dist_fes, fill = site_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = fill_vals) +
    labs(title = "Fe-S cluster proximity", y = "Distance (Å)") +
    theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Structure with Fe-S clusters highlighted
  p_ndh_fes <- ggplot() +
    geom_point(data = ndh_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = ndh_structure$ligands[ligand_type %in% fes_types],
               aes(x = lig_x, y = lig_z), color = "orange", size = 4, shape = 18) +
    geom_point(data = all_ndh_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
    geom_point(data = all_ndh_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
    coord_fixed() + labs(title = "NDH: sig sites vs Fe-S clusters (orange)") + theme_classic()
  
  print(p_box_fes + p_ndh_fes)
}

# ---- PER-GENE BREAKDOWN ----
message("\n=== SITE COUNTS BY GENE ===")
print(combined_ndh[, .(
  total = .N,
  with_structure = sum(!is.na(x_mean)),
  sig_no_ctrl = sum(site_class == "sig_no_control"),
  sig_with_ctrl = sum(site_class == "sig_with_control"),
  sig_both = sum(site_class == "sig_both")
), by = Gene][order(-sig_both)])

# ---- NDH-SPECIFIC: MEMBRANE ARM vs HYDROPHILIC ARM ----
# NDH has an L-shaped structure: membrane arm (ndhA-G) and hydrophilic arm (ndhH-K)
# Test if signals are enriched in one arm

combined_ndh[, ndh_arm := fifelse(Gene %in% c("ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG"),
                                  "membrane_arm", "hydrophilic_arm")]

message("\n=== NDH ARM ENRICHMENT ===")
for (arm in c("membrane_arm", "hydrophilic_arm")) {
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    sig_in <- sum(combined_ndh$site_class == sc & combined_ndh$ndh_arm == arm, na.rm = TRUE)
    sig_out <- sum(combined_ndh$site_class == sc & combined_ndh$ndh_arm != arm, na.rm = TRUE)
    bg_in <- sum(combined_ndh$site_class == "not_sig" & combined_ndh$ndh_arm == arm, na.rm = TRUE)
    bg_out <- sum(combined_ndh$site_class == "not_sig" & combined_ndh$ndh_arm != arm, na.rm = TRUE)
    
    if (sig_in + sig_out < 3) next
    mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
    ft <- fisher.test(mat)
    message(sprintf("  %s in %s: %d/%d (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                    sc, arm, sig_in, sig_in + sig_out,
                    100 * sig_in / (sig_in + sig_out),
                    100 * bg_in / (bg_in + bg_out),
                    ft$estimate, ft$p.value))
  }
}

# Visualize arms
ndh_structure$ca_df[, ndh_arm := fifelse(chain %in% c("A", "B", "C", "D", "E", "F", "G"),
                                         "membrane_arm", "hydrophilic_arm")]

p_ndh_arms <- ggplot(ndh_structure$ca_df, aes(x = x, y = z, color = ndh_arm)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("membrane_arm" = "coral", "hydrophilic_arm" = "turquoise")) +
  coord_fixed() + labs(title = "NDH arms") + theme_classic()

print(p_ndh_arms)

# Save results
saveRDS(ndh_result, "results/ndh_3d_analysis.rds")
saveRDS(combined_ndh, "results/ndh_merged_with_features.rds")


# What ligands do we have?
print(ndh_structure$ligands[, .N, by = ligand_type])

# Visualize surface exposure pattern on structure
all_ndh_coords <- merge(all_ndh_coords,
                        ndh_structure$ca_df[, .(chain, resno, neighbor_count)],
                        by.x = c("chain", "x", "y", "z"),
                        by.y = c("chain", "x", "y", "z"),
                        all.x = TRUE)

# Surface map with sig sites
p_ndh_surface <- ggplot() +
  geom_point(data = ndh_structure$ca_df, aes(x = x, y = z, color = neighbor_count),
             size = 0.8, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "Neighbors\n(low=exposed)") +
  geom_point(data = all_ndh_coords[site_class == "sig_both"], aes(x = x, y = z),
             color = "red", size = 3, shape = 1, stroke = 1.5) +
  coord_fixed() +
  labs(title = "NDH: sig_both sites (red circles) on surface map") +
  theme_classic()

print(p_ndh_surface)

# Scatter: neighbor count vs p-value
p_ndh_scatter <- ggplot(combined_ndh[!is.na(neighbor_count_min)],
                        aes(x = neighbor_count_min, y = -log10(P_aa_only), color = site_class)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("not_sig" = "gray50", "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  geom_smooth(aes(group = 1), method = "loess", color = "black", se = TRUE) +
  labs(title = "NDH: Surface exposure vs significance",
       x = "Neighbor count (lower = more exposed)", y = "-log10(P)") +
  theme_classic()

print(p_ndh_scatter)

# Correlation
cor_surf <- cor.test(combined_ndh$neighbor_count_min, -log10(combined_ndh$P_aa_only), 
                     method = "spearman", use = "complete")
message(sprintf("Spearman correlation (neighbors vs -log10P): rho = %.3f, p = %.2e",
                cor_surf$estimate, cor_surf$p.value))


# ---- rps: ribosome ssu ---- 
rps_structure <- load_pdb_structure("5MMJ")
print(rps_structure$chain_summary)

message("\nLigands:")
print(rps_structure$ligands[, .N, by = ligand_type])

# Quick structure visualization
p_rps_xy <- ggplot(rps_structure$ca_df, aes(x = x, y = y, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XY view") + theme_classic()

p_rps_xz <- ggplot(rps_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XZ view") + theme_classic()

p_rps_yz <- ggplot(rps_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - YZ view") + theme_classic()

print(p_rps_xy + p_rps_xz + p_rps_yz)
rps_def <- DEFAULT_COMPLEXES$rps_complex

rps_result <- analyze_complex(
  "rps_complex",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS,
  genes_to_analyze = names(rps_def$chain_map)
)

rps_result$gene_results

rps_structure$ca_df$neighbor_count <- with(rps_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

hist(rps_structure$ca_df$neighbor_count)
surface_features <- rbindlist(lapply(names(rps_result$gene_results), function(gene_name) {
  gr <- rps_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  pm <- merge(gr$position_maps,
              rps_structure$ca_df[, .(chain, resno, neighbor_count)],
              by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
  
  pm[, .(neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
         neighbor_count_min = min(neighbor_count, na.rm = TRUE)), by = .(gene, aln_pos)]
}))

combined_rps <- rps_result$combined_merged
combined_rps <- merge(combined_rps, surface_features,
                      by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)

median_neighbors <- median(combined_rps$neighbor_count_min, na.rm = TRUE)
combined_rps[, surface_exposed := neighbor_count_min < median_neighbors]

# ---- LIGAND DISTANCES ----
# rps has Fe-S clusters, quinone binding sites
message("\nLigand types present:")
print(rps_structure$ligands[, .N, by = ligand_type])

# Iron-sulfur clusters (electron transfer chain)
rna_nucs <- c("U", "G", "A", "C")  # different Fe-S nomenclature
rna_coords <- rps_structure$ligands[ligand_type %in% rna_nucs, .(x = lig_x, y = lig_y, z = lig_z)]

# Plastoquinone
mg_coords <- rps_structure$ligands[ligand_type %in% c("MG"), .(x = lig_x, y = lig_y, z = lig_z)]

if (nrow(rna_coords) > 0) {
  combined_rps[, dist_rna := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - rna_coords$x)^2 + (y_mean[i] - rna_coords$y)^2 + (z_mean[i] - rna_coords$z)^2))
  })]
  message("Added rna distance (", nrow(rna_coords), " clusters)")
}

if (nrow(mg_coords) > 0) {
  combined_rps[, dist_mg := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - mg_coords$x)^2 + (y_mean[i] - mg_coords$y)^2 + (z_mean[i] - mg_coords$z)^2))
  })]
  message("Added mg distance (", nrow(mg_coords), " sites)")
}

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("rps HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

message("\nSite class distribution:")
print(combined_rps[, .N, by = site_class])

# H1: Surface exposure
message("\nH1: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_rps[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_rps[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# H2: Interface proximity
message("\nH2: Interface proximity")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_rps[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_rps[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: Fe-S cluster proximity (if present)
if ("dist_rna" %in% names(combined_rps)) {
  message("\nH3: Distance to rna")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_rps[site_class == sc & !is.na(dist_rna), dist_rna]
    bg_d <- combined_rps[site_class == "not_sig" & !is.na(dist_rna), dist_rna]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# H4: Plastoquinone proximity (if present)
if ("dist_mg" %in% names(combined_rps)) {
  message("\nH4: Distance to mg")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_rps[site_class == sc & !is.na(dist_mg), dist_mg]
    bg_d <- combined_rps[site_class == "not_sig" & !is.na(dist_mg), dist_mg]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# ---- VISUALIZATIONS ----

# Build coordinate table for plotting
all_rps_coords <- rbindlist(lapply(rps_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
        gr$merged[, .(Gene, Position, site_class, P_aa_only)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)
all_rps_coords[, neglog_p := -log10(P_aa_only)]

# Structure with GWAS hits - all 3 projections
p_rps_gwas_xy <- ggplot() +
  geom_point(data = rps_structure$ca_df, aes(x = x, y = y), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rps_coords[site_class == "not_sig"], aes(x = x, y = y), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rps_coords[site_class == "sig_no_control"], aes(x = x, y = y), color = "gold", size = 2) +
  geom_point(data = all_rps_coords[site_class == "sig_with_control"], aes(x = x, y = y), color = "steelblue", size = 2) +
  geom_point(data = all_rps_coords[site_class == "sig_both"], aes(x = x, y = y), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rps - XY (top)") + theme_classic()

p_rps_gwas_xz <- ggplot() +
  geom_point(data = rps_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rps_coords[site_class == "not_sig"], aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rps_coords[site_class == "sig_no_control"], aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_rps_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_rps_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rps - XZ (side)") + theme_classic()

print(p_rps_gwas_xy + p_rps_gwas_xz)

# Color by -log10(P)
p_rps_pval <- ggplot() +
  geom_point(data = rps_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rps_coords[!is.na(neglog_p)], aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() + labs(title = "rps colored by significance") + theme_classic()

print(p_rps_pval)

# Boxplots
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_box_surf <- ggplot(combined_rps[!is.na(neighbor_count_min)],
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure", y = "Neighbor count (lower=exposed)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_box_int <- ggplot(combined_rps[!is.na(dist_interface_min)],
                    aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Interface proximity", y = "Distance to interface (Å)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_surf + p_box_int)

# Ligand-specific plots if available
if ("dist_rna" %in% names(combined_rps)) {
  p_box_rna <- ggplot(combined_rps[!is.na(dist_rna)],
                      aes(x = site_class, y = dist_rna, fill = site_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = fill_vals) +
    labs(title = "RNA proximity", y = "Distance (Å)") +
    theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Structure with Fe-S clusters highlighted
  p_rps_rna <- ggplot() +
    geom_point(data = rps_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = rps_structure$ligands[ligand_type %in% rna_nucs],
               aes(x = lig_x, y = lig_z), color = "orange", size = 4, shape = 18) +
    geom_point(data = all_rps_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
    geom_point(data = all_rps_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
    coord_fixed() + labs(title = "rps: sig sites vs rna (orange)") + theme_classic()
  
  print(p_box_rna + p_rps_rna)
}

plot(combined_rps$neighbor_count_min, -log10(combined_rps$P_aa_only))
plot(combined_rps$neighbor_count_min, -log10(combined_rps$P_aa_only))

df <- combined_rps
df$neglog10p <- -log10(df$P_aa_only)

ct <- cor.test(df$neighbor_count_min, df$neglog10p, method = "spearman")

p_overall <-
  ggplot(df, aes(neighbor_count_min, neglog10p)) +
  geom_point(color = "gray40", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("rho = %.3f\np = %.2g", ct$estimate, ct$p.value),
    hjust = 1.1, vjust = 1.2
  ) +
  theme_classic()

p_overall

p_by_class <-
  ggplot(df, aes(neighbor_count_min, neglog10p, color = site_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

p_by_class

# ---- rpl ---- 

rpl_structure <- load_pdb_structure("5H1S")
print(rpl_structure$chain_summary)

message("\nLigands:")
print(rpl_structure$ligands[, .N, by = ligand_type])

# Quick structure visualization
p_rpl_xy <- ggplot(rpl_structure$ca_df, aes(x = x, y = y, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XY view") + theme_classic()

p_rpl_xz <- ggplot(rpl_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XZ view") + theme_classic()

p_rpl_yz <- ggplot(rpl_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - YZ view") + theme_classic()

print(p_rpl_xy + p_rpl_xz + p_rpl_yz)
rpl_def <- DEFAULT_COMPLEXES$rpl_complex

rpl_result <- analyze_complex(
  "rpl_complex",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS,
  genes_to_analyze = names(rpl_def$chain_map)
)

rpl_result$gene_results

rpl_structure$ca_df$neighbor_count <- with(rpl_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

hist(rpl_structure$ca_df$neighbor_count)
surface_features <- rbindlist(lapply(names(rpl_result$gene_results), function(gene_name) {
  gr <- rpl_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  pm <- merge(gr$position_maps,
              rpl_structure$ca_df[, .(chain, resno, neighbor_count)],
              by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
  
  pm[, .(neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
         neighbor_count_min = min(neighbor_count, na.rm = TRUE)), by = .(gene, aln_pos)]
}))

combined_rpl <- rpl_result$combined_merged
combined_rpl <- merge(combined_rpl, surface_features,
                      by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)

median_neighbors <- median(combined_rpl$neighbor_count_min, na.rm = TRUE)
combined_rpl[, surface_exposed := neighbor_count_min < median_neighbors]

# ---- LIGAND DISTANCES ----
# rpl has Fe-S clusters, quinone binding sites
message("\nLigand types present:")
print(rpl_structure$ligands[, .N, by = ligand_type])

# Iron-sulfur clusters (electron transfer chain)
rna_nucs <- c("U", "G", "A", "C")  # different Fe-S nomenclature
rna_coords <- rpl_structure$ligands[ligand_type %in% rna_nucs, .(x = lig_x, y = lig_y, z = lig_z)]

if (nrow(rna_coords) > 0) {
  combined_rpl[, dist_rna := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - rna_coords$x)^2 + (y_mean[i] - rna_coords$y)^2 + (z_mean[i] - rna_coords$z)^2))
  })]
  message("Added rna distance (", nrow(rna_coords), " clusters)")
}

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("rpl HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

message("\nSite class distribution:")
print(combined_rpl[, .N, by = site_class])

# H1: Surface exposure
message("\nH1: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_rpl[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_rpl[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# H2: Interface proximity
message("\nH2: Interface proximity")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_rpl[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_rpl[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: rma cluster proximity (if present)
if ("dist_rna" %in% names(combined_rpl)) {
  message("\nH3: Distance to rna")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_rpl[site_class == sc & !is.na(dist_rna), dist_rna]
    bg_d <- combined_rpl[site_class == "not_sig" & !is.na(dist_rna), dist_rna]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# ---- VISUALIZATIONS ----

# Build coordinate table for plotting
all_rpl_coords <- rbindlist(lapply(rpl_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
        gr$merged[, .(Gene, Position, site_class, P_aa_only)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)
all_rpl_coords[, neglog_p := -log10(P_aa_only)]

# Structure with GWAS hits - all 3 projections
p_rpl_gwas_xy <- ggplot() +
  geom_point(data = rpl_structure$ca_df, aes(x = x, y = y), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpl_coords[site_class == "not_sig"], aes(x = x, y = y), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rpl_coords[site_class == "sig_no_control"], aes(x = x, y = y), color = "gold", size = 2) +
  geom_point(data = all_rpl_coords[site_class == "sig_with_control"], aes(x = x, y = y), color = "steelblue", size = 2) +
  geom_point(data = all_rpl_coords[site_class == "sig_both"], aes(x = x, y = y), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rpl - XY (top)") + theme_classic()

p_rpl_gwas_xz <- ggplot() +
  geom_point(data = rpl_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpl_coords[site_class == "not_sig"], aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rpl_coords[site_class == "sig_no_control"], aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_rpl_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_rpl_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rpl - XZ (side)") + theme_classic()

print(p_rpl_gwas_xy + p_rpl_gwas_xz)

# Color by -log10(P)
p_rpl_pval <- ggplot() +
  geom_point(data = rpl_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpl_coords[!is.na(neglog_p)], aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() + labs(title = "rpl colored by significance") + theme_classic()

print(p_rpl_pval)

# Boxplots
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_box_surf <- ggplot(combined_rpl[!is.na(neighbor_count_min)],
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure", y = "Neighbor count (lower=exposed)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_box_int <- ggplot(combined_rpl[!is.na(dist_interface_min)],
                    aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Interface proximity", y = "Distance to interface (Å)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_surf + p_box_int)

# Ligand-specific plots if available
if ("dist_rna" %in% names(combined_rpl)) {
  p_box_rna <- ggplot(combined_rpl[!is.na(dist_rna)],
                      aes(x = site_class, y = dist_rna, fill = site_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = fill_vals) +
    labs(title = "RNA proximity", y = "Distance (Å)") +
    theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Structure with Fe-S clusters highlighted
  p_rpl_rna <- ggplot() +
    geom_point(data = rpl_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = rpl_structure$ligands[ligand_type %in% rna_nucs],
               aes(x = lig_x, y = lig_z), color = "orange", size = 4, shape = 18) +
    geom_point(data = all_rpl_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
    geom_point(data = all_rpl_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
    coord_fixed() + labs(title = "rpl: sig sites vs rna (orange)") + theme_classic()
  
  print(p_box_rna + p_rpl_rna)
}

plot(combined_rpl$neighbor_count_min, -log10(combined_rpl$P_aa_only))
plot(combined_rpl$neighbor_count_min, -log10(combined_rpl$P_aa_only))

df <- combined_rpl
df$neglog10p <- -log10(df$P_aa_only)

ct <- cor.test(df$neighbor_count_min, df$neglog10p, method = "spearman")

p_overall <-
  ggplot(df, aes(neighbor_count_min, neglog10p)) +
  geom_point(color = "gray40", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("rho = %.3f\np = %.2g", ct$estimate, ct$p.value),
    hjust = 1.1, vjust = 1.2
  ) +
  theme_classic()

p_overall

p_by_class <-
  ggplot(df, aes(neighbor_count_min, neglog10p, color = site_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

p_by_class

# ---- rp ----

pdb <- get.pdb("5X8P", format = "cif")
pdb <- read.pdb(pdb)

rp_structure <- load_pdb_structure("5X8P")
print(rp_structure$chain_summary)

message("\nLigands:")
print(rp_structure$ligands[, .N, by = ligand_type])

p_rp_xy <- ggplot(rp_structure$ca_df, aes(x = x, y = y, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XY view") + theme_classic()

p_rp_xz <- ggplot(rp_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XZ view") + theme_classic()

p_rp_yz <- ggplot(rp_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - YZ view") + theme_classic()

print(p_rp_xy + p_rp_xz + p_rp_yz)
rp_def <- DEFAULT_COMPLEXES$rp_complex

rp_result <- analyze_complex(
  "rp_complex",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS,
  genes_to_analyze = names(rp_def$chain_map)
)

rp_result$gene_results

rp_structure$ca_df$neighbor_count <- with(rp_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

hist(rp_structure$ca_df$neighbor_count)
surface_features <- rbindlist(lapply(names(rp_result$gene_results), function(gene_name) {
  gr <- rp_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  pm <- merge(gr$position_maps,
              rp_structure$ca_df[, .(chain, resno, neighbor_count)],
              by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
  
  pm[, .(neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
         neighbor_count_min = min(neighbor_count, na.rm = TRUE)), by = .(gene, aln_pos)]
}))

combined_rp <- rp_result$combined_merged
combined_rp <- merge(combined_rp, surface_features,
                      by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)

median_neighbors <- median(combined_rp$neighbor_count_min, na.rm = TRUE)
combined_rp[, surface_exposed := neighbor_count_min < median_neighbors]

# ---- LIGAND DISTANCES ----
# rp has Fe-S clusters, quinone binding sites
message("\nLigand types present:")
print(rp_structure$ligands[, .N, by = ligand_type])

# Iron-sulfur clusters (electron transfer chain)
rna_nucs <- c("U", "G", "A", "C")  # different Fe-S nomenclature
rna_coords <- rp_structure$ligands[ligand_type %in% rna_nucs, .(x = lig_x, y = lig_y, z = lig_z)]

if (nrow(rna_coords) > 0) {
  combined_rp[, dist_rna := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - rna_coords$x)^2 + (y_mean[i] - rna_coords$y)^2 + (z_mean[i] - rna_coords$z)^2))
  })]
  message("Added rna distance (", nrow(rna_coords), " clusters)")
}

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("rp HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

message("\nSite class distribution:")
print(combined_rp[, .N, by = site_class])

# H1: Surface exposure
message("\nH1: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_rp[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_rp[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# H2: Interface proximity
message("\nH2: Interface proximity")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_rp[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_rp[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: rma cluster proximity (if present)
if ("dist_rna" %in% names(combined_rp)) {
  message("\nH3: Distance to rna")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_rp[site_class == sc & !is.na(dist_rna), dist_rna]
    bg_d <- combined_rp[site_class == "not_sig" & !is.na(dist_rna), dist_rna]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# ---- VISUALIZATIONS ----

# Build coordinate table for plotting
all_rp_coords <- rbindlist(lapply(rp_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
        gr$merged[, .(Gene, Position, site_class, P_aa_only)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)
all_rp_coords[, neglog_p := -log10(P_aa_only)]

# Structure with GWAS hits - all 3 projections
p_rp_gwas_xy <- ggplot() +
  geom_point(data = rp_structure$ca_df, aes(x = x, y = y), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rp_coords[site_class == "not_sig"], aes(x = x, y = y), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rp_coords[site_class == "sig_no_control"], aes(x = x, y = y), color = "gold", size = 2) +
  geom_point(data = all_rp_coords[site_class == "sig_with_control"], aes(x = x, y = y), color = "steelblue", size = 2) +
  geom_point(data = all_rp_coords[site_class == "sig_both"], aes(x = x, y = y), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rp - XY (top)") + theme_classic()

p_rp_gwas_xz <- ggplot() +
  geom_point(data = rp_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rp_coords[site_class == "not_sig"], aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rp_coords[site_class == "sig_no_control"], aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_rp_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_rp_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rp - XZ (side)") + theme_classic()

print(p_rp_gwas_xy + p_rp_gwas_xz)

# Color by -log10(P)
p_rp_pval <- ggplot() +
  geom_point(data = rp_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rp_coords[!is.na(neglog_p)], aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() + labs(title = "rp colored by significance") + theme_classic()

print(p_rp_pval)

# Boxplots
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_box_surf <- ggplot(combined_rp[!is.na(neighbor_count_min)],
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure", y = "Neighbor count (lower=exposed)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_box_int <- ggplot(combined_rp[!is.na(dist_interface_min)],
                    aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Interface proximity", y = "Distance to interface (Å)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_surf + p_box_int)

# Ligand-specific plots if available
if ("dist_rna" %in% names(combined_rp)) {
  p_box_rna <- ggplot(combined_rp[!is.na(dist_rna)],
                      aes(x = site_class, y = dist_rna, fill = site_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = fill_vals) +
    labs(title = "RNA proximity", y = "Distance (Å)") +
    theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Structure with Fe-S clusters highlighted
  p_rp_rna <- ggplot() +
    geom_point(data = rp_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = rp_structure$ligands[ligand_type %in% rna_nucs],
               aes(x = lig_x, y = lig_z), color = "orange", size = 4, shape = 18) +
    geom_point(data = all_rp_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
    geom_point(data = all_rp_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
    coord_fixed() + labs(title = "rp: sig sites vs rna (orange)") + theme_classic()
  
  print(p_box_rna + p_rp_rna)
}

plot(combined_rp$neighbor_count_min, -log10(combined_rp$P_aa_only))
plot(combined_rp$neighbor_count_min, -log10(combined_rp$P_aa_only))

df <- combined_rp
df$neglog10p <- -log10(df$P_aa_only)

ct <- cor.test(df$neighbor_count_min, df$neglog10p, method = "spearman")

p_overall <-
  ggplot(df, aes(neighbor_count_min, neglog10p)) +
  geom_point(color = "gray40", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    title = "ribosomal proteins: aa_only vs neighbor count",
    x = "Minimum neighbor count",
    y = expression(-log[10](p))
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("rho = %.3f\np = %.2g", ct$estimate, ct$p.value),
    hjust = 1.1, vjust = 1.2
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_overall

p_by_class <-
  ggplot(df, aes(neighbor_count_min, neglog10p, color = site_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

p_by_class


# RIBOSOME: ADJUSTED SURFACE ACCESSIBILITY INCLUDING rRNA

# Get all RNA atom coordinates (not just ligand centers - need full backbone)
# RNA atoms are stored differently - check what we have
message("Checking RNA in structure...")

# RNA residues in PDB are typically: A, G, C, U (or DA, DG, DC, DT for DNA)
rna_resids <- c("A", "G", "C", "U", "RA", "RG", "RC", "RU")

# Get RNA atoms from the full PDB (not just ligands table which has centers)
rna_atoms <- rp_structure$pdb$atom[rp_structure$pdb$atom$resid %in% rna_resids, ]
message("RNA atoms found: ", nrow(rna_atoms))

if (nrow(rna_atoms) > 0) {
  rna_coords <- data.table(x = rna_atoms$x, y = rna_atoms$y, z = rna_atoms$z)
  
  # Compute neighbor count INCLUDING RNA
  # Use P atoms (phosphate backbone) as representative - less dense than all atoms
  rna_p_atoms <- rna_atoms[rna_atoms$elety == "P", ]
  message("RNA phosphate atoms: ", nrow(rna_p_atoms))
  
  if (nrow(rna_p_atoms) > 0) {
    rna_p_coords <- data.table(x = rna_p_atoms$x, y = rna_p_atoms$y, z = rna_p_atoms$z)
    
    # Combined neighbor count: protein CA + RNA P atoms
    compute_neighbor_count_with_rna <- function(prot_x, prot_y, prot_z, 
                                                rna_x, rna_y, rna_z, 
                                                radius = 10) {
      n <- length(prot_x)
      sapply(1:n, function(i) {
        # Protein neighbors
        prot_dists <- sqrt((prot_x - prot_x[i])^2 + (prot_y - prot_y[i])^2 + (prot_z - prot_z[i])^2)
        n_prot <- sum(prot_dists > 0 & prot_dists <= radius)
        
        # RNA neighbors
        rna_dists <- sqrt((rna_x - prot_x[i])^2 + (rna_y - prot_y[i])^2 + (rna_z - prot_z[i])^2)
        n_rna <- sum(rna_dists <= radius)
        
        n_prot + n_rna
      })
    }
    
    # Compute adjusted neighbor count
    rp_structure$ca_df$neighbor_count_with_rna <- compute_neighbor_count_with_rna(
      rp_structure$ca_df$x, rp_structure$ca_df$y, rp_structure$ca_df$z,
      rna_p_coords$x, rna_p_coords$y, rna_p_coords$z
    )
    
    message("\nNeighbor count comparison:")
    message("  Protein-only: ", round(median(rp_structure$ca_df$neighbor_count), 1), " (median)")
    message("  With RNA:     ", round(median(rp_structure$ca_df$neighbor_count_with_rna), 1), " (median)")
    
    # Update surface features for combined_rp
    surface_features_rna <- rbindlist(lapply(names(rp_result$gene_results), function(gene_name) {
      gr <- rp_result$gene_results[[gene_name]]
      if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
      
      pm <- merge(gr$position_maps,
                  rp_structure$ca_df[, .(chain, resno, neighbor_count_with_rna)],
                  by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
      
      pm[, .(neighbor_rna_mean = mean(neighbor_count_with_rna, na.rm = TRUE),
             neighbor_rna_min = min(neighbor_count_with_rna, na.rm = TRUE)), by = .(gene, aln_pos)]
    }))
    
    combined_rp <- merge(combined_rp, surface_features_rna,
                         by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)
    
    # ---- RETEST WITH RNA-ADJUSTED COUNTS ----
    message("\n", paste(rep("=", 60), collapse = ""))
    message("SURFACE EXPOSURE - RNA-ADJUSTED")
    message(paste(rep("=", 60), collapse = ""))
    
    message("\nH1: Surface exposure (protein-only neighbor count)")
    for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
      test_n <- combined_rp[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
      bg_n <- combined_rp[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
      if (length(test_n) >= 3) {
        wt <- wilcox.test(test_n, bg_n)
        dir <- ifelse(median(test_n) < median(bg_n), "↓ exposed", "↑ buried")
        message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                        sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
      }
    }
    
    message("\nH1b: Surface exposure (RNA-adjusted neighbor count)")
    for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
      test_n <- combined_rp[site_class == sc & !is.na(neighbor_rna_min), neighbor_rna_min]
      bg_n <- combined_rp[site_class == "not_sig" & !is.na(neighbor_rna_min), neighbor_rna_min]
      if (length(test_n) >= 3) {
        wt <- wilcox.test(test_n, bg_n)
        dir <- ifelse(median(test_n) < median(bg_n), "↓ exposed", "↑ buried")
        message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                        sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
      }
    }
    
    # ---- CORRELATION COMPARISON ----
    ct_prot <- cor.test(combined_rp$neighbor_count_min, -log10(combined_rp$P_aa_only), 
                        method = "spearman", use = "complete")
    ct_rna <- cor.test(combined_rp$neighbor_rna_min, -log10(combined_rp$P_aa_only), 
                       method = "spearman", use = "complete")
    
    message("\nCorrelation with -log10(P):")
    message(sprintf("  Protein-only: rho = %.3f, p = %.2e", ct_prot$estimate, ct_prot$p.value))
    message(sprintf("  RNA-adjusted: rho = %.3f, p = %.2e", ct_rna$estimate, ct_rna$p.value))
    
    # ---- PLOTS ----
    # Compare the two metrics
    p_compare <- ggplot(combined_rp[!is.na(neighbor_rna_min)], 
                        aes(x = neighbor_count_min, y = neighbor_rna_min, color = site_class)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("not_sig" = "gray50", "sig_no_control" = "gold",
                                    "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
      labs(title = "Protein-only vs RNA-adjusted neighbor count",
           x = "Protein-only neighbors", y = "With RNA neighbors") +
      theme_classic()
    
    p_rna_scatter <- ggplot(combined_rp[!is.na(neighbor_rna_min)],
                            aes(x = neighbor_rna_min, y = -log10(P_aa_only), color = site_class)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("not_sig" = "gray50", "sig_no_control" = "gold",
                                    "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
      geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
               label = sprintf("rho = %.3f\np = %.2g", ct_rna$estimate, ct_rna$p.value)) +
      labs(title = "RNA-adjusted surface vs significance",
           x = "Neighbor count (protein + RNA)", y = "-log10(P)") +
      theme_classic()
    
    print(p_compare + p_rna_scatter)
    
    # Boxplot comparison
    p_box_prot <- ggplot(combined_rp[!is.na(neighbor_count_min)],
                         aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      scale_fill_manual(values = fill_vals) +
      labs(title = "Protein-only", y = "Neighbor count") +
      theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    p_box_rna <- ggplot(combined_rp[!is.na(neighbor_rna_min)],
                        aes(x = site_class, y = neighbor_rna_min, fill = site_class)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      scale_fill_manual(values = fill_vals) +
      labs(title = "RNA-adjusted", y = "Neighbor count") +
      theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_box_prot + p_box_rna)
    
    # Structure visualization - color by RNA-adjusted neighbors
    p_rp_surface_rna <- ggplot() +
      geom_point(data = rp_structure$ca_df, aes(x = x, y = z, color = neighbor_count_with_rna),
                 size = 0.8, alpha = 0.7) +
      scale_color_viridis_c(option = "plasma", direction = -1, name = "Neighbors\n(+RNA)") +
      geom_point(data = all_rp_coords[site_class == "sig_both"], aes(x = x, y = z),
                 color = "red", size = 3, shape = 1, stroke = 1.5) +
      coord_fixed() +
      labs(title = "Ribosome: RNA-adjusted surface (sig_both = red circles)") +
      theme_classic()
    
    print(p_rp_surface_rna)
    
  } else {
    message("No RNA P atoms found - check atom naming")
  }
} else {
  message("No RNA found in structure")
}

# ---- rpo ---- 

rpo_structure <- load_pdb_structure("8WA1")
print(rpo_structure$chain_summary)

message("\nLigands:")
print(rpo_structure$ligands[, .N, by = ligand_type])

# Quick structure visualization
p_rpo_xy <- ggplot(rpo_structure$ca_df, aes(x = x, y = y, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XY view") + theme_classic()

p_rpo_xz <- ggplot(rpo_structure$ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - XZ view") + theme_classic()

p_rpo_yz <- ggplot(rpo_structure$ca_df, aes(x = y, y = z, color = chain)) +
  geom_point(size = 0.5, alpha = 0.7) + coord_fixed() +
  labs(title = "NDH - YZ view") + theme_classic()

print(p_rpo_xy + p_rpo_xz + p_rpo_yz)
rpo_def <- DEFAULT_COMPLEXES$rpo_complex

rpo_result <- analyze_complex(
  "rpo_complex",
  ALIGNMENT_DIR,
  REFERENCE_PATTERN,
  sites_df,
  THRESH_AA_ONLY,
  THRESH_AA_PCS,
  genes_to_analyze = names(rpo_def$chain_map)
)

rpo_result$gene_results

rpo_structure$ca_df$neighbor_count <- with(rpo_structure$ca_df, 
                                           compute_neighbor_count(x, y, z))

hist(rpo_structure$ca_df$neighbor_count)
surface_features <- rbindlist(lapply(names(rpo_result$gene_results), function(gene_name) {
  gr <- rpo_result$gene_results[[gene_name]]
  if (is.null(gr) || is.null(gr$position_maps)) return(NULL)
  
  pm <- merge(gr$position_maps,
              rpo_structure$ca_df[, .(chain, resno, neighbor_count)],
              by.x = c("chain", "pdb_resno"), by.y = c("chain", "resno"), all.x = TRUE)
  
  pm[, .(neighbor_count_mean = mean(neighbor_count, na.rm = TRUE),
         neighbor_count_min = min(neighbor_count, na.rm = TRUE)), by = .(gene, aln_pos)]
}))

combined_rpo <- rpo_result$combined_merged
combined_rpo <- merge(combined_rpo, surface_features,
                      by.x = c("Gene", "Position"), by.y = c("gene", "aln_pos"), all.x = TRUE)

median_neighbors <- median(combined_rpo$neighbor_count_min, na.rm = TRUE)
combined_rpo[, surface_exposed := neighbor_count_min < median_neighbors]

# ---- LIGAND DISTANCES ----
# rpo has Fe-S clusters, quinone binding sites
message("\nLigand types present:")
print(rpo_structure$ligands[, .N, by = ligand_type])

# Iron-sulfur clusters (electron transfer chain)
rna_nucs <- c("U", "G", "A", "C")  # different Fe-S nomenclature
rna_coords <- rpo_structure$ligands[ligand_type %in% rna_nucs, .(x = lig_x, y = lig_y, z = lig_z)]

if (nrow(rna_coords) > 0) {
  combined_rpo[, dist_rna := sapply(1:.N, function(i) {
    min(sqrt((x_mean[i] - rna_coords$x)^2 + (y_mean[i] - rna_coords$y)^2 + (z_mean[i] - rna_coords$z)^2))
  })]
  message("Added rna distance (", nrow(rna_coords), " clusters)")
}

# ---- HYPOTHESIS TESTS ----
message("\n", paste(rep("=", 60), collapse = ""))
message("rpo HYPOTHESIS TESTS")
message(paste(rep("=", 60), collapse = ""))

message("\nSite class distribution:")
print(combined_rpo[, .N, by = site_class])

# H1: Surface exposure
message("\nH1: Surface exposure (neighbor count)")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_rpo[site_class == sc & !is.na(neighbor_count_min), neighbor_count_min]
  bg_n <- combined_rpo[site_class == "not_sig" & !is.na(neighbor_count_min), neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# H2: Interface proximity
message("\nH2: Interface proximity")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_d <- combined_rpo[site_class == sc & !is.na(dist_interface_min), dist_interface_min]
  bg_d <- combined_rpo[site_class == "not_sig" & !is.na(dist_interface_min), dist_interface_min]
  if (length(test_d) >= 3) {
    wt <- wilcox.test(test_d, bg_d)
    dir <- ifelse(median(test_d) < median(bg_d), "↓ at interface", "↑ away")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
  }
}

# H3: rma cluster proximity (if present)
if ("dist_rna" %in% names(combined_rpo)) {
  message("\nH3: Distance to rna")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_d <- combined_rpo[site_class == sc & !is.na(dist_rna), dist_rna]
    bg_d <- combined_rpo[site_class == "not_sig" & !is.na(dist_rna), dist_rna]
    if (length(test_d) >= 3) {
      wt <- wilcox.test(test_d, bg_d)
      dir <- ifelse(median(test_d) > median(bg_d), "↑ farther", "↓ closer")
      message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                      sc, length(test_d), median(test_d), dir, median(bg_d), wt$p.value))
    }
  }
}

# ---- VISUALIZATIONS ----

# Build coordinate table for plotting
all_rpo_coords <- rbindlist(lapply(rpo_result$gene_results, function(gr) {
  if (is.null(gr$position_maps)) return(NULL)
  merge(gr$position_maps[, .(gene, aln_pos, chain, x, y, z)],
        gr$merged[, .(Gene, Position, site_class, P_aa_only)],
        by.x = c("gene", "aln_pos"), by.y = c("Gene", "Position"), all.x = TRUE)
}), fill = TRUE)
all_rpo_coords[, neglog_p := -log10(P_aa_only)]

# Structure with GWAS hits - all 3 projections
p_rpo_gwas_xy <- ggplot() +
  geom_point(data = rpo_structure$ca_df, aes(x = x, y = y), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpo_coords[site_class == "not_sig"], aes(x = x, y = y), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rpo_coords[site_class == "sig_no_control"], aes(x = x, y = y), color = "gold", size = 2) +
  geom_point(data = all_rpo_coords[site_class == "sig_with_control"], aes(x = x, y = y), color = "steelblue", size = 2) +
  geom_point(data = all_rpo_coords[site_class == "sig_both"], aes(x = x, y = y), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rpo - XY (top)") + theme_classic()

p_rpo_gwas_xz <- ggplot() +
  geom_point(data = rpo_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpo_coords[site_class == "not_sig"], aes(x = x, y = z), color = "gray50", size = 1, alpha = 0.5) +
  geom_point(data = all_rpo_coords[site_class == "sig_no_control"], aes(x = x, y = z), color = "gold", size = 2) +
  geom_point(data = all_rpo_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
  geom_point(data = all_rpo_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
  coord_fixed() + labs(title = "rpo - XZ (side)") + theme_classic()

print(p_rpo_gwas_xy + p_rpo_gwas_xz)

# Color by -log10(P)
p_rpo_pval <- ggplot() +
  geom_point(data = rpo_structure$ca_df, aes(x = x, y = z), color = "gray60", size = 0.3, alpha = 0.3) +
  geom_point(data = all_rpo_coords[!is.na(neglog_p)], aes(x = x, y = z, color = neglog_p), size = 1.5) +
  scale_color_viridis_c(option = "inferno", name = "-log10(P)") +
  coord_fixed() + labs(title = "rpo colored by significance") + theme_classic()

print(p_rpo_pval)

# Boxplots
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_box_surf <- ggplot(combined_rpo[!is.na(neighbor_count_min)],
                     aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure", y = "Neighbor count (lower=exposed)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_box_int <- ggplot(combined_rpo[!is.na(dist_interface_min)],
                    aes(x = site_class, y = dist_interface_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Interface proximity", y = "Distance to interface (Å)") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p_box_surf + p_box_int)

# Ligand-specific plots if available
if ("dist_rna" %in% names(combined_rpo)) {
  p_box_rna <- ggplot(combined_rpo[!is.na(dist_rna)],
                      aes(x = site_class, y = dist_rna, fill = site_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = fill_vals) +
    labs(title = "RNA proximity", y = "Distance (Å)") +
    theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Structure with Fe-S clusters highlighted
  p_rpo_rna <- ggplot() +
    geom_point(data = rpo_structure$ca_df, aes(x = x, y = z), color = "gray80", size = 0.3, alpha = 0.3) +
    geom_point(data = rpo_structure$ligands[ligand_type %in% rna_nucs],
               aes(x = lig_x, y = lig_z), color = "orange", size = 4, shape = 18) +
    geom_point(data = all_rpo_coords[site_class == "sig_both"], aes(x = x, y = z), color = "darkred", size = 2.5) +
    geom_point(data = all_rpo_coords[site_class == "sig_with_control"], aes(x = x, y = z), color = "steelblue", size = 2) +
    coord_fixed() + labs(title = "rpo: sig sites vs rna (orange)") + theme_classic()
  
  print(p_box_rna + p_rpo_rna)
}

plot(combined_rpo$neighbor_count_min, -log10(combined_rpo$P_aa_only))
plot(combined_rpo$neighbor_count_min, -log10(combined_rpo$P_aa_only))

df <- combined_rpo
df$neglog10p <- -log10(df$P_aa_only)

ct <- cor.test(df$neighbor_count_min, df$neglog10p, method = "spearman")

p_overall <-
  ggplot(df, aes(neighbor_count_min, neglog10p)) +
  geom_point(color = "gray40", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf("rho = %.3f\np = %.2g", ct$estimate, ct$p.value),
    hjust = 1.1, vjust = 1.2
  ) +
  theme_classic()

p_overall

p_by_class <-
  ggplot(df, aes(neighbor_count_min, neglog10p, color = site_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

p_by_class

# ---- surface area ----


# COMBINED SURFACE ACCESSIBILITY ANALYSIS ACROSS ALL COMPLEXES
# Combine all complex data with complex labels
combined_all <- rbindlist(list(
  combined_atp[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class, neighbor_count_min, complex = "ATP synthase")],
  combined_b6f[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class,neighbor_count_min, complex = "Cyt b6f")],
  combined_psII[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class,neighbor_count_min, complex = "PSII")],
  combined_psI[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class,neighbor_count_min, complex = "PSI")]
  #combined_ndh[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class, neighbor_count_min, complex = "NDH")],
  #rbcL_merged[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class,neighbor_count_min, complex = "Rubisco")],
  #combined_rp[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class, neighbor_count_min, complex = "Ribosome")],
  #combined_rpo[, .(Gene, Position, P_aa_only, P_aa_with_pcs, site_class,neighbor_count_min, complex = "RNA Pol")]
), fill = TRUE)

combined_all[, neglog_p := -log10(P_aa_only)]
combined_all[, neglog_p_pcs := -log10(P_aa_with_pcs)]

# Filter to complete cases
combined_all <- combined_all[!is.na(neighbor_count_min) & !is.na(neglog_p)]

stopifnot(nrow(combined_all) > 0)

# ---- OVERALL CORRELATION ----
message("\n", paste(rep("=", 60), collapse = ""))
message("COMBINED SURFACE ACCESSIBILITY ANALYSIS")
message(paste(rep("=", 60), collapse = ""))

message("\nTotal residues: ", nrow(combined_all))
print(combined_all[, .N, by = complex])

ct_all <- cor.test(combined_all$neighbor_count_min, combined_all$neglog_p, 
                   method = "spearman")
message(sprintf("\nOVERALL correlation (neighbor count vs -log10 P_aa_only):"))
message(sprintf("  rho = %.4f, p = %.2e", ct_all$estimate, ct_all$p.value))

# ---- PER-COMPLEX CORRELATIONS ----
message("\n--- Per-complex correlations ---")
cor_by_complex <- combined_all[, {
  ct <- cor.test(neighbor_count_min, neglog_p, method = "spearman")
  .(rho = ct$estimate, p = ct$p.value, n = .N)
}, by = complex]

print(cor_by_complex[order(p)])

# ---- WILCOXON TESTS BY SITE CLASS ----
message("\n--- Surface exposure by site class (all complexes) ---")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_n <- combined_all[site_class == sc, neighbor_count_min]
  bg_n <- combined_all[site_class == "not_sig", neighbor_count_min]
  if (length(test_n) >= 3) {
    wt <- wilcox.test(test_n, bg_n)
    dir <- ifelse(median(test_n) < median(bg_n), "↓ MORE exposed", "↑ MORE buried")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.2e",
                    sc, length(test_n), median(test_n), dir, median(bg_n), wt$p.value))
  }
}

# ---- PLOTS ----

# 1. Overall scatter with regression
p_overall <- ggplot(combined_all, aes(neighbor_count_min, neglog_p)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
           label = sprintf("rho = %.3f\np = %.2g", ct_all$estimate, ct_all$p.value)) +
  labs(title = "Surface exposure vs significance (all complexes)",
       x = "Neighbor count (lower = more exposed)",
       y = expression(-log[10](P[aa_only]))) +
  theme_classic()

# 2. Faceted by complex
p_by_complex <- ggplot(combined_all, aes(neighbor_count_min, neglog_p)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  facet_wrap(~complex, scales = "free") +
  labs(title = "Surface exposure vs significance by complex",
       x = "Neighbor count", y = expression(-log[10](P))) +
  theme_classic()

# 3. Boxplot by site class
fill_vals <- c("not_sig" = "gray60", "sig_no_control" = "gold",
               "sig_with_control" = "steelblue", "sig_both" = "darkred")

p_boxplot <- ggplot(combined_all, aes(x = site_class, y = neighbor_count_min, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_fill_manual(values = fill_vals) +
  labs(title = "Surface exposure by significance class (all complexes)",
       y = "Neighbor count (lower = more exposed)", x = "") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Forest plot of per-complex correlations
cor_by_complex[, complex := factor(complex, levels = complex[order(rho)])]
p_forest <- ggplot(cor_by_complex, aes(x = rho, y = complex)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n, color = p < 0.05)) +
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "gray50"), 
                     name = "p < 0.05") +
  labs(title = "Per-complex correlation (neighbor count vs -log10 P)",
       x = "Spearman rho", y = "") +
  theme_classic()

print(p_overall)
print(p_by_complex)
print(p_boxplot)
print(p_forest)