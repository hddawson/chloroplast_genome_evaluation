library(data.table)
library(ape)
library(arrow)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD GWAS RESULTS ----
message("Loading GWAS results...")

base_out_dir <- "results/npc_experiment"

load_npc_results <- function(n_pcs) {
  dir_path <- file.path(base_out_dir, paste0("npc_", n_pcs))
  rds_files <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) return(NULL)
  
  all_res <- rbindlist(lapply(rds_files, function(f) {
    res_list <- readRDS(f)
    rbindlist(lapply(res_list, function(r) {
      # Extract effect sizes if available
      max_effect <- NA_real_
      if (!is.null(r$effects) && nrow(r$effects) > 0) {
        max_effect <- max(abs(r$effects$Effect), na.rm = TRUE)
      }
      data.table(
        Gene = r$Gene,
        Position = r$Aligned_Position,
        n_pcs = r$n_pcs,
        N = r$N,
        R2_aa = r$R2_aa,
        R2_full = r$R2_full,
        P_aa_only = r$P_aa_only,
        P_aa_with_pcs = r$P_aa_with_pcs,
        max_effect = max_effect
      )
    }))
  }), fill = TRUE)
  return(all_res)
}

gwas_0 <- load_npc_results(0)
gwas_1000 <- load_npc_results(1000)

stopifnot(nrow(gwas_0) > 0, nrow(gwas_1000) > 0)

# Create site IDs
gwas_0[, site := paste(Gene, Position, sep = "_")]
gwas_1000[, site := paste(Gene, Position, sep = "_")]

message("Loaded ", nrow(gwas_0), " sites from n_pcs=0")
message("Loaded ", nrow(gwas_1000), " sites from n_pcs=1000")

# ---- 2. CLASSIFY SITES ----
message("\nClassifying sites...")

# Top 20th percentile for n_pcs=0 (P_aa_only)
thresh_0 <- quantile(gwas_0$P_aa_only, 0.20, na.rm = TRUE)
sig_0 <- gwas_0[P_aa_only <= thresh_0, .(site, Gene, Position, P_0 = P_aa_only, effect_0 = max_effect)]

# Top 5th percentile for n_pcs=1000 (P_aa_with_pcs)
thresh_1000 <- quantile(gwas_1000$P_aa_with_pcs, 0.05, na.rm = TRUE)
sig_1000 <- gwas_1000[P_aa_with_pcs <= thresh_1000, .(site, Gene, Position, P_1000 = P_aa_with_pcs, effect_1000 = max_effect)]

message("Sites in top 20% (n_pcs=0): ", nrow(sig_0))
message("Sites in top 5% (n_pcs=1000): ", nrow(sig_1000))

# Classify
sites_0 <- sig_0$site
sites_1000 <- sig_1000$site

sig_both <- intersect(sites_0, sites_1000)
sig_no_control <- setdiff(sites_0, sites_1000)  # sig without PC control, not with
sig_control <- setdiff(sites_1000, sites_0)     # sig with PC control, not without

message("\nSite classification:")
message("  sig_both (top 20% @ 0pc AND top 5% @ 1000pc): ", length(sig_both))
message("  sig_no_control (top 20% @ 0pc only): ", length(sig_no_control))
message("  sig_control (top 5% @ 1000pc only): ", length(sig_control))

# Create lookup tables with gene/position info
site_info <- rbind(
  sig_0[, .(site, Gene, Position)],
  sig_1000[, .(site, Gene, Position)]
) |> unique()

site_class <- data.table(
  site = c(sig_both, sig_no_control, sig_control),
  class = c(rep("sig_both", length(sig_both)),
            rep("sig_no_control", length(sig_no_control)),
            rep("sig_control", length(sig_control)))
)
site_class <- merge(site_class, site_info, by = "site")

# Add effect sizes
site_class <- merge(site_class, sig_0[, .(site, effect_0)], by = "site", all.x = TRUE)
site_class <- merge(site_class, sig_1000[, .(site, effect_1000)], by = "site", all.x = TRUE)

message("\nSites by gene (sig_both):")
print(site_class[class == "sig_both", .N, by = Gene][order(-N)][1:10])

# ---- 3. LOAD HOLDOUT CLADES AND PREDICTION RESULTS ----
message("\nLoading holdout clades and prediction results...")

# Load tree and data
tree <- read.tree("results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree")
data <- as.data.table(read_parquet("data/processed_data.parquet"))

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
pheno <- setNames(data[[pheno_col]], data$ID)
pheno <- pheno[tree$tip.label]
pheno <- pheno[!is.na(pheno)]
tree <- keep.tip(tree, names(pheno))

# Recreate holdout clades with same seed
set.seed(42)

select_holdout_clades <- function(tree, target_frac = 0.15, min_clade = 10, max_clade = 30) {
  n_tips <- Ntip(tree)
  target_n <- floor(n_tips * target_frac)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  candidates <- list()
  for (node in internal_nodes) {
    tips <- extract.clade(tree, node)$tip.label
    if (length(tips) >= min_clade && length(tips) <= max_clade) {
      candidates[[as.character(node)]] <- tips
    }
  }
  selected <- list()
  used_tips <- character(0)
  candidate_order <- sample(names(candidates))
  for (node in candidate_order) {
    tips <- candidates[[node]]
    if (!any(tips %in% used_tips)) {
      selected[[node]] <- tips
      used_tips <- c(used_tips, tips)
      if (length(used_tips) >= target_n) break
    }
  }
  return(selected)
}

holdout_clades <- select_holdout_clades(tree)
message("Loaded ", length(holdout_clades), " holdout clades")

# Load prediction results
cv_results <- readRDS("results/phylo_cv_clade_results.rds")
clade_perf <- cv_results$clade_results
stopifnot("spearman_standard" %in% names(clade_perf))

message("Clade performance range: ", 
        round(min(clade_perf$spearman_standard, na.rm = TRUE), 3), " to ",
        round(max(clade_perf$spearman_standard, na.rm = TRUE), 3))

# ---- 4. CHECK SEGREGATION AT GWAS SITES FOR EACH CLADE ----
message("\nChecking segregation at GWAS sites...")

# Get unique genes with significant sites
sig_genes <- unique(site_class$Gene)
message("Genes with significant sites: ", length(sig_genes))

# Load alignments for these genes
aln_list <- list()
for (gene in sig_genes) {
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  if (file.exists(aln_file)) {
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (!is.null(aln)) {
      names(aln) <- sub("\\|.*", "", names(aln))
      aln_list[[gene]] <- as.matrix(aln)
    }
  }
}
message("Loaded alignments for ", length(aln_list), " genes")

# Function to check if a site segregates in a set of species
check_segregation <- function(aln_mat, position, species_ids) {
  # Get species present in alignment
  available <- intersect(species_ids, rownames(aln_mat))
  if (length(available) < 2) return(list(segregating = NA, n_alleles = NA, n_species = length(available)))
  
  # Get residues at this position
  residues <- aln_mat[available, position]
  residues <- residues[residues != "-"]  # exclude gaps
  
  if (length(residues) < 2) return(list(segregating = NA, n_alleles = NA, n_species = length(residues)))
  
  n_alleles <- length(unique(residues))
  segregating <- n_alleles >= 2
  
  return(list(segregating = segregating, n_alleles = n_alleles, n_species = length(residues)))
}

# For each clade, count segregating sites in each class
clade_segregation <- data.table(
  clade = character(),
  n_species = integer(),
  # sig_both
  n_sites_both = integer(),
  n_seg_both = integer(),
  prop_seg_both = numeric(),
  # sig_no_control
  n_sites_no_control = integer(),
  n_seg_no_control = integer(),
  prop_seg_no_control = numeric(),
  # sig_control
  n_sites_control = integer(),
  n_seg_control = integer(),
  prop_seg_control = numeric()
)

for (clade_name in names(holdout_clades)) {
  species_ids <- holdout_clades[[clade_name]]
  
  # Count segregation for each site class
  count_seg <- function(sites_in_class) {
    if (length(sites_in_class) == 0) return(list(n_sites = 0, n_seg = 0))
    
    class_sites <- site_class[site %in% sites_in_class]
    n_seg <- 0
    n_checked <- 0
    
    for (i in 1:nrow(class_sites)) {
      gene <- class_sites$Gene[i]
      pos <- class_sites$Position[i]
      
      if (!gene %in% names(aln_list)) next
      aln_mat <- aln_list[[gene]]
      if (pos > ncol(aln_mat)) next
      
      result <- check_segregation(aln_mat, pos, species_ids)
      if (!is.na(result$segregating)) {
        n_checked <- n_checked + 1
        if (result$segregating) n_seg <- n_seg + 1
      }
    }
    return(list(n_sites = n_checked, n_seg = n_seg))
  }
  
  seg_both <- count_seg(sig_both)
  seg_no_control <- count_seg(sig_no_control)
  seg_control <- count_seg(sig_control)
  
  clade_segregation <- rbind(clade_segregation, data.table(
    clade = clade_name,
    n_species = length(species_ids),
    n_sites_both = seg_both$n_sites,
    n_seg_both = seg_both$n_seg,
    prop_seg_both = if (seg_both$n_sites > 0) seg_both$n_seg / seg_both$n_sites else NA_real_,
    n_sites_no_control = seg_no_control$n_sites,
    n_seg_no_control = seg_no_control$n_seg,
    prop_seg_no_control = if (seg_no_control$n_sites > 0) seg_no_control$n_seg / seg_no_control$n_sites else NA_real_,
    n_sites_control = seg_control$n_sites,
    n_seg_control = seg_control$n_seg,
    prop_seg_control = if (seg_control$n_sites > 0) seg_control$n_seg / seg_control$n_sites else NA_real_
  ))
  
  if (which(names(holdout_clades) == clade_name) %% 10 == 0) {
    message("  Processed ", which(names(holdout_clades) == clade_name), "/", length(holdout_clades), " clades")
  }
}

message("Segregation analysis complete")

# ---- 5. MERGE WITH PREDICTION PERFORMANCE ----
clade_analysis <- merge(clade_segregation, clade_perf[, .(clade, spearman_standard, spearman_strict)], by = "clade")

message("\n=== CORRELATION: Segregation vs Prediction Performance ===")

# Test correlations
cor_both <- cor.test(clade_analysis$prop_seg_both, clade_analysis$spearman_standard, 
                     method = "spearman", use = "complete.obs")
cor_no_control <- cor.test(clade_analysis$prop_seg_no_control, clade_analysis$spearman_standard,
                           method = "spearman", use = "complete.obs")
cor_control <- cor.test(clade_analysis$prop_seg_control, clade_analysis$spearman_standard,
                        method = "spearman", use = "complete.obs")

message("\nProportion segregating (sig_both) vs within-clade Spearman:")
message("  rho = ", round(cor_both$estimate, 3), ", p = ", format.pval(cor_both$p.value, digits = 3))

message("\nProportion segregating (sig_no_control) vs within-clade Spearman:")
message("  rho = ", round(cor_no_control$estimate, 3), ", p = ", format.pval(cor_no_control$p.value, digits = 3))

message("\nProportion segregating (sig_control) vs within-clade Spearman:")
message("  rho = ", round(cor_control$estimate, 3), ", p = ", format.pval(cor_control$p.value, digits = 3))

# Also test raw counts
cor_count <- cor.test(clade_analysis$n_seg_both, clade_analysis$spearman_standard,
                      method = "spearman", use = "complete.obs")
message("\nN segregating sites (sig_both) vs within-clade Spearman:")
message("  rho = ", round(cor_count$estimate, 3), ", p = ", format.pval(cor_count$p.value, digits = 3))

# ---- 6. PLOTS ----
message("\nGenerating plots...")

# Main hypothesis plot
p1 <- ggplot(clade_analysis, aes(x = prop_seg_both, y = spearman_standard)) +
  geom_point(aes(size = n_species), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Within-Clade Prediction vs Segregation at GWAS Sites",
    subtitle = paste0("sig_both sites | rho = ", round(cor_both$estimate, 3), 
                      ", p = ", format.pval(cor_both$p.value, digits = 2)),
    x = "Proportion of sig_both sites segregating in clade",
    y = "Within-clade Spearman r",
    size = "Clade size"
  ) +
  theme_minimal()

# Compare all three site classes
plot_data <- melt(clade_analysis, 
                  id.vars = c("clade", "spearman_standard", "n_species"),
                  measure.vars = c("prop_seg_both", "prop_seg_no_control", "prop_seg_control"),
                  variable.name = "site_class", value.name = "prop_segregating")
plot_data[, site_class := factor(site_class, 
                                  levels = c("prop_seg_both", "prop_seg_no_control", "prop_seg_control"),
                                  labels = c("sig_both", "sig_no_control", "sig_control"))]

p2 <- ggplot(plot_data, aes(x = prop_segregating, y = spearman_standard)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~site_class, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Prediction Performance vs Segregation by Site Class",
    x = "Proportion of sites segregating",
    y = "Within-clade Spearman r"
  ) +
  theme_minimal()

# Distribution of segregation proportions
p3 <- ggplot(plot_data, aes(x = prop_segregating, fill = site_class)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  facet_wrap(~site_class, ncol = 1) +
  labs(
    title = "Distribution of Segregation Proportions Across Clades",
    x = "Proportion of sites segregating",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Scatter of segregation counts
p4 <- ggplot(clade_analysis, aes(x = n_seg_both, y = spearman_standard)) +
  geom_point(aes(color = prop_seg_both), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_viridis_c(name = "Prop.\nsegregating") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Prediction vs Number of Segregating sig_both Sites",
    x = "Number of segregating sites",
    y = "Within-clade Spearman r"
  ) +
  theme_minimal()

ggsave("results/segregation_vs_prediction_main.png", p1, width = 8, height = 6)
ggsave("results/segregation_vs_prediction_byclass.png", p2, width = 12, height = 5)
ggsave("results/segregation_distribution.png", p3, width = 8, height = 8)
ggsave("results/segregation_count_vs_prediction.png", p4, width = 8, height = 6)

message("Plots saved to results/")

# ---- 7. SUMMARY TABLE ----
message("\n=== CLADE SUMMARY (sorted by spearman) ===")
print(clade_analysis[order(-spearman_standard), 
                     .(clade, n_species, spearman_standard, 
                       n_seg_both, prop_seg_both,
                       n_seg_control, prop_seg_control)])

# ---- 8. SAVE RESULTS ----
results <- list(
  site_class = site_class,
  clade_segregation = clade_segregation,
  clade_analysis = clade_analysis,
  correlations = list(
    sig_both = cor_both,
    sig_no_control = cor_no_control,
    sig_control = cor_control,
    count_both = cor_count
  ),
  thresholds = list(
    thresh_0_pct20 = thresh_0,
    thresh_1000_pct5 = thresh_1000
  )
)

saveRDS(results, "results/segregation_analysis_results.rds")
message("\nResults saved to results/segregation_analysis_results.rds")