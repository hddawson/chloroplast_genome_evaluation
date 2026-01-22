# ==== PHYLOGENETIC AUTOCORRELATION OF AA VARIANTS BY PHENOTYPE ====
library(ape)
library(data.table)
library(ggplot2)

# ==== CONFIGURATION ====
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

# ==== LOAD DATA ====
cat("Loading data...\n")
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)

data <- as.data.table(arrow::read_parquet(data_file))
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID

tree <- read.tree(tree_file)
tree$node.label <- NULL

partition_map <- readRDS(partition_map_file)
setDT(partition_map)

# ==== LOAD GWAS EFFECTS ====
model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

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

# ==== PREPARE COMMON TAXA ====
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

# Subset alignment and phenotype
aln_sub <- aln_mat[common_taxa, ]
pheno_sub <- pheno[common_taxa]

# ==== MORAN'S I FOR PHYLOGENETIC AUTOCORRELATION ====
cat("\nCalculating phylogenetic weights matrix...\n")

# Inverse phylogenetic distance matrix as weights
phy_dist <- cophenetic.phylo(tree_pruned)
phy_dist <- phy_dist[common_taxa, common_taxa]

# Inverse distance weights (with small constant to avoid Inf)
W <- 1 / (phy_dist + 0.001)
diag(W) <- 0
# Row-standardize
W <- W / rowSums(W)

# Moran's I function
calc_moran_I <- function(x, W) {
  n <- length(x)
  x_centered <- x - mean(x, na.rm = TRUE)
  
  numerator <- sum(W * outer(x_centered, x_centered), na.rm = TRUE)
  denominator <- sum(x_centered^2, na.rm = TRUE)
  
  S0 <- sum(W, na.rm = TRUE)
  
  I <- (n / S0) * (numerator / denominator)
  
  # Expected value under null
  E_I <- -1 / (n - 1)
  
  # Variance (simplified)
  S1 <- 0.5 * sum((W + t(W))^2)
  S2 <- sum((rowSums(W) + colSums(W))^2)
  k <- (sum(x_centered^4) / n) / (sum(x_centered^2) / n)^2
  
  var_I <- (n * ((n^2 - 3*n + 3) * S1 - n * S2 + 3 * S0^2) -
              k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * S0^2)) /
    ((n - 1) * (n - 2) * (n - 3) * S0^2) - E_I^2
  
  z <- (I - E_I) / sqrt(var_I)
  p <- 2 * pnorm(-abs(z))
  
  list(I = I, E_I = E_I, z = z, p = p)
}

# ==== FUNCTION: ENCODE AA AS NUMERIC FOR MORAN'S I ====
# Use modal allele as reference (0), others as 1
encode_site <- function(site_col, taxa) {
  aa <- site_col[taxa]
  valid_aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  aa[!aa %in% valid_aa] <- NA
  
  if (sum(!is.na(aa)) < 50) return(rep(NA, length(aa)))
  
  # Modal allele
  tab <- table(aa)
  modal <- names(tab)[which.max(tab)]
  
  # Binary: 0 = modal, 1 = non-modal
  as.numeric(aa != modal)
}

# ==== TEST: MORAN'S I FOR RAW GENOTYPES ====
cat("\nCalculating Moran's I for raw genotypes...\n")

calc_site_moran <- function(pos, aln_sub, W, common_taxa) {
  site_encoded <- encode_site(aln_sub[, pos], common_taxa)
  if (all(is.na(site_encoded)) || var(site_encoded, na.rm = TRUE) == 0) {
    return(list(I = NA, p = NA))
  }
  calc_moran_I(site_encoded, W)
}

# Sample sites for speed (or do all if feasible)
set.seed(42)
n_sample <- 500

sig_both_pos <- gwas_effects[sig_class == "sig_both", GlobalPos]
sig_both_pos <- sig_both_pos[sig_both_pos <= ncol(aln_sub)]
not_sig_pos <- gwas_effects[sig_class == "not_sig", GlobalPos]
not_sig_pos <- not_sig_pos[not_sig_pos <= ncol(aln_sub)]

sample_sig <- sig_both_pos[sample(min(length(sig_both_pos), n_sample))]
sample_not <- not_sig_pos[sample(min(length(not_sig_pos), n_sample))]

cat("Testing", length(sample_sig), "sig_both and", length(sample_not), "not_sig sites\n")

moran_sig <- rbindlist(lapply(sample_sig, function(pos) {
  res <- calc_site_moran(pos, aln_sub, W, common_taxa)
  data.table(GlobalPos = pos, sig_class = "sig_both", I = res$I, p = res$p)
}))

moran_not <- rbindlist(lapply(sample_not, function(pos) {
  res <- calc_site_moran(pos, aln_sub, W, common_taxa)
  data.table(GlobalPos = pos, sig_class = "not_sig", I = res$I, p = res$p)
}))

moran_raw <- rbind(moran_sig, moran_not)
moran_raw <- moran_raw[!is.na(I)]

cat("\n=== Moran's I for Raw Genotypes ===\n")
print(moran_raw[, .(mean_I = mean(I), median_I = median(I), 
                    pct_sig = 100 * mean(p < 0.05)), by = sig_class])

# ==== TEST: MORAN'S I FOR RESIDUALS (GENOTYPE | PHENOTYPE) ====
cat("\nCalculating Moran's I for genotype residuals conditional on phenotype...\n")

calc_residual_moran <- function(pos, aln_sub, pheno_sub, W, common_taxa) {
  site_encoded <- encode_site(aln_sub[, pos], common_taxa)
  
  if (all(is.na(site_encoded)) || var(site_encoded, na.rm = TRUE) == 0) {
    return(list(I_raw = NA, I_resid = NA, p_raw = NA, p_resid = NA))
  }
  
  # Raw Moran's I
  raw <- calc_moran_I(site_encoded, W)
  
  # Residuals from phenotype regression
  valid_idx <- !is.na(site_encoded) & !is.na(pheno_sub)
  if (sum(valid_idx) < 50) {
    return(list(I_raw = raw$I, I_resid = NA, p_raw = raw$p, p_resid = NA))
  }
  
  fit <- lm(site_encoded ~ pheno_sub)
  resid <- residuals(fit)
  
  # Expand residuals to full vector
  resid_full <- rep(NA, length(site_encoded))
  resid_full[valid_idx] <- resid
  
  # Moran's I on residuals
  resid_moran <- calc_moran_I(resid_full, W)
  
  list(I_raw = raw$I, I_resid = resid_moran$I, 
       p_raw = raw$p, p_resid = resid_moran$p)
}

moran_resid_sig <- rbindlist(lapply(sample_sig, function(pos) {
  res <- calc_residual_moran(pos, aln_sub, pheno_sub, W, common_taxa)
  data.table(GlobalPos = pos, sig_class = "sig_both", 
             I_raw = res$I_raw, I_resid = res$I_resid,
             p_raw = res$p_raw, p_resid = res$p_resid)
}))

moran_resid_not <- rbindlist(lapply(sample_not, function(pos) {
  res <- calc_residual_moran(pos, aln_sub, pheno_sub, W, common_taxa)
  data.table(GlobalPos = pos, sig_class = "not_sig",
             I_raw = res$I_raw, I_resid = res$I_resid,
             p_raw = res$p_raw, p_resid = res$p_resid)
}))

moran_resid <- rbind(moran_resid_sig, moran_resid_not)
moran_resid <- moran_resid[!is.na(I_raw) & !is.na(I_resid)]

cat("\n=== Moran's I: Raw vs Residual (conditional on phenotype) ===\n")
summary_resid <- moran_resid[, .(
  mean_I_raw = mean(I_raw),
  mean_I_resid = mean(I_resid),
  delta_I = mean(I_raw - I_resid),
  pct_reduction = 100 * mean((I_raw - I_resid) / I_raw)
), by = sig_class]
print(summary_resid)

# Statistical test: is reduction greater for sig_both?
cat("\n=== Test: Does phenotype explain more autocorrelation in sig_both? ===\n")
moran_resid[, delta_I := I_raw - I_resid]

wt <- wilcox.test(delta_I ~ sig_class, data = moran_resid)
cat(sprintf("Wilcoxon test (delta_I): W=%.1f, p=%.4e\n", wt$statistic, wt$p.value))

tt <- t.test(delta_I ~ sig_class, data = moran_resid)
cat(sprintf("t-test (delta_I): t=%.2f, p=%.4e\n", tt$statistic, tt$p.value))

# ==== VISUALIZATIONS ====
cat("\nGenerating plots...\n")

# 1. Raw Moran's I by sig class
p1 <- ggplot(moran_raw, aes(x = sig_class, y = I, fill = sig_class)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sig Class", y = "Moran's I",
       title = "Phylogenetic Autocorrelation of Raw Genotypes") +
  theme_bw() + theme(legend.position = "none")

# 2. Raw vs Residual Moran's I
p2 <- ggplot(moran_resid, aes(x = I_raw, y = I_resid, color = sig_class)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Moran's I (raw genotype)", y = "Moran's I (residual | phenotype)",
       title = "Autocorrelation: Raw vs Phenotype-Conditioned") +
  theme_bw()

# 3. Delta I (reduction in autocorrelation)
p3 <- ggplot(moran_resid, aes(x = sig_class, y = delta_I, fill = sig_class)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sig Class", y = "ΔI (I_raw - I_resid)",
       title = "Reduction in Autocorrelation After Conditioning on Phenotype",
       subtitle = "Higher = phenotype explains more spatial pattern") +
  theme_bw() + theme(legend.position = "none")

# 4. Density of delta I
p4 <- ggplot(moran_resid, aes(x = delta_I, fill = sig_class, color = sig_class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "ΔI", y = "Density",
       title = "Distribution of Autocorrelation Reduction") +
  theme_bw()

library(patchwork)
p_combined <- (p1 | p2) / (p3 | p4)
ggsave("results/phylo_autocorr_moran.pdf", p_combined, width = 12, height = 10)

# ==== REPORT ====
cat("\n========================================\n")
cat("PHYLOGENETIC AUTOCORRELATION REPORT\n")
cat("========================================\n\n")

cat("HYPOTHESIS:\n")
cat("If sig_both sites are under spatially varying selection tracking climate,\n")
cat("then their phylogenetic autocorrelation should be EXPLAINED by phenotype.\n")
cat("i.e., Moran's I should decrease more when conditioning on phenotype.\n\n")

cat("RESULTS:\n")
cat("\nRaw Moran's I (phylogenetic clustering of genotypes):\n")
print(moran_raw[, .(mean_I = mean(I), sd_I = sd(I)), by = sig_class])

cat("\nAfter conditioning on phenotype:\n")
print(summary_resid)

cat("\nINTERPRETATION:\n")
sig_delta <- summary_resid[sig_class == "sig_both", delta_I]
not_delta <- summary_resid[sig_class == "not_sig", delta_I]

if (sig_delta > not_delta) {
  cat("*** sig_both sites show GREATER reduction in autocorrelation ***\n")
  cat("This supports spatially varying selection hypothesis:\n")
  cat("  - Genotypes cluster phylogenetically (shared ancestry)\n
")
  cat("  - But clustering is largely EXPLAINED by phenotype (climate)\n")
  cat("  - Consistent with local adaptation to temperature\n")
} else {
  cat("No evidence that phenotype explains more variance in sig_both sites.\n")
}

# ==== SAVE ====
saveRDS(list(
  moran_raw = moran_raw,
  moran_resid = moran_resid,
  summary_resid = summary_resid,
  tests = list(wilcox = wt, ttest = tt)
), "data/tmp/phylo_autocorr_analysis.rds")

cat("\nSaved: data/tmp/phylo_autocorr_analysis.rds\n")
cat("Saved: results/phylo_autocorr_moran.pdf\n")