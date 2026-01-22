# ==== SITE VARIANCE ANALYSIS BY SIG CLASS ====
library(ape)
library(data.table)
library(ggplot2)
library(patchwork)


# ==== CONFIGURATION ====
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

# ==== LOAD DATA ====
cat("Loading data...\n")
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)

partition_map <- readRDS(partition_map_file)
setDT(partition_map)

# ==== LOAD AND CLASSIFY GWAS EFFECTS ====
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

cat("Sites by class:\n")
print(gwas_effects[, .N, by = sig_class][order(-N)])

# ==== CALCULATE SITE-LEVEL VARIANCE METRICS ====
cat("\nCalculating site variance metrics...\n")

valid_aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

calc_site_stats <- function(site_col) {
  # Remove gaps/ambiguous
  aa <- site_col[site_col %in% valid_aa]
  n <- length(aa)
  if (n < 10) return(list(n_valid = n, pi = NA, n_segregating = NA, 
                          theta_w = NA, tajima_d_proxy = NA, 
                          heterozygosity = NA, gap_freq = NA))
  
  # Allele frequencies
  
  tab <- table(aa)
  freqs <- as.numeric(tab) / n
  k <- length(tab)  # number of alleles (unique AAs)
  
  # 1. Nucleotide diversity analog (pi) - expected heterozygosity
  # pi = 1 - sum(p_i^2), where p_i is frequency of allele i
  pi <- 1 - sum(freqs^2)
  
  
  # 2. Number of segregating sites (polymorphic = >1 allele)
  n_segregating <- as.integer(k > 1)
  
  # 3. Watterson's theta analog
  # theta_w = S / a_n, where a_n = sum(1/i) for i=1 to n-1
  # For a single site, S = 1 if polymorphic, 0 otherwise
  # We use k-1 as proxy for number of mutations
  a_n <- sum(1 / (1:(n-1)))
  theta_w <- (k - 1) / a_n
  
  # 4. Tajima's D proxy: (pi - theta_w) / SE
  # Simplified: just the difference, since SE calculation is complex
  tajima_d_proxy <- pi - theta_w
  
  # 5. Observed heterozygosity (proportion of pairwise differences)
  # Same as pi for haploid
  heterozygosity <- pi
  
  # 6. Gap frequency
  gap_freq <- 1 - (n / length(site_col))
  
  # 7. Minor allele frequency (MAF)
  if (k > 1) {
    maf <- sort(freqs, decreasing = TRUE)[2]  # second most common
  } else {
    maf <- 0
  }
  
  # 8. Effective number of alleles: n_e = 1 / sum(p_i^2)
  n_eff <- 1 / sum(freqs^2)
  
  list(n_valid = n, pi = pi, n_segregating = n_segregating,
       theta_w = theta_w, tajima_d_proxy = tajima_d_proxy,
       heterozygosity = heterozygosity, gap_freq = gap_freq,
       maf = maf, n_alleles = k, n_eff_alleles = n_eff)
}

# Calculate for all sites
site_stats <- rbindlist(lapply(1:nrow(gwas_effects), function(i) {
  if (i %% 1000 == 0) cat("  Site", i, "/", nrow(gwas_effects), "\n")
  
  pos <- gwas_effects$GlobalPos[i]
  if (pos > ncol(aln_mat)) return(NULL)
  
  stats <- calc_site_stats(aln_mat[, pos])
  data.table(
    GlobalPos = pos,
    Gene = gwas_effects$Gene[i],
    sig_class = gwas_effects$sig_class[i],
    n_valid = stats$n_valid,
    pi = stats$pi,
    theta_w = stats$theta_w,
    tajima_d_proxy = stats$tajima_d_proxy,
    n_segregating = stats$n_segregating,
    heterozygosity = stats$heterozygosity,
    maf = stats$maf,
    n_alleles = stats$n_alleles,
    n_eff_alleles = stats$n_eff_alleles,
    gap_freq = stats$gap_freq
  )
}))

site_stats <- site_stats[!is.na(pi)]
cat("Sites with valid stats:", nrow(site_stats), "\n")

# ==== STATISTICAL TESTS ====
cat("\n=== Population Genetic Metrics by Sig Class ===\n")

# Summary stats
summary_stats <- site_stats[, .(
  n_sites = .N,
  mean_pi = mean(pi),
  sd_pi = sd(pi),
  mean_theta_w = mean(theta_w),
  mean_tajima_d = mean(tajima_d_proxy),
  mean_maf = mean(maf),
  mean_n_alleles = mean(n_alleles),
  mean_n_eff = mean(n_eff_alleles),
  pct_polymorphic = 100 * mean(n_segregating)
), by = sig_class][order(sig_class)]

print(summary_stats)

# Kruskal-Wallis tests
cat("\n=== Kruskal-Wallis Tests ===\n")

kw_pi <- kruskal.test(pi ~ sig_class, data = site_stats)
kw_theta <- kruskal.test(theta_w ~ sig_class, data = site_stats)
kw_tajima <- kruskal.test(tajima_d_proxy ~ sig_class, data = site_stats)
kw_maf <- kruskal.test(maf ~ sig_class, data = site_stats)
kw_neff <- kruskal.test(n_eff_alleles ~ sig_class, data = site_stats)

cat("π (nucleotide diversity):", sprintf("χ²=%.2f, p=%.2e\n", kw_pi$statistic, kw_pi$p.value))
cat("θ_W (Watterson's theta):", sprintf("χ²=%.2f, p=%.2e\n", kw_theta$statistic, kw_theta$p.value))
cat("Tajima's D proxy:", sprintf("χ²=%.2f, p=%.2e\n", kw_tajima$statistic, kw_tajima$p.value))
cat("MAF:", sprintf("χ²=%.2f, p=%.2e\n", kw_maf$statistic, kw_maf$p.value))
cat("N_eff alleles:", sprintf("χ²=%.2f, p=%.2e\n", kw_neff$statistic, kw_neff$p.value))

# Pairwise Wilcoxon
cat("\n=== Pairwise Wilcoxon Tests (π) ===\n")
pw_pi <- pairwise.wilcox.test(site_stats$pi, site_stats$sig_class, 
                              p.adjust.method = "bonferroni")
print(pw_pi)

cat("\n=== Pairwise Wilcoxon Tests (Tajima's D proxy) ===\n")
pw_tajima <- pairwise.wilcox.test(site_stats$tajima_d_proxy, site_stats$sig_class,
                                  p.adjust.method = "bonferroni")
print(pw_tajima)

# ==== BY-GENE ANALYSIS ====
cat("\n=== π by Gene and Sig Class ===\n")

gene_pi <- site_stats[, .(
  mean_pi = mean(pi),
  mean_tajima = mean(tajima_d_proxy),
  n_sites = .N
), by = .(Gene, sig_class)]

gene_kw <- site_stats[, {
  if (length(unique(sig_class)) > 1 && .N > 10) {
    fit <- kruskal.test(pi ~ sig_class, data = .SD)
    .(kw_stat = fit$statistic, kw_p = fit$p.value)
  } else {
    .(kw_stat = NA_real_, kw_p = NA_real_)
  }
}, by = Gene]

cat("\nGenes with significant π differences by sig_class:\n")
print(gene_kw[kw_p < 0.05][order(kw_p)])

# ==== VISUALIZATIONS ====
cat("\nGenerating plots...\n")

# 1. Boxplots
p1 <- ggplot(site_stats, aes(x = sig_class, y = pi, fill = sig_class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  labs(x = "Sig Class", y = "π (nucleotide diversity analog)", 
       title = "Site Diversity (π) by Significance Class") +
  theme_bw() + theme(legend.position = "none")

p2 <- ggplot(site_stats, aes(x = sig_class, y = theta_w, fill = sig_class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  labs(x = "Sig Class", y = "θ_W (Watterson's theta analog)",
       title = "Mutation Rate Proxy (θ_W) by Significance Class") +
  theme_bw() + theme(legend.position = "none")

p3 <- ggplot(site_stats, aes(x = sig_class, y = tajima_d_proxy, fill = sig_class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Sig Class", y = "Tajima's D proxy (π - θ_W)",
       title = "Selection Signal (Tajima's D proxy)") +
  theme_bw() + theme(legend.position = "none")

p4 <- ggplot(site_stats, aes(x = sig_class, y = n_eff_alleles, fill = sig_class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  labs(x = "Sig Class", y = "Effective Number of Alleles",
       title = "Allelic Diversity by Significance Class") +
  theme_bw() + theme(legend.position = "none")

p_combined <- (p1 | p2) / (p3 | p4)
ggsave("results/popgen_variance_boxplots.pdf", p_combined, width = 10, height = 8)

# 2. Density plots
p_dens_pi <- ggplot(site_stats, aes(x = pi, fill = sig_class, color = sig_class)) +
  geom_density(alpha = 0.3) +
  labs(x = "π", y = "Density", title = "Distribution of Site Diversity (π)") +
  theme_bw()

p_dens_tajima <- ggplot(site_stats, aes(x = tajima_d_proxy, fill = sig_class, color = sig_class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Tajima's D proxy", y = "Density", title = "Distribution of Selection Signal") +
  theme_bw()

p_dens <- p_dens_pi / p_dens_tajima
ggsave("results/popgen_density.pdf", p_dens, width = 8, height = 8)

# 3. Pi vs Theta scatter (classic Tajima plot)
p_scatter <- ggplot(site_stats[n_segregating == 1], aes(x = theta_w, y = pi, color = sig_class)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ sig_class) +
  labs(x = "θ_W", y = "π", 
       title = "π vs θ_W (polymorphic sites only)\nAbove line = excess common variants, Below = excess rare variants") +
  theme_bw() + theme(legend.position = "none")

ggsave("results/pi_vs_theta_scatter.pdf", p_scatter, width = 10, height = 8)

# 4. By-gene faceted
p_gene <- ggplot(site_stats, aes(x = pi, fill = sig_class, color = sig_class)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(x = "π", y = "Density", title = "Site Diversity by Gene") +
  theme_bw() + theme(legend.position = "bottom")

ggsave("results/pi_by_gene.pdf", p_gene, width = 14, height = 12)

# 5. MAF spectrum
p_maf <- ggplot(site_stats[maf > 0], aes(x = maf, fill = sig_class)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity") +
  facet_wrap(~ sig_class, scales = "free_y") +
  labs(x = "Minor Allele Frequency", y = "Count",
       title = "Site Frequency Spectrum by Significance Class") +
  theme_bw() + theme(legend.position = "none")

ggsave("results/maf_spectrum.pdf", p_maf, width = 10, height = 8)

# ==== REPORT ====
cat("\n========================================\n")
cat("POPULATION GENETIC VARIANCE REPORT\n")
cat("========================================\n\n")

cat("METRICS:\n")
cat("  π: Nucleotide diversity (expected heterozygosity) = 1 - Σp_i²\n")
cat("  θ_W: Watterson's theta (mutation rate proxy) = (k-1)/a_n\n")
cat("  Tajima's D proxy: π - θ_W\n")
cat("    Positive = excess intermediate-frequency variants (balancing selection)\n")
cat("    Negative = excess rare variants (purifying selection or expansion)\n")
cat("  N_eff: Effective number of alleles = 1/Σp_i²\n\n")

cat("SUMMARY BY SIG CLASS:\n")
for (sc in c("sig_both", "sig_control", "sig_nocontrol", "not_sig")) {
  ss <- summary_stats[sig_class == sc]
  if (nrow(ss) > 0) {
    cat(sprintf("\n%s (n=%d sites, %.1f%% polymorphic):\n", 
                sc, ss$n_sites, ss$pct_polymorphic))
    cat(sprintf("  π: %.4f ± %.4f\n", ss$mean_pi, ss$sd_pi))
    cat(sprintf("  θ_W: %.4f\n", ss$mean_theta_w))
    cat(sprintf("  Tajima's D proxy: %.4f\n", ss$mean_tajima_d))
    cat(sprintf("  N_eff alleles: %.2f\n", ss$mean_n_eff))
  }
}

cat("\n\nINTERPRETATION:\n")
if (kw_tajima$p.value < 0.05) {
  cat("*** Significant difference in Tajima's D proxy between sig classes ***\n")
  tajima_means <- summary_stats[, .(sig_class, mean_tajima_d)]
  cat("\nTajima's D ranking (most positive = most balancing selection signal):\n")
  print(tajima_means[order(-mean_tajima_d)])
}

# ==== SAVE ====
saveRDS(list(
  site_stats = site_stats,
  summary_stats = summary_stats,
  gene_kw = gene_kw,
  tests = list(kw_pi = kw_pi, kw_theta = kw_theta, 
               kw_tajima = kw_tajima, kw_maf = kw_maf)
), "data/tmp/popgen_variance_analysis.rds")

cat("\nSaved: data/tmp/popgen_variance_analysis.rds\n")
cat("Saved: results/popgen_variance_boxplots.pdf\n")
cat("Saved: results/popgen_density.pdf\n")
cat("Saved: results/pi_vs_theta_scatter.pdf\n")
cat("Saved: results/pi_by_gene.pdf\n")
cat("Saved: results/maf_spectrum.pdf\n")