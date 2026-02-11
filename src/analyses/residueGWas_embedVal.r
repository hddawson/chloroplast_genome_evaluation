library(data.table)
library(ggplot2)
library(scales)

# ---- LOAD RESULTS ----
results <- fread("results/q3_gwas/q3_run_associations_nopca.csv")
stopifnot(nrow(results) > 0)

struct_runs <- readRDS("results/q3_gwas/struct_runs.rds")
stopifnot(nrow(struct_runs) > 0)

dir.create("results/q3_gwas/figures", showWarnings = FALSE, recursive = TRUE)

n_sites <- nrow(results)
bonf <- 0.05 / n_sites

message("Sites: ", n_sites, " | Bonferroni: ", signif(bonf, 3))

# =====================================================================
# 1. QQ PLOTS + GENOMIC INFLATION FACTOR (LAMBDA)
# =====================================================================

calc_lambda <- function(p) {
  chisq <- qchisq(1 - p[!is.na(p)], df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

lambda_corr <- calc_lambda(results$P_corrected)
lambda_uncorr <- calc_lambda(results$P_uncorrected)
message("Lambda corrected: ", round(lambda_corr, 3))
message("Lambda uncorrected: ", round(lambda_uncorr, 3))

make_qq_dt <- function(p, label) {
  p <- sort(p[!is.na(p)])
  n <- length(p)
  data.table(
    observed = -log10(p),
    expected = -log10(ppoints(n)),
    label = label
  )
}

qq_dt <- rbind(
  make_qq_dt(results$P_corrected, paste0("Corrected (λ=", round(lambda_corr, 3), ")")),
  make_qq_dt(results$P_uncorrected, paste0("Uncorrected (λ=", round(lambda_uncorr, 3), ")"))
)

p_qq <- ggplot(qq_dt, aes(x = expected, y = observed, color = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("steelblue", "firebrick")) +
  labs(x = expression(Expected ~ -log[10](P)),
       y = expression(Observed ~ -log[10](P)),
       title = "QQ Plot: Q3 Structural Run Associations",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = c(0.3, 0.85),
        legend.background = element_rect(fill = "white", color = NA))

ggsave("results/q3_gwas/figures/qq_plot.png", p_qq, width = 7, height = 7, dpi = 300)
ggsave("results/q3_gwas/figures/qq_plot.pdf", p_qq, width = 7, height = 7)

# =====================================================================
# 2. PER-GENE SIGNIFICANT SITE COUNTS BY Q3 TYPE
# =====================================================================

results[, sig_corr := P_corrected < bonf]

gene_q3_counts <- results[sig_corr == TRUE, .N, by = .(Gene, consensus_q3)]
gene_totals <- results[, .(n_tested = .N), by = Gene]
gene_sig_totals <- results[sig_corr == TRUE, .(n_sig = .N), by = Gene]

# Order genes by total significant sites
gene_order <- gene_sig_totals[order(-n_sig), Gene]
gene_q3_counts[, Gene := factor(Gene, levels = gene_order)]

p_gene_q3 <- ggplot(gene_q3_counts, aes(x = Gene, y = N, fill = consensus_q3)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set1", name = "Q3 Class") +
  labs(title = "Significant Q3 Runs per Gene (Bonferroni)",
       x = NULL, y = "Number of Significant Runs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("results/q3_gwas/figures/sig_sites_per_gene.png", p_gene_q3, width = 10, height = 5, dpi = 300)
ggsave("results/q3_gwas/figures/sig_sites_per_gene.pdf", p_gene_q3, width = 10, height = 5)

# =====================================================================
# 3. Q3 ENRICHMENT TEST (FISHER'S EXACT)
# =====================================================================

# Contingency: significant vs not, by Q3 class
q3_table <- results[, .(
  n_sig = sum(sig_corr),
  n_nonsig = sum(!sig_corr)
), by = consensus_q3]

message("\n=== Q3 class counts ===")
print(q3_table)

# Overall enrichment: chi-squared test
contingency_mat <- as.matrix(q3_table[, .(n_sig, n_nonsig)])
rownames(contingency_mat) <- q3_table$consensus_q3
chisq_test <- chisq.test(contingency_mat)
message("Chi-squared test for Q3 enrichment: p = ", signif(chisq_test$p.value, 3))

# Pairwise Fisher's exact: each Q3 class vs rest
q3_classes <- unique(results$consensus_q3)
enrichment_results <- rbindlist(lapply(q3_classes, function(q3) {
  is_class <- results$consensus_q3 == q3
  mat <- matrix(c(
    sum(is_class & results$sig_corr),
    sum(is_class & !results$sig_corr),
    sum(!is_class & results$sig_corr),
    sum(!is_class & !results$sig_corr)
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("this_q3", "other"), c("sig", "nonsig")))
  
  ft <- fisher.test(mat)
  data.table(
    Q3 = q3,
    n_sig = mat[1, 1],
    n_total = sum(is_class),
    prop_sig = mat[1, 1] / sum(is_class),
    odds_ratio = ft$estimate,
    p_value = ft$p.value,
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2]
  )
}))

message("\n=== Q3 Enrichment (Fisher's exact, each vs rest) ===")
print(enrichment_results)
fwrite(enrichment_results, "results/q3_gwas/q3_enrichment_fisher.csv")

# Forest plot of odds ratios
p_enrichment <- ggplot(enrichment_results, aes(x = Q3, y = odds_ratio, ymin = ci_low, ymax = ci_high)) +
  geom_pointrange(size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(title = "Q3 Class Enrichment Among Significant Sites",
       x = "Q3 Class", y = "Odds Ratio (vs all other classes)") +
  theme_minimal()

ggsave("results/q3_gwas/figures/q3_enrichment_or.png", p_enrichment, width = 6, height = 5, dpi = 300)
ggsave("results/q3_gwas/figures/q3_enrichment_or.pdf", p_enrichment, width = 6, height = 5)

# =====================================================================
# 4. EFFECT SIZE (R2) DISTRIBUTION BY Q3 CLASS
# =====================================================================

p_r2_q3 <- ggplot(results, aes(x = consensus_q3, y = R2_corrected, fill = consensus_q3)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "Effect Size by Q3 Class (Population-Corrected R²)",
       x = "Q3 Class", y = "R² (corrected)") +
  theme_minimal()

ggsave("results/q3_gwas/figures/r2_by_q3.png", p_r2_q3, width = 6, height = 5, dpi = 300)
ggsave("results/q3_gwas/figures/r2_by_q3.pdf", p_r2_q3, width = 6, height = 5)

# KW test for R2 differences across Q3
kw_r2 <- kruskal.test(R2_corrected ~ consensus_q3, data = results)
message("Kruskal-Wallis R2 ~ Q3: p = ", signif(kw_r2$p.value, 3))

# =====================================================================
# 5. GENE-LEVEL MANHATTAN WITH STRUCTURE TRACK (PANEL FIGURE)
# =====================================================================

# Recompute ordering
setorder(results, Gene, run_id)
results[, x_index := .I]
results[, gene_idx := as.integer(factor(Gene, levels = unique(Gene)))]
gene_labels <- results[, .(x_mid = median(x_index), x_min = min(x_index), x_max = max(x_index)), by = Gene]

# Top panel: Manhattan
p_top <- ggplot(results, aes(x = x_index, y = -log10(P_corrected), color = factor(gene_idx %% 2))) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_color_manual(values = c("0" = "grey30", "1" = "steelblue"), guide = "none") +
  scale_x_continuous(breaks = gene_labels$x_mid, labels = gene_labels$Gene) +
  labs(y = expression(-log[10](P)), x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(5, 5, 0, 5))

# Bottom panel: Q3 structure track
q3_colors <- c("C" = "#E41A1C", "E" = "#377EB8", "H" = "#4DAF4A")

p_bottom <- ggplot(results, aes(x = x_index, y = 1, fill = consensus_q3)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = q3_colors, name = "Q3") +
  scale_x_continuous(breaks = gene_labels$x_mid, labels = gene_labels$Gene) +
  labs(x = "Gene", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid = element_blank(),
        plot.margin = margin(0, 5, 5, 5))

library(patchwork)
p_combined <- p_top / p_bottom + plot_layout(heights = c(4, 1)) +
  plot_annotation(title = "Q3 Structural Run Associations with Temperature")

ggsave("results/q3_gwas/figures/manhattan_with_track.png", p_combined, width = 14, height = 7, dpi = 300)
ggsave("results/q3_gwas/figures/manhattan_with_track.pdf", p_combined, width = 14, height = 7)

# =====================================================================
# 6. OVERLAP WITH SITE-LEVEL GWAS
# =====================================================================

# Load site-level results (from npc_experiment at n_pcs=1000)
site_gwas_dir <- "results/npc_experiment/npc_1000"
if (dir.exists(site_gwas_dir)) {
  site_rds <- list.files(site_gwas_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(site_rds) > 0) {
    site_results <- rbindlist(lapply(site_rds, function(f) {
      res_list <- readRDS(f)
      rbindlist(lapply(res_list, function(r) {
        data.table(
          Gene = r$Gene,
          Position = r$Aligned_Position,
          P_site = r$P_aa_with_pcs,
          R2_site = r$R2_full
        )
      }))
    }), fill = TRUE)
    
    site_results[, sig_site := P_site < 0.05 / nrow(site_results)]
    
    # Map Q3 run membership to site positions using struct_runs
    site_in_runs <- merge(
      site_results,
      struct_runs[, .(Gene, Residue_Index, site_id)],
      by.x = c("Gene", "Position"),
      by.y = c("Gene", "Residue_Index"),
      all = FALSE
    )
    
    if (nrow(site_in_runs) > 0) {
      # Per-run: is any constituent site significant?
      run_site_overlap <- site_in_runs[, .(
        any_site_sig = any(sig_site),
        min_site_p = min(P_site, na.rm = TRUE),
        n_sites_in_run = .N
      ), by = site_id]
      
      overlap_dt <- merge(results, run_site_overlap, by = "site_id", all.x = TRUE)
      overlap_dt[is.na(any_site_sig), any_site_sig := FALSE]
      
      # Concordance table
      conc_table <- overlap_dt[, table(sig_corr, any_site_sig)]
      message("\n=== Concordance: Q3 run sig vs any site sig ===")
      print(conc_table)
      
      ft_conc <- fisher.test(conc_table)
      message("Fisher's exact OR: ", round(ft_conc$estimate, 2), 
              ", p = ", signif(ft_conc$p.value, 3))
      
      # Scatter: run p vs best site p within run
      overlap_dt <- overlap_dt[!is.na(min_site_p)]
      
      p_overlap <- ggplot(overlap_dt, aes(x = -log10(min_site_p), y = -log10(P_corrected), 
                                          color = consensus_q3)) +
        geom_point(alpha = 0.5, size = 1.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        scale_color_brewer(palette = "Set1", name = "Q3") +
        labs(title = "Q3 Run vs Best Site-Level P-value",
             x = expression(-log[10](P[site])),
             y = expression(-log[10](P[Q3~run]))) +
        theme_minimal()
      
      ggsave("results/q3_gwas/figures/run_vs_site_overlap.png", p_overlap, width = 7, height = 6, dpi = 300)
      ggsave("results/q3_gwas/figures/run_vs_site_overlap.pdf", p_overlap, width = 7, height = 6)
    } else {
      message("No position overlap between site GWAS and Q3 runs.")
    }
  }
} else {
  message("Site-level GWAS results not found at ", site_gwas_dir, "; skipping overlap.")
}

# =====================================================================
# 7. PERMUTATION TEST FOR Q3 ENRICHMENT
# =====================================================================

n_perms <- 1000L
set.seed(42)

# Observed: proportion significant per Q3 class
obs_prop <- results[, .(prop_sig = mean(sig_corr)), by = consensus_q3]
obs_stat <- var(obs_prop$prop_sig)  # variance across Q3 classes

perm_stats <- vapply(seq_len(n_perms), function(i) {
  perm_sig <- sample(results$sig_corr)  # shuffle significance labels
  perm_prop <- tapply(perm_sig, results$consensus_q3, mean)
  var(perm_prop)
}, numeric(1))

perm_p <- mean(perm_stats >= obs_stat)
message("\n=== Permutation test for Q3 enrichment ===")
message("Observed variance in sig proportion across Q3: ", signif(obs_stat, 4))
message("Permutation p-value (", n_perms, " perms): ", perm_p)

p_perm <- ggplot(data.table(stat = perm_stats), aes(x = stat)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40") +
  geom_vline(xintercept = obs_stat, color = "red", linewidth = 1) +
  labs(title = paste0("Permutation Test: Q3 Enrichment (p=", round(perm_p, 3), ")"),
       x = "Variance in Proportion Significant Across Q3 Classes",
       y = "Count") +
  theme_minimal()

ggsave("results/q3_gwas/figures/permutation_q3_enrichment.png", p_perm, width = 7, height = 5, dpi = 300)
ggsave("results/q3_gwas/figures/permutation_q3_enrichment.pdf", p_perm, width = 7, height = 5)

# =====================================================================
# 8. SUPPLEMENTARY TABLE
# =====================================================================

# Top hits table for supplement
top_hits <- results[P_corrected < bonf][order(P_corrected)]
top_hits[, P_corrected_sci := formatC(P_corrected, format = "e", digits = 2)]
top_hits[, P_uncorrected_sci := formatC(P_uncorrected, format = "e", digits = 2)]

fwrite(top_hits, "results/q3_gwas/significant_hits_table.csv")
message("\nSignificant hits saved: ", nrow(top_hits))

# =====================================================================
# SUMMARY STATS FOR MANUSCRIPT
# =====================================================================

message("\n========== MANUSCRIPT SUMMARY STATS ==========")
message("Total Q3 runs tested: ", n_sites)
message("Bonferroni threshold: ", signif(bonf, 3))
message("Lambda (corrected): ", round(lambda_corr, 3))
message("Lambda (uncorrected): ", round(lambda_uncorr, 3))
message("Significant runs (corrected): ", sum(results$sig_corr))
message("Significant runs (uncorrected): ", sum(results$P_uncorrected < bonf, na.rm = TRUE))
message("Chi-squared Q3 enrichment p: ", signif(chisq_test$p.value, 3))
message("Permutation Q3 enrichment p: ", perm_p)
message("KW R2 ~ Q3 class p: ", signif(kw_r2$p.value, 3))
message("Genes with significant hits: ", 
        paste(unique(top_hits$Gene), collapse = ", "))

# Save all stats
saveRDS(list(
  lambda_corrected = lambda_corr,
  lambda_uncorrected = lambda_uncorr,
  enrichment_fisher = enrichment_results,
  enrichment_chisq = chisq_test,
  enrichment_permutation_p = perm_p,
  kw_r2_test = kw_r2,
  n_sig_corrected = sum(results$sig_corr),
  n_sig_uncorrected = sum(results$P_uncorrected < bonf, na.rm = TRUE)
), "results/q3_gwas/manuscript_stats.rds")

message("\n=== Done ===")