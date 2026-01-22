# =============================================================================
# FIGURE 1: GWAS Results Overview
# Panels: A) Phenotype by order, B) Manhattan, C) Site classification, D) Complex enrichment
# =============================================================================

library(data.table)
library(arrow)
library(Biostrings)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

# Load GWAS results
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_results <- rbindlist(lapply(model_files, function(f) {
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
}))
stopifnot(nrow(gwas_results) > 0)

# -----------------------------------------------------------------------------
# CLASSIFY SITES
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# DEFINE COMPLEXES
# -----------------------------------------------------------------------------

complex_map <- c(
  # Photosystem I
  psaA = "PSI", psaB = "PSI", psaC = "PSI", psaI = "PSI", psaJ = "PSI",
  # Photosystem II
  psbA = "PSII", psbB = "PSII", psbC = "PSII", psbD = "PSII", psbE = "PSII",
  psbF = "PSII", psbH = "PSII", psbI = "PSII", psbJ = "PSII", psbK = "PSII",
  psbL = "PSII", psbM = "PSII", psbN = "PSII", psbT = "PSII", psbZ = "PSII",
  # Cytochrome b6f
  petA = "Cytb6f", petB = "Cytb6f", petD = "Cytb6f", petG = "Cytb6f",
  petL = "Cytb6f", petN = "Cytb6f",
  # ATP synthase
  atpA = "ATPase", atpB = "ATPase", atpE = "ATPase", atpF = "ATPase",
  atpH = "ATPase", atpI = "ATPase",
  # Rubisco
  rbcL = "rbcL",
  # NDH complex
  ndhA = "NDH", ndhB = "NDH", ndhC = "NDH", ndhD = "NDH", ndhE = "NDH",
  ndhF = "NDH", ndhG = "NDH", ndhH = "NDH", ndhI = "NDH", ndhJ = "NDH",
  ndhK = "NDH",
  # RNA polymerase
  rpoA = "RNAP", rpoB = "RNAP", rpoC1 = "RNAP", rpoC2 = "RNAP",
  # Ribosomal proteins
  rpl2 = "Ribosome", rpl14 = "Ribosome", rpl16 = "Ribosome", rpl20 = "Ribosome",

rpl22 = "Ribosome", rpl23 = "Ribosome", rpl32 = "Ribosome", rpl33 = "Ribosome",
  rpl36 = "Ribosome", rps2 = "Ribosome", rps3 = "Ribosome", rps4 = "Ribosome",
  rps7 = "Ribosome", rps8 = "Ribosome", rps11 = "Ribosome", rps12 = "Ribosome",
  rps14 = "Ribosome", rps15 = "Ribosome", rps16 = "Ribosome", rps18 = "Ribosome",
  rps19 = "Ribosome",
  # Other
  accD = "Other", cemA = "Other", clpP = "Other", ccsA = "Other",
  matK = "Other", ycf1 = "Other", ycf2 = "Other", ycf3 = "Other", ycf4 = "Other"
)

gwas_results[, Complex := complex_map[Gene]]

table(gwas_results$Complex, useNA = "ifany")
# -----------------------------------------------------------------------------
# PANEL A: Phenotype distribution by order
# -----------------------------------------------------------------------------

pheno_data <- data[, .(ID, pheno = get(pheno_col), Order)]
pheno_data <- pheno_data[!is.na(pheno) & !is.na(Order)]

order_stats <- pheno_data[, .(median_temp = median(pheno), n = .N), by = Order]
order_stats <- order_stats[n >= 20]
setorder(order_stats, median_temp)

pheno_data <- pheno_data[Order %in% order_stats$Order]
pheno_data[, Order := factor(Order, levels = order_stats$Order)]

length(unique(pheno_data$Order))

panel_A <- ggplot(pheno_data, aes(x = Order, y = pheno, fill = Order)) +
  geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = colorRampPalette(c("steelblue", "grey80", "firebrick"))(nrow(order_stats))) +
  geom_hline(yintercept = median(pheno_data$pheno), lty = 2, color = "grey40") +
  labs(x = NULL, y = "Mean temperature wettest quarter (°C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

panel_A
ggsave("results/figures/fig1_panel_A.pdf", panel_A, width = 8, height = 4)

# -----------------------------------------------------------------------------
# PANEL B: Manhattan plot (cleaner version)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PANEL B: Manhattan plot by index, colored by complex
# -----------------------------------------------------------------------------

gwas_results[, neglog10p := -log10(P_aa_with_pcs)]
gwas_results[is.infinite(neglog10p), neglog10p := NA]

# Index as x-axis
gwas_results[, idx := .I]

# Complex colors
complex_cols <- c(
  PSII = "#1b9e77",
  Cytb6f = "#d95f02", 
  PSI = "#7570b3",
  ATPase = "#e7298a",
  NDH = "#66a61e",
  rbcL = "#cornflowerblue",
  RNAP = "#a6761d",
  Ribosome = "#666666",
  Other = "#999999"
)

panel_B <- ggplot(gwas_results, aes(x = idx, y = neglog10p, color = Complex)) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_hline(yintercept = -log10(thresh_control), lty = 2, color = "firebrick", linewidth = 0.5) +
  scale_color_manual(values = complex_cols) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Site index", y = "-log10(P)", color = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.3),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1)))

ggsave("results/figures/figure1_panel_B.pdf", panel_B, width = 8, height = 3)

# -----------------------------------------------------------------------------
# PANEL C: Site classification scatter
# -----------------------------------------------------------------------------
gwas_results[, neglog10p_only := -log10(P_aa_only)]
gwas_results[, neglog10p_ctrl := -log10(P_aa_with_pcs)]
gwas_results[is.infinite(neglog10p_only), neglog10p_only := NA]
gwas_results[is.infinite(neglog10p_ctrl), neglog10p_ctrl := NA]

gwas_results[, sig_class := factor(sig_class, 
  levels = c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))]

cols <- c(not_sig = "grey80", sig_nocontrol = "steelblue", 
          sig_control = "darkorange", sig_both = "firebrick")
labels <- c(not_sig = "not_sig", sig_nocontrol = "sig_nocontrol",
            sig_control = "sig_control", sig_both = "sig_both")

counts <- gwas_results[, .N, by = sig_class]
counts[, label := paste0(labels[sig_class], " (n=", N, ")")]

panel_C <- ggplot(gwas_results, aes(x = neglog10p_only, y = neglog10p_ctrl, color = sig_class)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = -log10(thresh_control), lty = 2, color = "darkorange") +
  geom_vline(xintercept = -log10(thresh_nocontrol), lty = 2, color = "steelblue") +
  scale_color_manual(values = cols, labels = counts$label, name = NULL) +
  labs(x = "-log10(P) without pop. structure control",
       y = "-log10(P) with pop. structure control") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = NA))

ggsave("results/figures/figure1_panel_C.pdf", panel_C, width = 5, height = 5)
panel_C
# -----------------------------------------------------------------------------
# PANEL D: Complex enrichment
# -----------------------------------------------------------------------------

gwas_results[, neglog10p_ctrl := -log10(P_aa_with_pcs)]
gwas_results[, neglog10p_only := -log10(P_aa_only)]

# Quick test
kruskal.test(neglog10p_ctrl ~ Complex, data = gwas_results)

sig_ctrl <- gwas_results[P_aa_with_pcs < thresh_control]
sig_only <- gwas_results[P_aa_only < thresh_nocontrol]

# Tests on significant sites
kruskal.test(neglog10p_ctrl ~ Complex, data = sig_ctrl)
kruskal.test(neglog10p_only ~ Complex, data = sig_only)

# Boxplot - significant with control
panel_D <- ggplot(sig_ctrl, aes(x = reorder(Complex, -neglog10p_ctrl, FUN = median), 
                                 y = neglog10p_ctrl, fill = Complex)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = complex_cols, guide = "none") +
  labs(x = NULL, y = "-log10(P) with control", 
       subtitle = paste0("n = ", nrow(sig_ctrl), " significant sites")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot - significant without control
panel_D.2 <- ggplot(sig_only, aes(x = reorder(Complex, -neglog10p_only, FUN = median), 
                                   y = neglog10p_only, fill = Complex)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = complex_cols, guide = "none") +
  labs(x = NULL, y = "-log10(P) without control",
       subtitle = paste0("n = ", nrow(sig_only), " significant sites")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

panel_D <- panel_D / panel_D.2 + plot_annotation(title = "Significant sites by complex")
panel_D
ggsave("results/figures/figure1_panel_D.pdf", panel_D, width = 5, height = 4)



# -----------------------------------------------------------------------------
# PANEL D: Site classification by complex (Other expanded)
# -----------------------------------------------------------------------------

gwas_results[, display_group := ifelse(Complex == "Other", paste0("other_", Gene), Complex)]

group_stats <- gwas_results[, .N, by = .(display_group, sig_class)]
group_totals <- gwas_results[, .(n_total = .N), by = display_group]
group_stats <- merge(group_stats, group_totals, by = "display_group")

# Get complex for each display_group (for coloring labels)
group_complex <- gwas_results[, .(Complex = Complex[1]), by = display_group]

prop_order <- gwas_results[, .(
  prop_sig = mean(sig_class %in% c("sig_both", "sig_control")),
  is_complex = !grepl("^other_", display_group)
), by = display_group]
prop_order <- prop_order[!duplicated(display_group)]
setorder(prop_order, -is_complex, -prop_sig)

group_stats[, display_group := factor(display_group, levels = prop_order$display_group)]
group_stats[, sig_class := factor(sig_class, 
  levels = c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))]

# Clean labels
clean_labels <- gsub("^other_", "", prop_order$display_group)

# Get colors for each label based on complex
group_complex <- merge(group_complex, prop_order[, .(display_group)], by = "display_group")
group_complex[, display_group := factor(display_group, levels = prop_order$display_group)]
setorder(group_complex, display_group)

complex_cols <- c(
  PSII = "#1b9e77",
  Cytb6f = "#d95f02", 
  PSI = "#7570b3",
  ATPase = "#e7298a",
  NDH = "#66a61e",
  rbcL = "#6495ED",
  RNAP = "#a6761d",
  Ribosome = "#666666",
  Other = "#999999"
)

label_colors <- complex_cols[group_complex$Complex]

panel_D <- ggplot(group_stats, aes(x = display_group, y = N, fill = sig_class)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols, labels = labels, name = NULL) +
  scale_x_discrete(labels = clean_labels) +
  labs(x = NULL, y = "Number of sites") +
  theme_minimal() +
  facet_wrap(~ sig_class, scales="free") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, 
                                   color = label_colors),
        panel.grid.major.x = element_blank())

ggsave("results/figures/figure1_panel_D_barplot.pdf", panel_D, width = 6, height = 4)
panel_D

# -----------------------------------------------------------------------------
# PANEL D: Heatmap of sig_class by group, colored labels by complex
# -----------------------------------------------------------------------------

# Create display group: use Complex for major groups, Gene for Other
gwas_results[, display_group := ifelse(Complex == "Other", Gene, Complex)]

# Proportion of each sig_class within each display_group
heat_data <- gwas_results[, .N, by = .(display_group, Complex, sig_class)]
heat_data[, prop := N / sum(N), by = display_group]

# Order: complexes first by prop sig, then Other genes by prop sig
group_order <- gwas_results[, .(
  prop_sig = mean(sig_class %in% c("sig_both", "sig_control")),
  Complex = Complex[1],
  is_complex = Complex[1] != "Other"
), by = display_group]
group_order <- group_order[!duplicated(display_group)]
setorder(group_order, -is_complex, -prop_sig)

heat_data[, display_group := factor(display_group, levels = group_order$display_group)]
heat_data[, sig_class := factor(sig_class, 
  levels = c("not_sig", "sig_nocontrol", "sig_control", "sig_both"))]

# Label colors by complex
complex_cols <- c(
  PSII = "#1b9e77",
  Cytb6f = "#d95f02", 
  PSI = "#7570b3",
  ATPase = "#e7298a",
  NDH = "#66a61e",
  rbcL = "#6495ED",
  RNAP = "#a6761d",
  Ribosome = "#666666",
  Other = "#999999"
)
label_colors <- complex_cols[group_order$Complex]

panel_D <- ggplot(heat_data, aes(x = sig_class, y = display_group, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = N), size = 2.5) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Proportion") +
  scale_x_discrete(labels = labels) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = label_colors))

ggsave("results/figures/figure1_panel_D.pdf", panel_D, width = 5, height = 5)
# -----------------------------------------------------------------------------
# SUPPLEMENTARY: QQ plot
# -----------------------------------------------------------------------------

p_vals <- gwas_results$P_aa_with_pcs
p_vals <- p_vals[!is.na(p_vals) & p_vals > 0]
lambda <- median(qchisq(1 - p_vals, df = 1)) / qchisq(0.5, df = 1)

qq_data <- data.table(
  expected = -log10(ppoints(length(p_vals))),
  observed = -log10(sort(p_vals))
)

panel_qq <- ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 0.5, color = "steelblue", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "firebrick") +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)",
       title = sprintf("λ = %.2f", lambda)) +
  theme_minimal()

ggsave("results/figures/fig1_panel_qq.pdf", panel_qq, width = 4, height = 4)
panel_qq


# -----------------------------------------------------------------------------
# GENERATE FIGURE
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# SUMMARY STATISTICS FOR TEXT
# -----------------------------------------------------------------------------

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total sites tested:", nrow(gwas_results), "\n")
cat("Sites by class:\n")
print(table(gwas_results$sig_class))
cat("\nThresholds used:\n")
cat("  P_aa_with_pcs 5%:", thresh_control, "\n")
cat("  P_aa_only 20%:", thresh_nocontrol, "\n")
cat("\nGenomic inflation factor (λ):", 
    round(median(qchisq(1 - gwas_results$P_aa_with_pcs, df = 1), na.rm = TRUE) / 
            qchisq(0.5, df = 1), 2), "\n")