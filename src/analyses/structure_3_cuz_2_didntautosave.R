# structural_enrichment_analysis.R
# Analyze secondary structure enrichment in GWAS hits

library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)


# Compare distributions across site classes
par(mfrow = c(2, 3))

boxplot(mean_rsa.x ~ site_class, data = high_cons, 
        main = "Relative Solvent Accessibility", ylab = "RSA",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_asa.x ~ site_class, data = high_cons,
        main = "Absolute Solvent Accessibility", ylab = "ASA (Å²)",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_disorder.x ~ site_class, data = high_cons,
        main = "Disorder Probability", ylab = "P(disorder)",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_phi ~ site_class, data = high_cons,
        main = "Phi Angle", ylab = "φ (degrees)",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_psi ~ site_class, data = high_cons,
        main = "Psi Angle", ylab = "ψ (degrees)",
        col = c("gray70", "gold", "steelblue", "darkred"))

par(mfrow = c(1, 1))

# Statistical tests (Wilcoxon vs not_sig background)
cat("\n=== Wilcoxon tests vs not_sig ===\n")
for (feat in c("mean_rsa.x", "mean_asa.x", "mean_disorder.x", "mean_phi", "mean_psi")) {
  cat("\n", feat, ":\n", sep = "")
  bg <- high_cons[site_class == "not_sig", get(feat)]
  
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- high_cons[site_class == sc, get(feat)]
    if (length(test_vals) < 5) next
    
    wt <- wilcox.test(test_vals, bg)
    direction <- ifelse(median(test_vals, na.rm=T) > median(bg, na.rm=T), "↑", "↓")
    cat(sprintf("  %s: median=%.3f %s (bg=%.3f), p=%.2e\n", 
                sc, median(test_vals, na.rm=T), direction, median(bg, na.rm=T), wt$p.value))
  }
}

summary(high_cons)

# ---- pca ---- 

# Cleaner PCA - just biophysical features
pca_cols2 <- c("mean_rsa.x", "mean_asa.x", "mean_disorder.x", "mean_phi", "mean_psi")

pca_data2 <- high_cons[, ..pca_cols2]
pca_data2 <- pca_data2[complete.cases(pca_data2)]

pca_res2 <- prcomp(pca_data2, scale. = TRUE)
summary(pca_res2)

# Loadings
loadings2 <- as.data.frame(pca_res2$rotation[, 1:3])
loadings2$feature <- rownames(loadings2)
print(loadings2[order(abs(loadings2$PC1), decreasing = TRUE), ])

# Add PCs
high_cons_pca2 <- high_cons[complete.cases(high_cons[, ..pca_cols2])]
high_cons_pca2$PC1 <- pca_res2$x[, 1]
high_cons_pca2$PC2 <- pca_res2$x[, 2]

high_cons_pca2$site_class <- factor(high_cons_pca2$site_class,
                                    levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))

# Plot with density contours instead of ellipses for clarity
p2 <- ggplot(high_cons_pca2, aes(x = PC1, y = PC2, color = site_class)) +
  geom_point(alpha = 0.2, size = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("not_sig" = "gray50",
                                "sig_no_control" = "gold3",
                                "sig_with_control" = "steelblue",
                                "sig_both" = "darkred")) +
  labs(title = "PCA of Biophysical Features",
       subtitle = "RSA, ASA, disorder, phi, psi",
       x = paste0("PC1 (", round(summary(pca_res2)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res2)$importance[2,2]*100, 1), "%)")) +
  theme_classic() +
  theme(legend.position = "top")

p2

# Test if PC1 differs by site class
cat("\n=== PC1 by site class (Wilcoxon vs not_sig) ===\n")
bg_pc1 <- high_cons_pca2[site_class == "not_sig", PC1]
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_pc1 <- high_cons_pca2[site_class == sc, PC1]
  wt <- wilcox.test(test_pc1, bg_pc1)
  cat(sprintf("%s: median=%.2f (bg=%.2f), p=%.2e\n", 
              sc, median(test_pc1), median(bg_pc1), wt$p.value))
}
# Select numeric columns for PCA
pca_cols <- c("mean_rsa", "mean_asa", "mean_disorder", "mean_phi", "mean_psi",
              "prop_H", "prop_E", "prop_C", "entropy", "consensus_freq")

pca_data <- high_cons[, ..pca_cols]
pca_data <- pca_data[complete.cases(pca_data)]

# Scale and run PCA
pca_res <- prcomp(pca_data, scale. = TRUE, center = T)

# Variance explained
summary(pca_res)

# Add PCs back to data
high_cons_pca <- high_cons[complete.cases(high_cons[, ..pca_cols])]
high_cons_pca$PC1 <- pca_res$x[, 1]
high_cons_pca$PC2 <- pca_res$x[, 2]
high_cons_pca$PC3 <- pca_res$x[, 3]

# Plot
library(ggplot2)

high_cons_pca$site_class <- factor(high_cons_pca$site_class,
                                   levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))

p1 <- ggplot(high_cons_pca, aes(x = PC1
                                , y = PC2, color = site_class)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_color_manual(values = c("not_sig" = "gray70",
                                "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue",
                                "sig_both" = "darkred")) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  labs(title = "PCA of Structural Features",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")) +
  theme_classic() +
  theme(legend.position = "top")

p1

# Loadings - what drives each PC?
loadings <- as.data.frame(pca_res$rotation[, 1:3])
loadings$feature <- rownames(loadings)
print(loadings[order(abs(loadings$PC1), decreasing = TRUE), ])

# ---- 9. SAVE ----
saveRDS(list(
  struct_var = struct_var,
  merged = merged,
  high_cons = high_cons,
  enrichment = results
), "results/structural_enrichment.rds")

message("\nSaved results to results/structural_enrichment.rds")
message("Saved plots to results/structure_enrichment_plot.pdf/png")