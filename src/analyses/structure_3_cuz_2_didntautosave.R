# structural_enrichment_analysis.R
# Analyze secondary structure enrichment in GWAS hits

library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD NETSURF DATA ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))

# Parse IDs: ">MK637829.1_Gene_atpF_Taxonomy_..."
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

stopifnot(all(!is.na(netsurf_df$ID)))
stopifnot(all(!is.na(netsurf_df$Gene)))
message("Loaded ", nrow(netsurf_df), " netsurf rows, ", 
        length(unique(netsurf_df$Gene)), " genes, ",
        length(unique(netsurf_df$ID)), " sequences")

# ---- 2. LOAD GWAS RESULTS ----
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
stopifnot(nrow(sites_df) > 0)
message("Loaded ", nrow(sites_df), " GWAS sites across ", length(unique(sites_df$Gene)), " genes")

# ---- 3. MAP NETSURF POSITIONS TO ALIGNMENT ----
aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

netsurf_mapped_list <- list()

for (aln_file in aln_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  
  netsurf_gene <- netsurf_df[Gene == gene]
  if (nrow(netsurf_gene) == 0) next
  
  aln <- readAAStringSet(aln_file)
  aln_ids <- sub("\\|.*", "", names(aln))
  names(aln) <- aln_ids
  
  common_ids <- intersect(unique(netsurf_gene$ID), aln_ids)
  if (length(common_ids) == 0) {
    message("Gene ", gene, ": no overlap")
    next
  }
  
  for (acc in common_ids) {
    chars <- strsplit(as.character(aln[[acc]]), "")[[1]]
    
    # ungapped position -> alignment position
    ungapped <- 0L
    pos_map <- integer(sum(chars != "-"))
    for (i in seq_along(chars)) {
      if (chars[i] != "-") {
        ungapped <- ungapped + 1L
        pos_map[ungapped] <- i
      }
    }
    
    ns_acc <- netsurf_gene[ID == acc]
    ns_acc[, Position := pos_map[n]]
    
    # Filter positions that exceed alignment (truncated seqs)
    ns_acc <- ns_acc[!is.na(Position)]
    if (nrow(ns_acc) > 0) {
      netsurf_mapped_list[[paste(gene, acc, sep = "_")]] <- ns_acc
    }
  }
  message("Gene ", gene, ": mapped ", length(common_ids), " sequences")
}

netsurf_mapped <- rbindlist(netsurf_mapped_list, fill = TRUE)
stopifnot(nrow(netsurf_mapped) > 0)
message("Total mapped: ", nrow(netsurf_mapped), " residue predictions")

# ---- 4. CALCULATE STRUCTURAL CONSENSUS BY POSITION ----
pred_confidence <- netsurf_mapped[, .(
  mean.p.q3.H = mean(p.q3_H.),
  mean.p.q3.C = mean(p.q3_C.),
  mean.p.q3.E = mean(p.q3_E.)
), by = .(Gene, Position)]

summary(pred_confidence)

hist(pred_confidence$mean.p.q3.H)
plot(pred_confidence$mean.p.q3.H, pred_confidence$mean.p.q3.C)
plot(pred_confidence$mean.p.q3.H, pred_confidence$mean.p.q3.E)
plot(pred_confidence$mean.p.q3.C, pred_confidence$mean.p.q3.E)

struct_var <- netsurf_mapped[, .(
  n_seqs = .N,
  prop_H = mean(q3 == "H"),
  prop_E = mean(q3 == "E"),
  prop_C = mean(q3 == "C"),
  entropy = {
    props <- c(mean(q3 == "H"), mean(q3 == "E"), mean(q3 == "C"))
    props <- props[props > 0]
    if (length(props) > 1) -sum(props * log2(props)) else 0
  },
  consensus = names(which.max(table(q3))),
  consensus_freq = max(table(q3)) / .N
), by = .(Gene, Position)]

struct_var[, conservation := cut(consensus_freq,
                                 breaks = c(0, 0.5, 0.95, 1),
                                 labels = c("variable", "moderate", "high"),
                                 include.lowest = TRUE)]

message("\n=== Structural Conservation ===")
print(table(struct_var$conservation))

# ---- 5. MERGE WITH GWAS ----
merged <- merge(sites_df, struct_var, by = c("Gene", "Position"), all.x = TRUE)
message("Merged: ", sum(!is.na(merged$consensus)), " of ", nrow(merged), " GWAS sites have structure data")

# Filter to high-confidence structural calls
high_cons <- merged[!is.na(consensus_freq) & consensus_freq >= 0.8]
stopifnot(nrow(high_cons) > 0)

# ---- 6. CLASSIFY SITES ----
high_cons[, site_class := "not_sig"]
high_cons[P_aa_only < quantile(P_aa_only, 0.25), site_class := "sig_no_control"]
high_cons[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05), site_class := "sig_with_control"]
high_cons[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05) & 
            P_aa_only < quantile(P_aa_only, 0.25), site_class := "sig_both"]

message("\n=== Site Classification ===")
print(table(high_cons$site_class))

# ---- 7. ENRICHMENT TESTS ----
results <- NULL
site_classes <- c("sig_no_control", "sig_with_control", "sig_both", "not_sig")
q3_states <- c("H", "E", "C")
background <- high_cons[site_class == "not_sig"]

for (sc in site_classes) {
  test_sites <- high_cons[site_class == sc]
  if (nrow(test_sites) == 0) next
  
  for (state in q3_states) {
    obs_count <- sum(test_sites$consensus == state)
    bg_count <- sum(background$consensus == state)
    
    cont_table <- matrix(c(obs_count, nrow(test_sites) - obs_count,
                           bg_count, nrow(background) - bg_count),
                         nrow = 2, byrow = TRUE)
    
    ft <- fisher.test(cont_table)
    
    results <- rbind(results, data.frame(
      site_class = sc,
      q3_state = state,
      n_sites = nrow(test_sites),
      obs_prop = obs_count / nrow(test_sites),
      bg_prop = bg_count / nrow(background),
      fold_enrichment = (obs_count / nrow(test_sites)) / (bg_count / nrow(background)),
      p_value = ft$p.value,
      odds_ratio = as.numeric(ft$estimate)
    ))
  }
}

results$p_adj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$p_adj), ]

message("\n=== Enrichment Results ===")
print(results[, c("site_class", "q3_state", "fold_enrichment", "p_adj")])

# ---- 8. PLOT ----
results$sig_label <- ifelse(results$p_adj < 0.001, "***",
                            ifelse(results$p_adj < 0.01, "**",
                                   ifelse(results$p_adj < 0.05, "*", "")))

results$site_class <- factor(results$site_class, 
                             levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))

p <- ggplot(results, aes(x = q3_state, y = fold_enrichment, fill = site_class)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = sig_label, y = fold_enrichment + 0.05), 
            position = position_dodge(0.8), vjust = 0, size = 5) +
  scale_fill_manual(values = c("not_sig" = "gray70", 
                               "sig_no_control" = "gold", 
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred"),
                    name = "Site Class") +
  labs(x = "Secondary Structure", 
       y = "Fold Enrichment vs Background",
       title = "Structure Enrichment in GWAS Hits") +
  theme_classic() +
  theme(legend.position = "top")

p

ggsave("results/structure_enrichment_plot.pdf", p, width = 8, height = 6)
ggsave("results/structure_enrichment_plot.png", p, width = 8, height = 6, dpi = 150)

# ---- q8 ----- 

# Q8 has finer categories: G (3-10 helix), H (alpha helix), I (pi helix), 
# B (beta bridge), E (strand), S (bend), T (turn), C (coil)

struct_var_q8 <- netsurf_mapped[, .(
  n_seqs = .N,
  consensus_q8 = names(which.max(table(q8))),
  consensus_freq_q8 = max(table(q8)) / .N
), by = .(Gene, Position)]

merged_q8 <- merge(high_cons, struct_var_q8, by = c("Gene", "Position"), all.x = TRUE)
table(netsurf_mapped$q8)
# Enrichment test
results_q8 <- NULL
q8_states <- c("H", "E", "C", "T", "S", "G", "B", "I")
background_q8 <- merged_q8[site_class == "not_sig" & !is.na(consensus_q8)]

for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_sites <- merged_q8[site_class == sc & !is.na(consensus_q8)]
  if (nrow(test_sites) == 0) next
  
  for (state in q8_states) {
    obs_count <- sum(test_sites$consensus_q8 == state)
    bg_count <- sum(background_q8$consensus_q8 == state)
    
    # Skip if no observations
    if (obs_count + bg_count == 0) next
    
    cont_table <- matrix(c(obs_count, nrow(test_sites) - obs_count,
                           bg_count, nrow(background_q8) - bg_count),
                         nrow = 2, byrow = TRUE)
    
    ft <- fisher.test(cont_table)
    
    results_q8 <- rbind(results_q8, data.frame(
      site_class = sc,
      q8_state = state,
      n_sites = nrow(test_sites),
      obs_count = obs_count,
      obs_prop = obs_count / nrow(test_sites),
      bg_prop = bg_count / nrow(background_q8),
      fold_enrichment = (obs_count / nrow(test_sites)) / (bg_count / nrow(background_q8)),
      p_value = ft$p.value,
      odds_ratio = as.numeric(ft$estimate)
    ))
  }
}

results_q8$p_adj <- p.adjust(results_q8$p_value, method = "BH")
results_q8 <- results_q8[order(results_q8$p_adj), ]
print(results_q8[, c("site_class", "q8_state", "obs_count", "fold_enrichment", "p_adj")])

# Q8 annotations
q8_labels <- c(
  "H" = "H (α-helix)",
  "G" = "G (3₁₀-helix)",
  "I" = "I (π-helix)",
  "E" = "E (β-strand)",
  "B" = "B (β-bridge)",
  "T" = "T (turn)",
  "S" = "S (bend)",
  "C" = "C (coil)"
)

# Color by structural category
q8_colors <- c(
  "H" = "#E41A1C",  # red - helix
  "G" = "#FC8D62",  # light red - helix  
  "I" = "#FDAE6B",  # orange - helix
  "E" = "#377EB8",  # blue - sheet
  "B" = "#7FCDBB",  # light blue - sheet
  "T" = "#4DAF4A",  # green - loop
  "S" = "#984EA3",  # purple - loop
  "C" = "#FF7F00"   # orange - loop
)

results_q8$q8_label <- q8_labels[results_q8$q8_state]
results_q8$sig_label <- ifelse(results_q8$p_adj < 0.001, "***",
                               ifelse(results_q8$p_adj < 0.01, "**",
                                      ifelse(results_q8$p_adj < 0.05, "*", "")))

results_q8$site_class <- factor(results_q8$site_class,
                                levels = c("sig_no_control", "sig_with_control", "sig_both"))

# Filter to states with actual observations
plot_data <- results_q8[results_q8$obs_count > 0, ]

# Order q8 by structure type: helices, sheets, loops
plot_data$q8_label <- factor(plot_data$q8_label,
                             levels = c("H (α-helix)", "G (3₁₀-helix)", "I (π-helix)",
                                        "E (β-strand)", "B (β-bridge)",
                                        "T (turn)", "S (bend)", "C (coil)"))

# Add SE for fold enrichment (using delta method approximation)
# SE(fold) ≈ fold * sqrt(1/obs + 1/bg)
plot_data$bg_count <- plot_data$bg_prop * nrow(background_q8)
plot_data$se_fold <- plot_data$fold_enrichment * sqrt(1/plot_data$obs_count + 1/plot_data$bg_count)
plot_data$se_fold
# Label with n
plot_data$bar_label <- paste0("n=", plot_data$obs_count, " ", plot_data$sig_label)

p <- ggplot(plot_data, aes(x = q8_label, y = fold_enrichment, fill = q8_state)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = fold_enrichment - se_fold, 
                    ymax = fold_enrichment + se_fold),
                position = position_dodge(0.8), width = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_text(aes(label = bar_label, y = fold_enrichment + se_fold + 0.08),
            position = position_dodge(0.8), vjust = 0, size = 2.8) +
  scale_fill_manual(values = q8_colors, guide = "none") +
  facet_wrap(~site_class, ncol = 1,
             labeller = labeller(site_class = c(
               "sig_no_control" = "Significant (no pop. control)",
               "sig_with_control" = "Significant (with pop. control)",
               "sig_both" = "Significant (both)"
             ))) +
  labs(x = NULL,
       y = "Fold Enrichment vs Background",
       title = "Q8 Secondary Structure Enrichment in GWAS Hits",
       subtitle = "Helices (red) | Sheets (blue) | Loops (green/purple/orange)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"))

p

ggsave("results/q8_structure_enrichment.pdf", p, width = 10, height = 10)
ggsave("results/q8_structure_enrichment.png", p, width = 10, height = 10, dpi = 150)

ggsave("results/q8_structure_enrichment.pdf", p, width = 10, height = 10)
ggsave("results/q8_structure_enrichment.png", p, width = 10, height = 10, dpi = 150)

head(netsurf_mapped)

# ---- other stuff ---- 

# Summarize continuous features by position
cont_features <- netsurf_mapped[, .(
  mean_rsa = mean(rsa),
  mean_asa = mean(asa),
  mean_phi = mean(phi),
  mean_psi = mean(psi),
  mean_disorder = mean(disorder)
), by = .(Gene, Position)]

# Merge with high_cons
high_cons <- merge(high_cons, cont_features, by = c("Gene", "Position"), all.x = TRUE)

# Compare distributions across site classes
par(mfrow = c(2, 3))

boxplot(mean_rsa ~ site_class, data = high_cons, 
        main = "Relative Solvent Accessibility", ylab = "RSA",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_asa ~ site_class, data = high_cons,
        main = "Absolute Solvent Accessibility", ylab = "ASA (Å²)",
        col = c("gray70", "gold", "steelblue", "darkred"))

boxplot(mean_disorder ~ site_class, data = high_cons,
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
for (feat in c("mean_rsa", "mean_asa", "mean_disorder", "mean_phi", "mean_psi")) {
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
pca_cols2 <- c("mean_rsa", "mean_asa", "mean_disorder", "mean_phi", "mean_psi")

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