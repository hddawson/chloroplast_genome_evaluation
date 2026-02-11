library(arrow)
library(data.table)
library(ggplot2)
library(stringr)
library(Biostrings)
library(viridis)
library(patchwork)



# ---- 1. LOAD GWAS ----
model_files <- list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

sites_df <- rbindlist(lapply(model_files, function(f) {
  rbindlist(lapply(readRDS(f), function(m) {
    data.table(Gene = m$Gene, Position = m$Aligned_Position,
               P_aa_only = m$P_aa_only, P_aa_with_pcs = m$P_aa_with_pcs, N = m$N)
  }))
}))
message("GWAS sites: ", nrow(sites_df))

# ---- 2. LOAD NETSURF ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

# ---- 3. MAP NETSURF TO ALIGNMENT ----
aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

netsurf_mapped <- rbindlist(lapply(aln_files, function(aln_file) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  netsurf_gene <- netsurf_df[Gene == gene]
  if (nrow(netsurf_gene) == 0) return(NULL)
  
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*", "", names(aln))
  common_ids <- intersect(unique(netsurf_gene$ID), names(aln))
  if (length(common_ids) == 0) return(NULL)
  
  rbindlist(lapply(common_ids, function(acc) {
    chars <- strsplit(as.character(aln[[acc]]), "")[[1]]
    pos_map <- cumsum(chars != "-")
    pos_map[chars == "-"] <- NA
    pos_lookup <- setNames(seq_along(chars), pos_map)
    
    ns_acc <- netsurf_gene[ID == acc]
    ns_acc[, Position := as.integer(pos_lookup[as.character(n)])]
    ns_acc[!is.na(Position)]
  }))
}), fill = TRUE)

stopifnot(nrow(netsurf_mapped) > 0)

# ---- 4. STRUCTURAL CONSENSUS BY POSITION ----
struct_var <- netsurf_mapped[, .(
  n_seqs = .N,
  consensus_q3 = names(which.max(table(q3))),
  consensus_freq_q3 = max(table(q3)) / .N,
  mean_rsa = mean(rsa),
  mean_asa = mean(asa),
  mean_disorder = mean(disorder)
), by = .(Gene, Position)]

# ---- 5. MERGE GWAS + STRUCTURE ----
high_cons <- merge(sites_df, struct_var, by = c("Gene", "Position"), all.x = FALSE)
high_cons <- high_cons[consensus_freq_q3 >= 0.8]
stopifnot(nrow(high_cons) > 0)
message("high_cons: ", nrow(high_cons), " sites")

# ---- 6. POSITION MAPPING ----
pos_maps <- netsurf_mapped[, .(Aligned_Position = Position, Ungapped_Position = n), by = .(Gene, ID)]
pos_consensus <- pos_maps[, .(Residue_Index = as.integer(median(Ungapped_Position))), by = .(Gene, Aligned_Position)]

high_cons <- merge(high_cons, pos_consensus, 
                   by.x = c("Gene", "Position"), 
                   by.y = c("Gene", "Aligned_Position"), 
                   all.x = TRUE)

# Prep features
high_cons[, log_p_aa := -log10(P_aa_only + 1e-300)]
high_cons[, log_p_controlled := -log10(P_aa_with_pcs + 1e-300)]

# ---- 7. LOAD PCA SCORES ----
pca_files <- list.files("data/tmp/", pattern = "_pca_scores\\.parquet$", full.names = TRUE)
stopifnot(length(pca_files) > 0)
message("Found ", length(pca_files), " PCA score files")

pca_data <- rbindlist(lapply(pca_files, function(f) {
  message("Loading: ", basename(f))
  as.data.table(read_parquet(f))
}), fill = TRUE)

message("Total PCA rows: ", nrow(pca_data))
message("Genes in PCA: ", length(unique(pca_data$Gene)))

# ---- 8. MERGE PCA WITH FEATURES ----
merged <- merge(pca_data, high_cons, 
                by = c("Gene", "Residue_Index"), 
                all.x = FALSE)

message("Merged dataset: ", nrow(merged), " rows")
message("Genes in merged: ", length(unique(merged$Gene)))
stopifnot(nrow(merged) > 0)

# Add GWAS classification
merged[, gwas_class := "not_sig"]
merged[P_aa_only < quantile(P_aa_only, 0.25, na.rm = TRUE), gwas_class := "sig_no_control"]
merged[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05, na.rm = TRUE), gwas_class := "sig_with_control"]
merged[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05, na.rm = TRUE) & 
         P_aa_only < quantile(P_aa_only, 0.25, na.rm = TRUE), gwas_class := "sig_both"]
merged[, gwas_class := factor(gwas_class, levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

# ---- 9. PLOT FUNCTIONS ----
plot_pca_hex <- function(data, color_var, title, bins = 80, is_numeric = TRUE) {
  if (is_numeric) {
    ggplot(data, aes(x = PC1, y = PC2, z = .data[[color_var]])) +
      stat_summary_hex(fun = mean, bins = bins) +
      scale_fill_viridis_c(name = color_var) +
      labs(title = title, x = "PC1", y = "PC2") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")
  } else {
    ggplot(data, aes(x = PC1, y = PC2)) +
      geom_hex(bins = bins) +
      scale_fill_viridis_c(name = "count", trans = "log10") +
      facet_wrap(as.formula(paste("~", color_var))) +
      labs(title = title, x = "PC1", y = "PC2") +
      theme_minimal(base_size = 10)
  }
}

# ---- 10. GENE LOCALIZATION PLOTS ----
message("\nCreating gene localization plots...")

# Overall density
p_density <- ggplot(merged, aes(x = PC1, y = PC2)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  labs(title = "PCA: All Genes Combined", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 14)
ggsave("results/pca_overall_density.png", p_density, width = 8, height = 6, dpi = 150)

# By gene (faceted)
genes_sorted <- merged[, .N, by = Gene][order(-N)][, Gene]
n_plot <- min(16, length(genes_sorted))

p_genes <- ggplot(merged[Gene %in% genes_sorted[1:n_plot]], aes(x = PC1, y = PC2)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  facet_wrap(~ Gene, ncol = 4) +
  labs(title = paste0("PCA by Gene (Top ", n_plot, " by N)"), x = "PC1", y = "PC2") +
  theme_minimal(base_size = 10)
ggsave("results/pca_by_gene.png", p_genes, width = 14, height = 10, dpi = 150)

# ---- 11. FEATURE COLORING PLOTS ----
message("Creating feature-colored plots...")

# GWAS p-values
p_pval_aa <- plot_pca_hex(merged, "log_p_aa", "PCA colored by -log10(P_aa_only)")
ggsave("results/pca_colored_log_p_aa.png", p_pval_aa, width = 8, height = 6, dpi = 150)

p_pval_ctrl <- plot_pca_hex(merged, "log_p_controlled", "PCA colored by -log10(P_aa_with_pcs)")
ggsave("results/pca_colored_log_p_controlled.png", p_pval_ctrl, width = 8, height = 6, dpi = 150)

# Structural features
p_rsa <- plot_pca_hex(merged, "mean_rsa", "PCA colored by Mean RSA")
ggsave("results/pca_colored_rsa.png", p_rsa, width = 8, height = 6, dpi = 150)

p_asa <- plot_pca_hex(merged, "mean_asa", "PCA colored by Mean ASA")
ggsave("results/pca_colored_asa.png", p_asa, width = 8, height = 6, dpi = 150)

p_disorder <- plot_pca_hex(merged, "mean_disorder", "PCA colored by Mean Disorder")
ggsave("results/pca_colored_disorder.png", p_disorder, width = 8, height = 6, dpi = 150)

# Q3 structure
p_q3 <- plot_pca_hex(merged, "consensus_q3", "PCA by Q3 Secondary Structure", is_numeric = FALSE)
ggsave("results/pca_by_q3.png", p_q3, width = 10, height = 8, dpi = 150)

# GWAS class
p_gwas <- plot_pca_hex(merged, "gwas_class", "PCA by GWAS Significance Class", is_numeric = FALSE)
ggsave("results/pca_by_gwas_class.png", p_gwas, width = 10, height = 8, dpi = 150)

# ---- 12. PC3 vs PC4 ----
p_pc34_density <- ggplot(merged, aes(x = PC3, y = PC4)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  labs(title = "PC3 vs PC4: Overall Density", x = "PC3", y = "PC4") +
  theme_minimal(base_size = 14)
ggsave("results/pca_pc3_pc4_density.png", p_pc34_density, width = 8, height = 6, dpi = 150)

p_pc34_pval <- ggplot(merged, aes(x = PC3, y = PC4, z = log_p_aa)) +
  stat_summary_hex(fun = mean, bins = 80) +
  scale_fill_viridis_c(name = "log_p_aa") +
  labs(title = "PC3 vs PC4 colored by -log10(P)", x = "PC3", y = "PC4") +
  theme_minimal(base_size = 12)
ggsave("results/pca_pc3_pc4_pval.png", p_pc34_pval, width = 8, height = 6, dpi = 150)

# ---- 13. SUMMARY STATS ----
message("\n=== PCA SUMMARY ===")
message("Total residues: ", nrow(merged))
message("Genes: ", length(unique(merged$Gene)))
message("PC1 range: [", round(min(merged$PC1), 2), ", ", round(max(merged$PC1), 2), "]")
message("PC2 range: [", round(min(merged$PC2), 2), ", ", round(max(merged$PC2), 2), "]")

gene_counts <- merged[, .N, by = Gene][order(-N)]
message("\nTop 5 genes by residue count:")
print(head(gene_counts, 5))

message("\nAll plots saved to results/")
message("Done.")