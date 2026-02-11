library(arrow)
library(data.table)
library(ggplot2)
library(stringr)
library(Biostrings)
library(viridis)
library(patchwork)

# ---- 1. LOAD PCA PROJECTED DATA ----
message("Loading PCA projected embeddings...")
pca_file <- "results/embeddings_full_pca.parquet"
stopifnot(file.exists(pca_file))

pca_data <- as.data.table(read_parquet(pca_file))
message("Total PCA rows: ", nrow(pca_data))
message("Genes in PCA: ", length(unique(pca_data$Gene)))

# ---- 2. LOAD GWAS ----
model_files <- list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

sites_df <- rbindlist(lapply(model_files, function(f) {
  rbindlist(lapply(readRDS(f), function(m) {
    data.table(Gene = m$Gene, Position = m$Aligned_Position,
               P_aa_only = m$P_aa_only, P_aa_with_pcs = m$P_aa_with_pcs, N = m$N)
  }))
}))
message("GWAS sites: ", nrow(sites_df))

# ---- 3. LOAD NETSURF ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

# ---- 4. MAP NETSURF TO ALIGNMENT ----
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

# ---- 5. STRUCTURAL CONSENSUS BY POSITION ----
struct_var <- netsurf_mapped[, .(
  n_seqs = .N,
  consensus_q3 = names(which.max(table(q3))),
  consensus_freq_q3 = max(table(q3)) / .N,
  mean_rsa = mean(rsa),
  mean_asa = mean(asa),
  mean_disorder = mean(disorder)
), by = .(Gene, Position)]

# ---- 6. MERGE GWAS + STRUCTURE ----
high_cons <- merge(sites_df, struct_var, by = c("Gene", "Position"), all.x = FALSE)
high_cons <- high_cons[consensus_freq_q3 >= 0.8]
stopifnot(nrow(high_cons) > 0)
message("high_cons: ", nrow(high_cons), " sites")

# ---- 7. POSITION MAPPING ----
pos_maps <- netsurf_mapped[, .(Aligned_Position = Position, Ungapped_Position = n), by = .(Gene, ID)]
pos_consensus <- pos_maps[, .(Residue_Index = as.integer(median(Ungapped_Position))), by = .(Gene, Aligned_Position)]

high_cons <- merge(high_cons, pos_consensus, 
                   by.x = c("Gene", "Position"), 
                   by.y = c("Gene", "Aligned_Position"), 
                   all.x = TRUE)

# Prep features
high_cons[, log_p_aa := -log10(P_aa_only + 1e-300)]
high_cons[, log_p_controlled := -log10(P_aa_with_pcs + 1e-300)]

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

# ---- 9. EXTRACT AMINO ACIDS FROM ALIGNMENTS ----
message("\nExtracting amino acids from alignments...")

aa_data <- rbindlist(lapply(aln_files, function(aln_file) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  message("  Processing alignments for: ", gene)
  
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*", "", names(aln))
  
  rbindlist(lapply(names(aln), function(acc) {
    seq_char <- strsplit(as.character(aln[[acc]]), "")[[1]]
    pos_map <- cumsum(seq_char != "-")
    
    # Create mapping of ungapped position to amino acid
    ungapped_pos <- pos_map[seq_char != "-"]
    ungapped_aa <- seq_char[seq_char != "-"]
    
    data.table(
      ID = acc,
      Gene = gene,
      Residue_Index = ungapped_pos,
      Amino_Acid = ungapped_aa
    )
  }))
}))

message("Extracted ", nrow(aa_data), " amino acid annotations")

# Add amino acids to merged data
merged <- merge(merged, aa_data, 
                by = c("Gene", "Residue_Index", "ID"),
                all.x = TRUE)
message("Amino acids available: ", sum(!is.na(merged$Amino_Acid)))

# ---- 10. PLOT FUNCTIONS ----
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

# ---- 11. GENE LOCALIZATION PLOTS ----
message("\nCreating gene localization plots...")

# Overall density
p_density <- ggplot(merged, aes(x = PC1, y = PC2)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  labs(title = "PCA: All Genes Combined", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 14)
ggsave("results/pca_overall_density.png", p_density, width = 8, height = 6, dpi = 150)

# By gene (ALL GENES)
genes_sorted <- merged[, .N, by = Gene][order(-N)][, Gene]
n_genes <- length(genes_sorted)
message("Plotting all ", n_genes, " genes...")

# Calculate grid dimensions
ncols <- ceiling(sqrt(n_genes))
nrows <- ceiling(n_genes / ncols)

p_genes_all <- ggplot(merged, aes(x = PC1, y = PC2)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  facet_wrap(~ Gene, ncol = ncols) +
  labs(title = paste0("PCA by Gene (All ", n_genes, " genes)"), x = "PC1", y = "PC2") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(size = 6))

# Save with appropriate dimensions
plot_width <- max(16, ncols * 2.5)
plot_height <- max(12, nrows * 2.5)
ggsave("results/pca_by_gene_all.png", p_genes_all, 
       width = plot_width, height = plot_height, dpi = 150, limitsize = FALSE)

# ---- 12. FEATURE COLORING PLOTS ----
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

# ---- 13. Q3 STRUCTURE PLOTS ----
message("Creating Q3 structure plots...")

# Q3 faceted
p_q3_facet <- ggplot(merged, aes(x = PC1, y = PC2)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  facet_wrap(~ consensus_q3) +
  labs(title = "PCA by Q3 Secondary Structure", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12)
ggsave("results/pca_by_q3_facet.png", p_q3_facet, width = 10, height = 8, dpi = 150)

# Q3 colored points (sample for visibility)
set.seed(42)
sample_size <- min(50000, nrow(merged))
merged_sample <- merged[sample(.N, sample_size)]

p_q3_points <- ggplot(merged_sample, aes(x = PC1, y = PC2, color = consensus_q3)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Q3") +
  labs(title = paste0("PCA colored by Q3 (", sample_size, " sampled points)"), 
       x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12)
ggsave("results/pca_by_q3_points.png", p_q3_points, width = 8, height = 6, dpi = 150)

# ---- 14. AMINO ACID PLOTS ----
message("Creating amino acid plots...")

if (sum(!is.na(merged$Amino_Acid)) > 0) {
  # Amino acid faceted
  p_aa_facet <- ggplot(merged[!is.na(Amino_Acid)], aes(x = PC1, y = PC2)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(name = "count", trans = "log10") +
    facet_wrap(~ Amino_Acid) +
    labs(title = "PCA by Amino Acid", x = "PC1", y = "PC2") +
    theme_minimal(base_size = 10)
  ggsave("results/pca_by_amino_acid.png", p_aa_facet, width = 14, height = 12, dpi = 150)
  
  # Amino acid colored points (sample)
  merged_aa_sample <- merged[!is.na(Amino_Acid)][sample(.N, min(50000, .N))]
  
  p_aa_points <- ggplot(merged_aa_sample, aes(x = PC1, y = PC2, color = Amino_Acid)) +
    geom_point(alpha = 0.3, size = 0.5) +
    scale_color_discrete(name = "AA") +
    labs(title = paste0("PCA colored by Amino Acid (", nrow(merged_aa_sample), " sampled points)"), 
         x = "PC1", y = "PC2") +
    theme_minimal(base_size = 12)
  ggsave("results/pca_by_amino_acid_points.png", p_aa_points, width = 9, height = 6, dpi = 150)
  
  # Amino acid counts
  aa_counts <- merged[!is.na(Amino_Acid), .N, by = Amino_Acid][order(-N)]
  message("\nAmino acid distribution:")
  print(aa_counts)
} else {
  message("No amino acid data available for plotting")
}

# GWAS class
p_gwas <- plot_pca_hex(merged, "gwas_class", "PCA by GWAS Significance Class", is_numeric = FALSE)
ggsave("results/pca_by_gwas_class.png", p_gwas, width = 10, height = 8, dpi = 150)

# ---- 15. PC3 vs PC4 ----
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

p_pc34_pval <- ggplot(merged, aes(x = PC3, y = PC4, z = log_p_controlled)) +
  stat_summary_hex(fun = mean, bins = 80) +
  scale_fill_viridis_c(name = "log_p_controlled") +
  labs(title = "PC3 vs PC4 colored by -log10(P)", x = "PC3", y = "PC4") +
  theme_minimal(base_size = 12)
ggsave("results/pca_pc3_pc4_pval_controlled.png", p_pc34_pval, width = 8, height = 6, dpi = 150)

# ---- 16. SUMMARY STATS ----
message("\n=== PCA SUMMARY ===")
message("Total residues: ", nrow(merged))
message("Genes: ", length(unique(merged$Gene)))
message("PC1 range: [", round(min(merged$PC1), 2), ", ", round(max(merged$PC1), 2), "]")
message("PC2 range: [", round(min(merged$PC2), 2), ", ", round(max(merged$PC2), 2), "]")

gene_counts <- merged[, .N, by = Gene][order(-N)]
message("\nGene counts:")
print(gene_counts)

message("\nAll plots saved to results/")
message("Done.")