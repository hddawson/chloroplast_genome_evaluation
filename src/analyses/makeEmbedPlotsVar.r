library(arrow)
library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(Biostrings)
library(ape)

# ---- 1. LOAD PCA DATA ----
message("Loading PCA projected embeddings...")
pca_data <- as.data.table(read_parquet("results/embeddings_full_pca.parquet"))
message("Total PCA rows: ", nrow(pca_data))

# ---- 2. ADD TAXONOMY (ORDER) ----
message("\nAdding taxonomy information...")
processed_data <- as.data.table(read_parquet("data/processed_data.parquet"))

# Extract Order for each ID
tax_lookup <- unique(processed_data[, .(ID, Order)])
message("Taxonomy lookup: ", nrow(tax_lookup), " IDs")

pca_data <- merge(pca_data, tax_lookup, by = "ID", all.x = TRUE)
message("PCA data with Order: ", sum(!is.na(pca_data$Order)), " / ", nrow(pca_data), " rows")

# ---- 3. COMPUTE NUCLEOTIDE DIVERSITY (PI) PER POSITION ----
message("\nComputing nucleotide diversity per position...")

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)
stopifnot(length(aln_files) > 0)

compute_pi <- function(bases) {
  # Remove gaps
  bases <- bases[bases != "-"]
  if (length(bases) < 2) return(NA)
  
  # Count alleles
  freqs <- table(bases) / length(bases)
  
  # Pi = expected heterozygosity = 1 - sum(p^2)
  pi <- 1 - sum(freqs^2)
  return(pi)
}

diversity_data <- rbindlist(lapply(aln_files, function(aln_file) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  message("  Processing: ", gene)
  
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*", "", names(aln))
  
  # Convert to matrix
  aln_mat <- as.matrix(aln)
  aln_len <- ncol(aln_mat)
  
  # For each aligned position, compute diversity
  pos_diversity <- sapply(1:aln_len, function(pos) {
    compute_pi(aln_mat[, pos])
  })
  
  # Map aligned positions to ungapped positions
  # Use first sequence as reference for mapping
  ref_seq <- strsplit(as.character(aln[[1]]), "")[[1]]
  ungapped_pos <- cumsum(ref_seq != "-")
  ungapped_pos[ref_seq == "-"] <- NA
  
  # Create data.table
  dt <- data.table(
    Gene = gene,
    Aligned_Position = 1:aln_len,
    Ungapped_Position = ungapped_pos,
    Pi = pos_diversity,
    n_seqs = nrow(aln_mat)
  )
  
  # Remove gap positions and aggregate by ungapped position
  dt_ungapped <- dt[!is.na(Ungapped_Position)]
  dt_ungapped[, .(
    Pi = mean(Pi, na.rm = TRUE),
    n_seqs = first(n_seqs)
  ), by = .(Gene, Residue_Index = Ungapped_Position)]
}))

message("Computed diversity for ", nrow(diversity_data), " positions")

# ---- 4. MERGE DIVERSITY WITH PCA ----
pca_data <- merge(pca_data, diversity_data, by = c("Gene", "Residue_Index"), all.x = TRUE)
message("Merged diversity: ", sum(!is.na(pca_data$Pi)), " / ", nrow(pca_data), " rows")

# ---- 5. COMPUTE POSITION VARIABILITY ACROSS TAXA ----
message("\nComputing PC variability per position...")
pos_variability <- pca_data[, .(
  n_taxa = .N,
  PC1_sd = sd(PC1, na.rm = TRUE),
  PC1_mean = mean(PC1, na.rm = TRUE),
  PC1_range = max(PC1, na.rm = TRUE) - min(PC1, na.rm = TRUE),
  PC2_sd = sd(PC2, na.rm = TRUE)
), by = .(Gene, Residue_Index)]

pca_data <- merge(pca_data, pos_variability, by = c("Gene", "Residue_Index"))

# ---- 6. PLOT: TAXONOMY ----
message("\nCreating taxonomy plots...")

# Get top orders by count
top_orders <- pca_data[!is.na(Order), .N, by = Order][order(-N)][1:12, Order]

# Overall by Order
p_order_overall <- ggplot(pca_data[Order %in% top_orders], aes(x = PC1, y = PC2)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  facet_wrap(~ Order, ncol = 4) +
  labs(title = "PC1 vs PC2 by Order (Top 12)") +
  theme_minimal(base_size = 10)
ggsave("results/pca_by_order_facet.png", p_order_overall, width = 14, height = 10, dpi = 150)

# Colored points (sampled)
set.seed(42)
sample_size <- min(100000, nrow(pca_data[Order %in% top_orders]))
pca_sample <- pca_data[Order %in% top_orders][sample(.N, sample_size)]

p_order_colored <- ggplot(pca_sample, aes(x = PC1, y = PC2, color = Order)) +
  geom_point(alpha = 0.2, size = 0.3) +
  scale_color_brewer(palette = "Paired", name = "Order") +
  labs(title = paste0("PC1 vs PC2 colored by Order (", sample_size, " sampled points)")) +
  theme_minimal(base_size = 12)
ggsave("results/pca_by_order_colored.png", p_order_colored, width = 10, height = 7, dpi = 150)

# Density by major orders
major_orders <- pca_data[!is.na(Order), .N, by = Order][order(-N)][1:4, Order]

for (ord in major_orders) {
  p_dens <- ggplot(pca_data[Order == ord], aes(x = PC1, y = PC2)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(name = "count", trans = "log10") +
    labs(title = paste0("PC1 vs PC2: ", ord, " only")) +
    theme_minimal(base_size = 12)
  
  fname <- paste0("results/pca_order_", gsub("[^A-Za-z0-9]", "_", ord), ".png")
  ggsave(fname, p_dens, width = 8, height = 6, dpi = 150)
}

# ---- 7. PLOT: NUCLEOTIDE DIVERSITY (PI) ----
message("\nCreating diversity plots...")

p_pi <- ggplot(pca_data[!is.na(Pi)], aes(x = PC1, y = PC2, z = Pi)) +
  stat_summary_hex(fun = mean, bins = 80) +
  scale_fill_viridis_c(name = "Mean Pi") +
  labs(title = "PC1 vs PC2 colored by nucleotide diversity (Pi)") +
  theme_minimal(base_size = 12)
ggsave("results/pca_colored_by_pi.png", p_pi, width = 9, height = 7, dpi = 150)

# Pi distribution across PC1
p_pi_dist <- ggplot(pca_data[!is.na(Pi)], aes(x = PC1, y = Pi)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  geom_smooth(color = "red", se = TRUE) +
  labs(title = "Nucleotide diversity (Pi) vs PC1", 
       x = "PC1", y = "Nucleotide Diversity (Pi)") +
  theme_minimal(base_size = 12)
ggsave("results/pi_vs_pc1.png", p_pi_dist, width = 8, height = 6, dpi = 150)

# ---- 8. PLOT: POSITION VARIABILITY ----
message("\nCreating position variability plots...")

p_var <- ggplot(pca_data, aes(x = PC1, y = PC2, z = PC1_sd)) +
  stat_summary_hex(fun = mean, bins = 80) +
  scale_fill_viridis_c(name = "PC1 SD\n(variability)") +
  labs(title = "PC1 vs PC2 colored by position variability across taxa") +
  theme_minimal(base_size = 12)
ggsave("results/pca_colored_by_variability.png", p_var, width = 9, height = 7, dpi = 150)

# Variability vs PC1
p_var_dist <- ggplot(unique(pca_data[, .(Gene, Residue_Index, PC1_mean, PC1_sd)]), 
                     aes(x = PC1_mean, y = PC1_sd)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(name = "count", trans = "log10") +
  geom_smooth(color = "red", se = TRUE) +
  labs(title = "Position variability (PC1 SD) vs mean PC1", 
       x = "Mean PC1", y = "PC1 SD across taxa") +
  theme_minimal(base_size = 12)
ggsave("results/variability_vs_pc1.png", p_var_dist, width = 8, height = 6, dpi = 150)

# ---- 9. CORRELATION ANALYSIS ----
message("\nComputing correlations...")

cor_data <- unique(pca_data[, .(Gene, Residue_Index, PC1_mean, PC1_sd, Pi)])
cor_data <- cor_data[complete.cases(cor_data)]

cor_pc1_pi <- cor(cor_data$PC1_mean, cor_data$Pi)
cor_pc1_var <- cor(cor_data$PC1_mean, cor_data$PC1_sd)
cor_pi_var <- cor(cor_data$Pi, cor_data$PC1_sd)

message("\n=== CORRELATIONS ===")
message("PC1 vs Pi: ", round(cor_pc1_pi, 3))
message("PC1 vs Variability (PC1_sd): ", round(cor_pc1_var, 3))
message("Pi vs Variability: ", round(cor_pi_var, 3))

# ---- 10. CLUSTER IDENTIFICATION ----
message("\nDefining clusters...")

pca_data[, cluster := ifelse(PC1 > 0, "Right_cluster", "Left_fan")]

cluster_summary <- pca_data[, .(
  n = .N,
  mean_Pi = mean(Pi, na.rm = TRUE),
  mean_PC1_sd = mean(PC1_sd, na.rm = TRUE),
  n_orders = length(unique(Order[!is.na(Order)]))
), by = .(Gene, cluster)]

message("\nCluster summary:")
print(cluster_summary[order(Gene, cluster)])

# Plot clusters
p_clusters <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.05, size = 0.3) +
  scale_color_manual(values = c("Right_cluster" = "red", "Left_fan" = "blue")) +
  labs(title = "Right cluster vs Left fan") +
  theme_minimal()
ggsave("results/pca_clusters_identified.png", p_clusters, width = 8, height = 6, dpi = 150)

# ---- 11. GENE-SPECIFIC: MatK ----
message("\nAnalyzing MatK specifically...")

matk_data <- pca_data[Gene == "MatK"]

p_matk_pi <- ggplot(matk_data[!is.na(Pi)], aes(x = PC1, y = PC2, color = Pi)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_viridis_c(name = "Pi") +
  labs(title = "MatK: PC1 vs PC2 colored by Pi") +
  theme_minimal()
ggsave("results/matk_colored_by_pi.png", p_matk_pi, width = 8, height = 6, dpi = 150)

p_matk_order <- ggplot(matk_data[!is.na(Order)], aes(x = PC1, y = PC2, color = Order)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = "MatK: PC1 vs PC2 colored by Order") +
  theme_minimal()
ggsave("results/matk_colored_by_order.png", p_matk_order, width = 9, height = 6, dpi = 150)

message("\nDone! Check results/ for diagnostic plots.")
message("\nKey findings:")
message("  - Correlation PC1 vs Pi: ", round(cor_pc1_pi, 3))
message("  - Correlation PC1 vs Variability: ", round(cor_pc1_var, 3))