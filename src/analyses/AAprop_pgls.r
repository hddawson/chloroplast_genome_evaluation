library(arrow)
library(data.table)
library(ggplot2)
library(stringr)
library(Biostrings)
library(pheatmap)
library(viridis)

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
  consensus_q8 = names(which.max(table(q8))),
  consensus_freq_q8 = max(table(q8)) / .N,
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

feature_vars <- c("log_p_aa", "log_p_controlled", "mean_rsa", "mean_asa", "mean_disorder", "consensus_freq_q3")

# ---- 7. LOAD ALL EMBEDDINGS ONCE ----
embed_files <- list.files("data/embeddings/", pattern = "_residue_embeddings\\.parquet$", full.names = TRUE)
stopifnot(length(embed_files) > 0)
message("Found ", length(embed_files), " embedding files")

message("Loading all embeddings into memory...")
embeds_list <- lapply(embed_files, function(f) {
  message("  Loading: ", basename(f))
  as.data.table(read_parquet(f))
})
embeds_data <- rbindlist(embeds_list, fill = TRUE)
message("Total embeddings: ", nrow(embeds_data), " rows")

embed_cols <- grep("^embedding_", names(embeds_data), value = TRUE)
n_dims <- length(embed_cols)
message("Embedding dimensions: ", n_dims)
stopifnot(n_dims > 0)

# ---- 8. MERGE ONCE WITH FEATURES ----
merged <- merge(high_cons, 
                embeds_data[, c("ID", "Gene", "Residue_Index", embed_cols), with = FALSE],
                by = c("Gene", "Residue_Index"), 
                all.x = FALSE)
message("Merged dataset: ", nrow(merged), " rows")
stopifnot(nrow(merged) > 0)

# ---- 9. DIMENSION-BY-DIMENSION CORRELATION ----
cor_results <- data.table(dimension = integer(), variable = character(), 
                         correlation = numeric(), n_obs = integer())

for (i in seq_along(embed_cols)) {
  dim_name <- embed_cols[i]
  if (i %% 50 == 0) message("Processing dimension ", i, "/", n_dims)
  
  for (fvar in feature_vars) {
    complete_idx <- complete.cases(merged[[dim_name]], merged[[fvar]])
    if (sum(complete_idx) > 10) {
      cor_val <- cor(merged[[dim_name]][complete_idx], merged[[fvar]][complete_idx])
      cor_results <- rbind(cor_results, data.table(
        dimension = i,
        variable = fvar,
        correlation = cor_val,
        n_obs = sum(complete_idx)
      ))
    }
  }
}

# ---- 10. SAVE RESULTS ----
saveRDS(cor_results, "results/embedding_dimension_correlations.rds")
fwrite(cor_results, "results/embedding_dimension_correlations.csv")

message("\nSaved correlation results: ", nrow(cor_results), " rows")

# ---- 11. VISUALIZE TOP CORRELATIONS ----
cor_results[, abs_cor := abs(correlation)]
top_cors <- cor_results[order(-abs_cor)][1:50]

p_top <- ggplot(top_cors, aes(x = dimension, y = correlation, color = variable)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Top 50 Embedding-Feature Correlations", 
       x = "Embedding Dimension", y = "Correlation") +
  theme_minimal()
ggsave("results/top_embedding_correlations.png", p_top, width = 10, height = 6, dpi = 150)

# Heatmap of all correlations
cor_wide <- dcast(cor_results, dimension ~ variable, value.var = "correlation")
cor_matrix <- as.matrix(cor_wide[, -1])
rownames(cor_matrix) <- cor_wide$dimension

png("results/embedding_feature_correlations_heatmap.png", width = 8, height = 12, units = "in", res = 150)
pheatmap(cor_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "Embedding Dimensions vs Features")
dev.off()

message("Done.")