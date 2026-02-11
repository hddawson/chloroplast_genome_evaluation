library(arrow)
library(data.table)
library(ggplot2)
library(stringr)
library(Biostrings)


embeds_data <- read_parquet("data/all_residue_embeddings.parquet")

unique(embeds_data$Gene)
head(embeds_data$)
# ---- 1. LOAD NETSURF ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

# ---- 2. LOAD GWAS ----
model_files <- list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

sites_df <- rbindlist(lapply(model_files, function(f) {
  cat(f)
  rbindlist(lapply(readRDS(f), function(m) {
    data.table(Gene = m$Gene, Position = m$Aligned_Position,
               P_aa_only = m$P_aa_only, P_aa_with_pcs = m$P_aa_with_pcs, N = m$N)
  }))
}))

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

# ---- 5. MERGE & FILTER ----
high_cons <- merge(sites_df, struct_var, by = c("Gene", "Position"), all.x = FALSE)
high_cons <- high_cons[consensus_freq_q3 >= 0.8]
stopifnot(nrow(high_cons) > 0)

message("high_cons: ", nrow(high_cons), " sites with structure data")

# ---- 6. MERGE WITH EMBEDDINGS ----
# embeds_data uses Residue_Index (ungapped), high_cons uses Position (aligned)
# Need to map aligned positions back to ungapped for each sequence

# Get one representative sequence per gene from netsurf_mapped to build position map
pos_maps <- netsurf_mapped[, .(
  Aligned_Position = Position,
  Ungapped_Position = n
), by = .(Gene, ID)]

# For each gene/position, get the consensus ungapped position (should be consistent)
pos_consensus <- pos_maps[, .(
  Residue_Index = as.integer(median(Ungapped_Position))
), by = .(Gene, Aligned_Position)]

# Add ungapped position to high_cons
high_cons <- merge(high_cons, pos_consensus, 
                   by.x = c("Gene", "Position"), 
                   by.y = c("Gene", "Aligned_Position"), 
                   all.x = TRUE)

stopifnot(sum(!is.na(high_cons$Residue_Index)) > 0)

# Merge with embeddings (keeping only embedding columns + keys)
embed_cols <- grep("^embedding_", names(embeds_data), value = TRUE)
embeds_dt <- as.data.table(embeds_data)[, c("ID", "Gene", "Residue_Index", embed_cols), with = FALSE]

# Merge on Gene + Residue_Index (will get multiple rows per position if multiple sequences)
merged_embeds <- merge(high_cons, embeds_dt, by = c("Gene", "Residue_Index"), all.x = TRUE)
stopifnot(nrow(merged_embeds) > 0)
message("Merged: ", nrow(merged_embeds), " rows with embeddings")

# ---- 7. PCA ON EMBEDDINGS ----
embed_matrix <- as.matrix(merged_embeds[, ..embed_cols])
complete_rows <- complete.cases(embed_matrix)
stopifnot(sum(complete_rows) > 100)

#uncomment this and run as a script for a sizeable dataset
#pca_result <- prcomp(embed_matrix[complete_rows, ], center = TRUE, scale. = TRUE)
#saveRDS(pca_result, "data/tmp/pca.rds")
pca <- readRDS("data/tmp/pca.rds")
pca_result <- pca
#quit()# Variance explained
var_explained <- summary(pca_result)$importance[2, 1:100]
barplot(var_explained, main="variance explained")
message("Top 100 PCs explain: ", round(sum(var_explained) * 100, 1), "% variance")

# Add PCs to data
n_pcs <- 100  # adjust as needed
pc_scores <- as.data.table(pca_result$x[, 1:n_pcs])
setnames(pc_scores, paste0("PC", 1:n_pcs))

# Rebuild with PCs
merged_final <- cbind(merged_embeds[complete_rows, ], pc_scores)

message("Final dataset: ", nrow(merged_final), " rows, ", n_pcs, " PCs added")
message("Ready to analyze PCs vs structural properties")

library(ggplot2)
library(patchwork)

# ---- PREP VARIABLES ----
# Extract species from Taxonomy if available, or from ID
if ("Taxonomy" %in% names(merged_final)) {
  merged_final[, Species := sub(".*Taxonomy:", "", Taxonomy)]
  merged_final[, Species := sub("_.*", "", Species)]
} else {
  # Try to get from embeds_data
  tax_lookup <- unique(as.data.table(embeds_data)[, .(ID, Taxonomy)])
  merged_final <- merge(merged_final, tax_lookup, by = "ID", all.x = TRUE)
  merged_final[, Species := sub(".*:", "", sub("_.*", "", Taxonomy))]
}

# GWAS class (reuse logic from original)
merged_final[, gwas_class := "not_sig"]
merged_final[P_aa_only < quantile(P_aa_only, 0.25, na.rm = TRUE), gwas_class := "sig_no_control"]
merged_final[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05, na.rm = TRUE), gwas_class := "sig_with_control"]
merged_final[P_aa_with_pcs < quantile(P_aa_with_pcs, 0.05, na.rm = TRUE) & 
               P_aa_only < quantile(P_aa_only, 0.25, na.rm = TRUE), gwas_class := "sig_both"]

merged_final[, gwas_class := factor(gwas_class, 
                                    levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

# Association strength (-log10 p)
merged_final[, log_p_aa := -log10(P_aa_only + 1e-300)]
merged_final[, log_p_controlled := -log10(P_aa_with_pcs + 1e-300)]

# Residue as factor
#merged_final[, Residue := factor(Residue)]

# ---- PLOTTING FUNCTION ----
plot_pca <- function(data, color_var, title, point_size = 0.8, alpha = 0.6) {
  p <- ggplot(data, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = point_size, alpha = alpha) +
    labs(title = title, color = color_var) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10, face = "bold"))
  
  # Use viridis for continuous, default for discrete
  if (is.numeric(data[[color_var]])) {
    p <- p + scale_color_viridis_c()
  }
  p
}

# ---- PANEL PLOTS ----
p_gene <- plot_pca(merged_final, "Gene", "By Gene")
p_gene
p_species <- plot_pca(merged_final, "Species", "By Species")
#p_residue <- plot_pca(merged_final, "Residue", "By Residue") + 
#  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))
p_q3 <- plot_pca(merged_final, "consensus_q3", "By Q3 Structure") +
  scale_color_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A"))
p_q8 <- plot_pca(merged_final, "consensus_q8", "By Q8 Structure")
p_gwas <- plot_pca(merged_final, "gwas_class", "By GWAS Class") +
  scale_color_manual(values = c("not_sig" = "gray70", "sig_no_control" = "gold", 
                                "sig_with_control" = "steelblue", "sig_both" = "darkred"))
p_pval <- plot_pca(merged_final, "log_p_aa", "By -log10(P) uncorrected")
p_pval_ctrl <- plot_pca(merged_final, "log_p_controlled", "By -log10(P) controlled")
p_rsa <- plot_pca(merged_final, "mean_rsa", "By RSA")
p_disorder <- plot_pca(merged_final, "mean_disorder", "By Disorder")



# ---- COMBINE ----
combined <- (p_gene | p_species ) /
  (p_q3 | p_q8 | p_gwas) /
  (p_pval | p_pval_ctrl | p_rsa)

combined <- combined + plot_annotation(
  title = "Embedding PCA colored by biological features",
  theme = theme(plot.title = element_text(size = 14, face = "bold"))
)


# ---- CORRELATION MATRIX ----
cor_vars <- c(paste0("PC", 1:10), 
              "log_p_aa", "log_p_controlled", 
              "mean_rsa", "mean_asa", "mean_disorder", "consensus_freq_q3")

cor_matrix <- cor(merged_final[, ..cor_vars], use = "pairwise.complete.obs")

# Round for readability
cor_round <- round(cor_matrix, 3)

# Print full matrix
options(width = 150)
print(cor_round)

library(pheatmap)
pheatmap(cor_round)

pheatmap(cor_round,
         display_numbers = cor_round,
         number_format = "%.2f",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "PC vs Feature Correlations")

# ---- TOP CORRELATIONS (excluding self) ----
cor_dt <- as.data.table(as.table(cor_matrix))
setnames(cor_dt, c("Var1", "Var2", "r"))
cor_dt <- cor_dt[Var1 != Var2]
cor_dt <- cor_dt[!duplicated(paste(pmin(Var1, Var2), pmax(Var1, Var2)))]  # remove duplicates
cor_dt[, abs_r := abs(r)]
setorder(cor_dt, -abs_r)

message("\n=== Top 20 Correlations ===")
print(cor_dt[1:20], digits = 3)

# ---- PCs vs FEATURES ONLY ----
pc_names <- paste0("PC", 1:10)
feature_names <- c("log_p_aa", "log_p_controlled", "mean_rsa", "mean_asa", "mean_disorder", "consensus_freq_q3")

pc_vs_features <- cor_matrix[pc_names, feature_names]
message("\n=== PCs vs Features ===")
print(round(pc_vs_features, 3))



ggsave("results/pca_panels.pdf", combined, width = 30, height = 30)
ggsave("results/pca_panels.png", combined, width = 16, height = 14, dpi = 150)

message("Saved: results/pca_panels.pdf")
combined


unique(merged_final$Gene)
cor(merged_final$log_p_controlled, merged_final$log_p_aa)
ggplot(merged_final, aes(x = PC3, y = log_p_controlled)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(trans = "log10") +
  labs(x = "PC3", y = "-log10(P) uncorrected",
       title = paste("r =", round(cor(merged_final$PC3, merged_final$log_p_aa, use = "complete.obs"), 3))) +
  theme_minimal()

nrow(merged_final)

# Check before merge
message("high_cons genes: ", paste(unique(high_cons$Gene), collapse = ", "))
message("embeds_dt genes: ", paste(unique(embeds_dt$Gene), collapse = ", "))

# The issue: merge keeps embeds_dt$Gene if there's a mismatch or duplicate column handling
# Fix: rename or drop Gene from embeds_dt before merge

embed_cols <- grep("^embedding_", names(embeds_data), value = TRUE)
embeds_dt <- as.data.table(embeds_data)[, c("ID", "Gene", "Residue_Index", embed_cols), with = FALSE]

# Check if Gene values match
setnames(embeds_dt, "Gene", "Gene_embed")  # rename to avoid collision

merged_embeds <- merge(high_cons, embeds_dt, 
                       by.x = c("Gene", "Residue_Index"), 
                       by.y = c("Gene_embed", "Residue_Index"), 
                       all.x = TRUE)

# Verify
stopifnot(length(unique(merged_embeds$Gene)) > 1)
message("Genes after merge: ", length(unique(merged_embeds$Gene)))

# Diagnostic
message("high_cons genes: ", paste(sort(unique(high_cons$Gene)), collapse = ", "))
message("embeds genes: ", paste(sort(unique(embeds_dt$Gene)), collapse = ", "))
message("Overlap: ", length(intersect(unique(high_cons$Gene), unique(embeds_dt$Gene))))
message("Unique genes in parquet: ", paste(unique(embeds_data$Gene), collapse = ", "))
message("Total rows: ", nrow(embeds_data))

ggsave(filename="results/aaa_pc3_pvals_cor.py")
