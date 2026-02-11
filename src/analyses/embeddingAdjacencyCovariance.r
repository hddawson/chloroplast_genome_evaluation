library(arrow)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Load mean-pooled embeddings (already have)
embeds <- read_parquet("data/all_gene_embeddings.parquet")
stopifnot(nrow(embeds) > 0)
cat("Loaded", nrow(embeds), "protein embeddings\n")

# Aggregate to gene-level (mean across species)
embed_cols <- grep("^embedding_", names(embeds), value = TRUE)

gene_embeddings <- embeds |>
  group_by(Gene) |>
  summarise(across(all_of(embed_cols), \(x) mean(x, na.rm = TRUE))) |>
  ungroup()

stopifnot(nrow(gene_embeddings) > 0)
cat("Aggregated to", nrow(gene_embeddings), "genes\n")

# Build adjacency matrix (same 3-letter prefix = same complex)
genes <- gene_embeddings$Gene
prefixes <- substr(genes, 1, 3)

adj_matrix <- outer(prefixes, prefixes, FUN = "==") * 1
diag(adj_matrix) <- NA
rownames(adj_matrix) <- colnames(adj_matrix) <- genes

cat("Unique complexes:", length(unique(prefixes)), "\n")

# Cosine similarity matrix
embed_mat <- gene_embeddings |>
  select(all_of(embed_cols)) |>
  as.matrix()
rownames(embed_mat) <- genes

norm_mat <- embed_mat / sqrt(rowSums(embed_mat^2))
cos_sim <- tcrossprod(norm_mat)
diag(cos_sim) <- NA

stopifnot(all(dim(cos_sim) == dim(adj_matrix)))

# Correlation test
upper_idx <- upper.tri(adj_matrix)
adj_vec <- adj_matrix[upper_idx]
sim_vec <- cos_sim[upper_idx]

cor_test <- cor.test(sim_vec, adj_vec)
cat("\n=== Cosine Similarity vs PPI ===\n")
cat("Pearson r:", round(cor_test$estimate, 4), "\n")
cat("p-value:", format.pval(cor_test$p.value, digits = 3), "\n")

# Summary stats
cat("\nWithin-complex similarity:", round(mean(sim_vec[adj_vec == 1]), 4), "\n")
cat("Between-complex similarity:", round(mean(sim_vec[adj_vec == 0]), 4), "\n")


# Order genes by complex for visualization
gene_order <- genes[order(prefixes)]

# 1. Heatmaps side by side
pdf("ppi_heatmaps.pdf", width = 14, height = 6)

par(mfrow = c(1, 2))

# Similarity matrix ordered by complex
pheatmap(cos_sim[gene_order, gene_order],
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Cosine Similarity (ordered by complex)",
         fontsize_row = 5, fontsize_col = 5) -> p1

# Clustered - does it recover complexes?
pheatmap(cos_sim,
         main = "Cosine Similarity (clustered)",
         fontsize_row = 5, fontsize_col = 5,
         annotation_row = data.frame(Complex = prefixes, row.names = genes),
         annotation_col = data.frame(Complex = prefixes, row.names = genes)) -> p2

dev.off()

# 2. Boxplot
pdf("ppi_boxplot.pdf", width = 5, height = 4)

boxplot(sim_vec[adj_vec == 0], sim_vec[adj_vec == 1],
        names = c("Between", "Within"),
        ylab = "Cosine Similarity",
        main = sprintf("r = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value),
        col = c("gray80", "steelblue"))

dev.off()

cat("Saved: ppi_heatmaps.pdf, ppi_boxplot.pdf\n")

# Create pairwise dataframe with cosine similarity and ground truth
pairs_df <- expand.grid(gene_a = genes, gene_b = genes, stringsAsFactors = FALSE) |>
  filter(gene_a < gene_b)  # upper triangle only, no duplicates

pairs_df$cosine_sim <- mapply(function(a, b) {
  cos_sim[a, b]
}, pairs_df$gene_a, pairs_df$gene_b)

pairs_df$same_complex <- mapply(function(a, b) {
  adj_matrix[a, b]
}, pairs_df$gene_a, pairs_df$gene_b)

stopifnot(!any(is.na(pairs_df$cosine_sim)))
stopifnot(all(pairs_df$same_complex %in% c(0, 1)))

write.csv(pairs_df, "data/gene_pairs_cosine.csv", row.names = FALSE)
cat("Exported", nrow(pairs_df), "pairs\n")
head(pairs_df)

# python ../agent/src/agentExperiment...

library(ggplot2)
library(pROC)

results <- read.csv("../agent/data/agent_experiment_results.csv")
stopifnot(nrow(results) > 0)
table(results$same_complex)
# ROC curves
roc_without <- roc(results$same_complex, results$pred_without)
roc_with <- roc(results$same_complex, results$pred_with)

pdf("agent_roc_comparison.pdf", width = 6, height = 5)

plot(roc_without, col = "gray50", main = "Agent PPI Prediction: With vs Without Cosine Sim")
plot(roc_with, col = "steelblue", add = TRUE)
legend("bottomright", 
       legend = c(sprintf("Without cosine (AUC=%.2f)", auc(roc_without)),
                  sprintf("With emb. cosine sim (AUC=%.2f)", auc(roc_with))),
       col = c("gray50", "steelblue"), lwd=1)

dev.off()


plot(roc_without, col = "gray50", main = "PPI Prediction Comparison")
plot(roc_with, col = "steelblue", add = TRUE)
legend("bottomright", 
       legend = c(sprintf("Agent without cosine (AUC=%.2f)", auc(roc_without)),
                  sprintf("Agent with cosine (AUC=%.2f)", auc(roc_with))),
       col = c("gray50", "steelblue"), 
       lty = c(1, 1, 2), lwd = 2)

cat("AUC without:", round(auc(roc_without), 3), "\n")
cat("AUC with:", round(auc(roc_with), 3), "\n")

# Check prediction distributions
cat("Without cosine:\n")
print(table(results$pred_without))

cat("\nWith cosine:\n")
print(table(round(results$pred_with, 1)))

# Is the model just using cosine as a threshold?
cat("\nCorrelation: pred_with vs cosine_sim:", 
    round(cor(results$pred_with, results$cosine_sim), 3), "\n")

# Baseline: cosine sim alone as predictor
roc_cosine <- roc(results$same_complex, results$cosine_sim)
cat("AUC cosine alone:", round(auc(roc_cosine), 3), "\n")

