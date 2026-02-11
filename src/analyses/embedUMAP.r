library(arrow)
library(data.table)
library(dplyr)
library(uwot)
library(ggplot2)

# Load data
pca <- read_parquet("results/embeddings_full_pca.parquet")
data <- as.data.table(read_parquet("data/processed_data.parquet"))

# Extract PC columns
colnames(pca)
pc_cols <- grep("^PC[0-9]+$", names(pca), value = TRUE)
print(pc_cols)
stopifnot(length(pc_cols) > 0)

# Run UMAP on PC coordinates
set.seed(42)
umap_res <- umap(as.matrix(pca[, pc_cols]), n_neighbors = 15, min_dist = 0.1, n_threads = 4)

# Add UMAP coords to data
pca$UMAP1 <- umap_res[, 1]
pca$UMAP2 <- umap_res[, 2]

# Merge metadata
pca <- pca %>% left_join(data[, .(ID, Organism, Order)], by = "ID")

# Plot by Gene
p_gene <- ggplot(pca, aes(UMAP1, UMAP2, color = Gene)) +
  geom_point(size = 0.1, alpha = 0.3) +
  theme_minimal() +
  ggtitle("UMAP colored by Gene")

# Plot by Residue (site class)
p_residue <- ggplot(pca, aes(UMAP1, UMAP2, color = Residue)) +
  geom_point(size = 0.1, alpha = 0.3) +
  theme_minimal() +
  ggtitle("UMAP colored by Residue")

# Save
ggsave("umap_by_gene.png", p_gene, width = 10, height = 8)
ggsave("umap_by_residue.png", p_residue, width = 10, height = 8)

# If you have significance column, add:
# p_sig <- ggplot(pca, aes(UMAP1, UMAP2, color = significance)) + ...