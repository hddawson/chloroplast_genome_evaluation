library(arrow)
library(data.table)

# 1. Load the supermatrix
supermatrix <- as.data.table(read_parquet("data/mergedAln.pq"))

# Optional: keep sample IDs as rownames if parquet saved them in an index column
# If your Python code stored IDs as the index, arrow might load them as the first column
# Assuming first column is sample ID:
#setnames(supermatrix, index, "sample_id")
#setkey(supermatrix, sample_id)

# 2. Convert each column to binary: major allele OR gap = 1, else 0
for (col in names(supermatrix)[-1]) {
  vals <- supermatrix[[col]]
  tab <- table(vals)
  major <- names(sort(tab[names(tab) != "-"], decreasing = TRUE))[1]
  if (length(major) == 0) major <- "-"  # all gaps case
  supermatrix[[col]] <- as.integer(vals == major | vals == "-")
}

# 3. Prepare numeric matrix for PCA
nzv_cols <- supermatrix[, which(apply(.SD, 2, function(x) var(as.numeric(factor(x))) > 0))]
supermatrix <- supermatrix[, ..nzv_cols]
mat <- as.matrix(supermatrix[, -1, with = FALSE])
rownames(mat) <- supermatrix$index



# 4. PCA (center & scale optional for binary)
pca_res <- prcomp(mat, center = TRUE, scale. = TRUE, rank.=2)
saveRDS(pca_res,"data/pca.RDS")
# 5. Plot PCA (base R)
plot(
  pca_res$x[, 1], pca_res$x[, 2],
  xlab = "PC1", ylab = "PC2",
  main = "PCA on Binary Encoded Supermatrix"
)
text(
  pca_res$x[, 1], pca_res$x[, 2],
  labels = rownames(mat),
  pos = 3, cex = 0.7
)
