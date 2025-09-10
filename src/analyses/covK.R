library(data.table)
library(BGLR)
library(parallel)
library(arrow)

df        <- as.data.table(read_parquet("data/processed_data.parquet"))
y         <- df$pheno_Topt_site_p50
data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X <- as.matrix(data[,-1])
pca_X <- prcomp(X, scale. = TRUE, rank. = 100)

# ─── calculate kinship matrix ────────────────────────────────────────────────
cat("Computing kinship matrix K...\n")
X_centered <- scale(X, center=TRUE, scale=FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

cat("Computing Cov(K)...\n")
K_centered <- scale(K, center=TRUE, scale=FALSE)
CovK <- tcrossprod(K_centered) / (nrow(K_centered) - 1)

# PCA of K
cat("Performing PCA of K...\n")
pca_K <- eigen(K_centered)
eigenvectors <- pca_K$vectors
eigenvalues  <- pca_K$values

save(K, CovK, pca_K, pca_X, eigenvalues, file="data/kinship_analysis.RData")
load("data/kinship_analysis.RData")

library(ggplot2)
library(viridis)
library(corrplot)

# Load required libraries
library(ggplot2)
library(viridis)
library(corrplot)

# Load data
load("data/kinship_analysis.RData")


# 1. PCA Eigenvalue plots
png("plots/pca_eigenvalues.png", width=800, height=600)
par(mfrow=c(2,2))

# Scree plot
plot(eigenvalues[1:min(50, length(eigenvalues))], 
     type="b", pch=19, 
     xlab="Principal Component", 
     ylab="Eigenvalue",
     main="Scree Plot")

# Cumulative variance explained
cum_var <- cumsum(eigenvalues) / sum(eigenvalues)
plot(cum_var[1:min(50, length(cum_var))], 
     type="b", pch=19,
     xlab="Principal Component", 
     ylab="Cumulative Proportion",
     main="Cumulative Variance Explained")


par(mfrow=c(1,2))

# Create PDF
pdf("plots/pca_analysis_report.pdf", width=12, height=8)
var_explained <- eigenvalues / sum(eigenvalues) * 100
cum_var <- cumsum(var_explained)
# Page 1: Variance explained
par(mfrow=c(1,2))
barplot(var_explained[1:10], 
        names.arg=paste0("PC", 1:10),
        xlab="Principal Component", 
        ylab="Variance Explained (%)",
        main="Individual Variance Explained",
        col="steelblue")

plot(cum_var[1:min(20, length(cum_var))], 
     type="b", pch=19,
     xlab="Principal Component", 
     ylab="Cumulative Variance (%)",
     main="Cumulative Variance Explained",
     col="darkred")
grid()

# Pages 2-6: First 10 PC pairs (2 pairs per page)
pc_pairs <- list(
  c(1,2), c(3,4), c(5,6), c(7,8), c(9,10),
  c(1,3), c(2,4), c(1,5), c(2,6), c(3,7)
)

for(i in seq(1, 10, 2)) {
  par(mfrow=c(1,2))
  
  # First pair
  pc1 <- pc_pairs[[i]][1]
  pc2 <- pc_pairs[[i]][2]
  plot(eigenvectors[,pc1], eigenvectors[,pc2],
       xlab=paste0("PC", pc1, " (", round(var_explained[pc1], 1), "%)"), 
       ylab=paste0("PC", pc2, " (", round(var_explained[pc2], 1), "%)"),
       main=paste0("PC", pc1, " vs PC", pc2),
       pch=19, cex=0.4, col="steelblue")
  grid()
  
  # Second pair (if exists)
  if(i+1 <= 10) {
    pc1 <- pc_pairs[[i+1]][1]
    pc2 <- pc_pairs[[i+1]][2]
    plot(eigenvectors[,pc1], eigenvectors[,pc2],
         xlab=paste0("PC", pc1, " (", round(var_explained[pc1], 1), "%)"), 
         ylab=paste0("PC", pc2, " (", round(var_explained[pc2], 1), "%)"),
         main=paste0("PC", pc1, " vs PC", pc2),
         pch=19, cex=0.4, col="darkgreen")
    grid()
  }
}

dev.off()





# First few PCs
if(ncol(eigenvectors) >= 2) {
  plot(eigenvectors[,1], eigenvectors[,2],
       xlab="PC1", ylab="PC2", 
       main="PC1 vs PC2",
       pch=19, cex=0.5)
}
dev.off()

plot(eigenvectors[,2], eigenvectors[,3],
     xlab="PC2", ylab="PC3", 
     main="PC2 vs PC3",
     pch=19, cex=0.5)

# 2. Kinship matrix visualization (efficient for large matrices)
n <- nrow(K)
cat(sprintf("Kinship matrix size: %d x %d (%d elements)\n", n, n, n^2))

# Sample subset for visualization if too large
max_display <- 2000
if(n > max_display) {
  idx <- sample(n, max_display)
  K_sub <- K[idx, idx]
  cat(sprintf("Sampling %d x %d subset for visualization\n", max_display, max_display))
} else {
  K_sub <- K
}

# Kinship matrix heatmap
png("plots/kinship_heatmap.png", width=800, height=800)
image(K_sub, 
      col=viridis(100), 
      main="Kinship Matrix",
      xlab="Individual", ylab="Individual")
dev.off()

# Kinship distribution
png("plots/kinship_distribution.png", width=800, height=600)
par(mfrow=c(2,2))

# Overall distribution
hist(K[upper.tri(K)], 
     breaks=50, 
     main="Kinship Coefficient Distribution",
     xlab="Kinship Coefficient")

# Diagonal (self-kinship)
hist(diag(K), 
     breaks=30,
     main="Self-Kinship (Diagonal)",
     xlab="Self-Kinship Coefficient")

# Off-diagonal distribution
off_diag <- K[upper.tri(K)]
hist(off_diag, 
     breaks=50,
     main="Pairwise Kinship",
     xlab="Kinship Coefficient")

# Summary statistics
boxplot(list(Diagonal=diag(K), 
             OffDiagonal=off_diag),
        main="Kinship Summary")
dev.off()

# 3. Combined PCA-Kinship summary
png("plots/pca_kinship_summary.png", width=1200, height=800)
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=TRUE))

# Top eigenvalues
barplot(eigenvalues[1:min(10, length(eigenvalues))],
        names.arg=1:min(10, length(eigenvalues)),
        main="Top 10 Eigenvalues",
        xlab="PC", ylab="Eigenvalue")

# PC loadings heatmap (first few PCs)
if(ncol(eigenvectors) >= 5) {
  image(t(eigenvectors[,1:5]), 
        col=colorRampPalette(c("blue","white","red"))(50),
        main="PC Loadings (1-5)",
        xlab="PC", ylab="Individual")
}

# Kinship matrix corner
corner_size <- min(100, nrow(K))
image(K[1:corner_size, 1:corner_size],
      col=viridis(50),
      main=sprintf("Kinship Matrix (%dx%d corner)", corner_size, corner_size))









