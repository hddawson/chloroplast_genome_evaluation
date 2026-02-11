library(arrow)
library(data.table)

embeds <- read_parquet("data/all_gene_embeddings.parquet")


embed_cols <- grep("embe", colnames(embeds))

#pca <- prcomp(embeds[,embed_cols])
#saveRDS(pca, "data/tmp/pca_results.rds")
pca <- readRDS("data/tmp/pca_results.rds")

pvar <- round(pca$sdev^2 * 100 / sum(pca$sdev^2),3)
barplot(pvar[1:10], main="% var explained")

plot(pca$x[,1],pca$x[,2])

cols <- as.factor(embeds$Gene)
plot(pca$x[,1], pca$x[,2], col=cols, pch=16,
     xlab=paste0("PCA1, ",pvar[1],"%"),
     ylab=paste0("PCA2, ",pvar[2],"%"),
     main="PCA on gene embeddings |esmc_300m|61 genes, present in>95% species|Length-based outlier filtration|")
#legend("topright", legend=levels(cols), col=1:length(levels(cols)), pch=16)
plot(pca$x[,2], pca$x[,3], col=cols, pch=16,
     xlab=paste0("PCA2, ",pvar[2],"%"),
     ylab=paste0("PCA3, ",pvar[3],"%"),
     main="PCA on gene embeddings |esmc_300m|61 genes, present in>95% species|Length-based outlier filtration|")


prefix <- substr(embeds$Gene, 1, 3)
prefix_factor <- as.factor(prefix)

shapes <- c(16, 17, 15, 3, 4, 8, 18, 1)
shapes <- shapes[(as.numeric(prefix_factor) - 1) %% length(shapes) + 1]

cols <- rainbow(length(levels(prefix_factor)))[as.numeric(prefix_factor)]


plot(pca$x[,1], pca$x[,2], col=cols, pch=shapes,
     xlab=paste0("PCA1, ",pvar[1],"%"),
     ylab=paste0("PCA2, ",pvar[2],"%"),
     main="PCA on gene embeddings |esmc_300m|61 genes, present in>95% species|Length-based outlier filtration|")

# Create 61 plots - one for each gene highlighted in red
genes <- unique(embeds$Gene)

# Remove missing values
valid_idx <- complete.cases(pca$x[,1], pca$x[,2], embeds$Gene)
pc1 <- pca$x[valid_idx, 1]
pc2 <- pca$x[valid_idx, 2]
gene_names <- embeds$Gene[valid_idx]

cat("runnin!")
# Create plots for each gene
for(i in 1:length(genes)) {
  gene <- genes[i]
  cat(gene, "\n")
  
  # Color: red for target gene, gray for all others
  png(paste0("plots/embeds_pca_", gene, ".png"))
  # Plot gray points first, then red on top
  plot(pc1, pc2,
       col = "gray70",
       pch = 16,
       cex = 1.2,
       xlab = paste0("PC1 (", round(pvar[1], 1), "%)"),
       ylab = paste0("PC2 (", round(pvar[2], 1), "%)"),
       main = paste("Gene:", gene))
  
  # Add red points for target gene on top
  target_idx <- gene_names == gene
  points(pc1[target_idx], pc2[target_idx], col = "red", pch = 16, cex = 1.2)
  
  # Optional: save plots
  dev.off()
}
