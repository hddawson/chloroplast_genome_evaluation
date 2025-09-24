library(arrow)
library(data.table)

# Load data
embeds <- as.data.table(read_parquet("data/all_gene_embeddings.parquet"))

# Get embedding columns and genes
embed_cols <- grep("embed", colnames(embeds), value=TRUE)
genes <- unique(embeds$Gene)

# Initialize columns for MDS results
embeds[, c("MDS1", "MDS2", "MDS_zscore") := list(NA_real_, NA_real_, NA_real_)]

# Initialize outlier tracking
all_outliers <- c()

# Process each gene for MDS outlier detection
cat("Processing", length(genes), "genes for MDS outlier detection...\n")

for(gene in genes) {
  cat("Processing gene:", gene, "\n")
  
  # Filter for specific gene
  X <- embeds[Gene == gene]
  X_embeds <- X[, ..embed_cols]
  X_matrix <- as.matrix(X_embeds)
  
  # MDS Analysis
  dist_matrix <- dist(X_matrix)
  mds_result <- cmdscale(dist_matrix, k=2)
  
  # Calculate outliers
  mds_df <- data.frame(
    MDS1 = mds_result[,1],
    MDS2 = mds_result[,2]
  )
  
  d <- sqrt(rowSums(scale(mds_df[,c("MDS1","MDS2")], center=TRUE, scale=FALSE)^2))
  z <- scale(d)
  is_outlier <- abs(z) > 3
  
  # Add MDS results back to original dataframe
  gene_rows <- which(embeds$Gene == gene)
  embeds[gene_rows, MDS1 := mds_result[,1]]
  embeds[gene_rows, MDS2 := mds_result[,2]]
  embeds[gene_rows, MDS_zscore := as.numeric(z)]
  
  # Track outliers with original row indices
  gene_outlier_rows <- gene_rows[is_outlier]
  all_outliers <- c(all_outliers, gene_outlier_rows)
  
  # Create MDS plot
  png(paste0("plots/mds_", gene, ".png"))
  plot(mds_result[,1], mds_result[,2],
       col = ifelse(is_outlier, "red", "gray70"),
       pch = 16,
       cex = 1.2,
       xlab = "MDS1",
       ylab = "MDS2",
       main = paste("MDS -", gene, "| Outliers (z>3) in red"))
  dev.off()
  
  cat("  Found", sum(is_outlier), "outliers\n")
}

#okay, 3 works well as a default, but the following genes should have 
#let's go through each gene
genes <- unique(embeds$Gene)
# Manual cutoffs (fill in as needed)
cutoffs <- data.table(
  Gene = genes,
  Cutoff = rep(3, length(genes)) # default 3, adjust manually
)

# Apply manual cutoffs and plot
for (gene in genes) {
  cval <- cutoffs[Gene == gene, Cutoff]
  X <- embeds[Gene == gene]
  z <- X$MDS_zscore
  
  col_assign <- ifelse(abs(z) > 6, "black",
                       ifelse(abs(z) > 5, "blue",
                              ifelse(abs(z) > 4, "purple",
                                     ifelse(abs(z) > 3, "red",
                                            ifelse(abs(z) > 2, "orange",
                                                   ifelse(abs(z) > 1, "gold", "gray70"))))))
  
  plot(X$MDS1, X$MDS2,
       col = col_assign,
       pch = 16,
       cex = 1.2,
       xlab = "MDS1",
       ylab = "MDS2",
       main = paste("MDS -", gene, "| Manual cutoff =", cval))
}
cutoffs[(which(cutoffs$Gene=="ycf4")), "Cutoff"] <- 2
cutoffs[(which(cutoffs$Gene=="ycf3")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps8")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps7")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps4")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps3")), "Cutoff"] <- 2
cutoffs[(which(cutoffs$Gene=="rps19")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps18")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rps14")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="rps11")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpoC1")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpoB")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpoB")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpl36")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpl33")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpl2")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpl23")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="rpl20")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rpl16")), "Cutoff"] <- 1
cutoffs[(which(cutoffs$Gene=="rpl14")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="rbcL")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="psbZ")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="psbT")), "Cutoff"] <- 6
cutoffs[(which(cutoffs$Gene=="psbN")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="psbM")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="psbK")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbJ")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbI")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbH")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbF")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbE")), "Cutoff"] <- 6
cutoffs[(which(cutoffs$Gene=="psbD")), "Cutoff"] <- 6
cutoffs[(which(cutoffs$Gene=="psbB")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psbA")), "Cutoff"] <- 6
cutoffs[(which(cutoffs$Gene=="psaJ")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psaC")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psaB")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="psaA")), "Cutoff"] <- 6
cutoffs[(which(cutoffs$Gene=="petN")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="petG")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="petA")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhK")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhJ")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhI")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhH")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhG")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhE")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhD")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="ndhC")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ndhB")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="ndhA")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="matK")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="cemA")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="ccsA")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="atpI")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="atpH")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="atpF")), "Cutoff"] <- 5
cutoffs[(which(cutoffs$Gene=="atpE")), "Cutoff"] <- 4
cutoffs[(which(cutoffs$Gene=="atpB")), "Cutoff"] <- 3
cutoffs[(which(cutoffs$Gene=="atpA")), "Cutoff"] <- 4


embeds <- merge(embeds, cutoffs, by="Gene", all.x=TRUE)
hist(embeds$MDS_zscore)
hist(embeds$Cutoff)
embeds[, ManualOutlier := abs(MDS_zscore) >= Cutoff] 

all_outliers <- which(embeds$ManualOutlier)

cat("Total outliers found:", length(all_outliers), "\n")
clean_embeds <- embeds[ManualOutlier == FALSE]

cat("Samples remaining:", nrow(clean_embeds), "of", nrow(embeds), "\n")

saveRDS(embeds, "data/tmp/embeds_with_mds.rds")


# Remove outliers from dataset
cat("Total outliers found:", length(all_outliers), "\n")
cat("Samples remaining:", nrow(clean_embeds), "of", nrow(embeds), "\n")
saveRDS(embeds, "data/tmp/embeds_with_mds.rds")
# PCA on cleaned data
embeds <- readRDS("data/tmp/embeds_with_mds.rds")
clean_embeds <- embeds[ManualOutlier == FALSE]
cat("Performing PCA on cleaned data...\n")
pca_embed_cols <- grep("embe", colnames(clean_embeds))
pca <- prcomp(clean_embeds[, ..pca_embed_cols])
saveRDS(pca, "data/tmp/pca_results_cleaned.rds")
pca <- readRDS("data/tmp/pca_results_cleaned.rds")
pvar <- round(pca$sdev^2 * 100 / sum(pca$sdev^2), 3)
dim(pca$x)
# Save PCA results

# PCA visualizations
png("plots/pca_variance.png")
barplot(pvar[1:10], main="% var explained (cleaned data)")
dev.off()

# Overall PCA plot
png("plots/pca_overall.png")
cols <- as.factor(clean_embeds$Gene)
par(mfrow=c(1,2))
plot(pca$x[,1], pca$x[,2], col=cols, pch=16,
     xlab=paste0("PC1, ",pvar[1],"%"),
     ylab=paste0("PC2, ",pvar[2],"%"),
     main="PCA on cleaned gene embeddings")
plot(pca$x[,2], pca$x[,3], col=cols, pch=16,
     xlab=paste0("PC2, ",pvar[2],"%"),
     ylab=paste0("PC3, ",pvar[3],"%"),
     main="")
dev.off()
# Gene-specific PCA plots
valid_idx <- complete.cases(pca$x[,1], pca$x[,2], clean_embeds$Gene)
pc1 <- pca$x[valid_idx, 1]
pc2 <- pca$x[valid_idx, 2]
gene_names <- clean_embeds$Gene[valid_idx]
clean_genes <- unique(clean_embeds$Gene)

cat("Creating", length(clean_genes), "gene-specific PCA plots...\n")

for(gene in clean_genes) {
  png(paste0("plots/pca_", gene, ".png"))
  
  plot(pc1, pc2,
       col = "gray70",
       pch = 16,
       cex = 1.2,
       xlab = paste0("PC1 (", round(pvar[1], 1), "%)"),
       ylab = paste0("PC2 (", round(pvar[2], 1), "%)"),
       main = paste("PCA - Gene:", gene, "(cleaned)"))
  
  target_idx <- gene_names == gene
  points(pc1[target_idx], pc2[target_idx], col = "red", pch = 16, cex = 1.2)
  
  dev.off()
}

cat("Workflow complete. Plots saved in plots/\n")
cat("Cleaned PCA results saved to data/tmp/pca_results_cleaned.rds\n")
cat("Original embeds dataframe now contains MDS1, MDS2, and MDS_zscore columns\n")