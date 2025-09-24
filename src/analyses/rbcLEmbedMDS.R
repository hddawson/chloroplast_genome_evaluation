library(data.table)
library(parallel)
library(arrow)

setwd("/workdir/hdd29/chloroplast_genome_evaluation")
df        <- as.data.table(read_parquet("data/processed_data.parquet"))
par(mfrow=c(1,2))
hist(df$Total_Amino_Acids, main="Total AA per plastome")
hist(df$geno_genomeLength)
par(mfrow=c(1,1))

y <- df$pheno_Topt_site_p50
X <- as.data.table(read_parquet("data/rbcL_embeddings.parquet"))

# Define outliers
outliers <- c("JN032131.1", "MZ475300.1", "MH325132.1", "MK645903.1")

# Merge data
merged <- merge(df[, .(ID, pheno_Topt_site_p50)], X, by.x="ID", by.y="sample_id")

# Linear models
#fit <- lm(pheno_Topt_site_p50 ~ . -ID, data=merged)
#summary(fit)
#plot(fit$fitted.values, merged$pheno_Topt_site_p50)
# 
# aa_cols <- intersect(c("A","R","N","D","C","Q","E","G","H","I",
#                        "L","K","M","F","P","S","T","W","Y","V"),
#                      names(df))
# form <- as.formula(paste("pheno_Topt_site_p50 ~", paste(aa_cols, collapse="+")))
# #fit2 <- lm(form, data=df)
# #summary(fit2)
# #plot(fit2$fitted.values, df$pheno_Topt_site_p50)
# 
# par(mfrow=c(1,2), mar=c(4,4,3,1))
# plot(fit$fitted.values, merged$pheno_Topt_site_p50,
#      xlab="Fitted T_opt (rbcL embeddings)",
#      ylab="Observed T_opt",
#      main="rbcL Embeddings Model")
# plot(fit2$fitted.values, df$pheno_Topt_site_p50,
#      xlab="Fitted T_opt (AA proportions)",
#      ylab="Observed T_opt",
#      main="Amino Acid Composition Model")
# 
# # Correlation analysis
# y <- merged$pheno_Topt_site_p50
# cors <- sapply(merged[, !c("ID","pheno_Topt_site_p50"), with=FALSE], function(x) cor(x, y))
# par(mfrow=c(1,1))
# hist(cors, main="correlations of embeddings with T_opt_site")
# high_cor <- cors[abs(cors) > 0.7]
# high_cor

# MDS Analysis
embedding_cols <- names(X)[names(X) != "sample_id"]
embedding_matrix <- as.matrix(X[, ..embedding_cols])

# Ensure we have embeddings for all samples
dim(embedding_matrix)

# Calculate distance matrix
dist_matrix <- dist(embedding_matrix)

# Perform MDS
mds_result <- cmdscale(dist_matrix, k=2)
saveRDS(mds_result, "results/mds_result.rds")
cat("SAVED!!!")
mds_result <- readRDS("results/mds_result.rds")
# Create color vector for outliers
sample_ids <- X$sample_id
is_outlier <- sample_ids %in% outliers
colors <- ifelse(is_outlier, "red", "black")
point_sizes <- ifelse(is_outlier, 1.5, 1)
png("plots/MDS_rbcL.png")
# Plot MDS with outliers highlighted
plot(mds_result[,1], mds_result[,2], 
     col = colors, 
     pch = 19,
     cex = point_sizes,
     xlab = "MDS Dimension 1", 
     ylab = "MDS Dimension 2",
     main = "MDS of rbcL Embeddings")
points(mds_result[is_outlier,1], mds_result[is_outlier,2],
       col="red", pch=19, cex=1.5)
# Add legend
legend("bottomright", 
       legend = c("Normal", "Outlier"), 
       col = c("black", "red"), 
       pch = 19,
       cex = 0.8)
dev.off()


mds_df <- data.frame(
  MDS1 = mds_result[,1],
  MDS2 = mds_result[,2],
  sample_id = X$sample_id,
  is_outlier = X$sample_id %in% outliers
)

ggplot(mds_df, aes(x = MDS1, y = MDS2, label = sample_id)) +
  geom_point(aes(color = is_outlier), size = 1.5) +
  geom_text(data = subset(mds_df, is_outlier),
            aes(label = sample_id),
            hjust = 0, nudge_x = 0.01, size = 2.5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "MDS of rbcL Embeddings", color = "Outlier") +
  theme_minimal()


subset_df <- subset(mds_df, MDS1 > 0.120)
ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = is_outlier), size = 1.5) +
  geom_text(data = subset_df, aes(label = sample_id),
            hjust = 0, nudge_x = 0.01, size = 4) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "MDS of rbcL Embeddings", color = "Outlier") +
  theme_minimal()

library(dbscan)

mds_df <- data.frame(
  MDS1 = mds_result[,1],
  MDS2 = mds_result[,2],
  sample_id = X$sample_id,
  is_outlier = X$sample_id %in% outliers
)
d <- sqrt(rowSums(scale(mds_df[,c("MDS1","MDS2")], center=TRUE, scale=FALSE)^2))
z <- scale(d)
outliers <- mds_df$sample_id[abs(z) > 3]
sum(z>3)
mds_df$is_outlier <- abs(z) > 4

ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = is_outlier), size = 1.5) +
  geom_text(data = subset(mds_df, is_outlier),
            aes(label = sample_id),
            hjust = 0, nudge_x = 0.01, size = 3) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "MDS of rbcL Embeddings, outliers have z-score of dist from centroid >5", color = "Outlier") +
  theme_minimal()
sum(mds_df$is_outlier)
mds_df[mds_df$is_outlier, "sample_id"]
# "MZ475300.1" "MZ160994.1" "MT876494.1" "JN032131.1" "MK645903.1" "MW817015.1" "MW829763.1" "MH325132.1"

keep_idx <- !X$sample_id %in% c("MZ475300.1","MZ160994.1","MT876494.1","JN032131.1","MK645903.1","MW817015.1","MW829763.1","MH325132.1")
X_clean <- X[keep_idx, ]
embedding_matrix_clean <- as.matrix(X_clean[, ..embedding_cols])
dim(embedding_matrix_clean)




#pca_result <- prcomp(embedding_matrix_clean, scale.=TRUE)
saveRDS(pca_result,"data/tmp/rbcL_PCA.rds")
par(mfrow=c(2,2))
barplot(pca_result$sdev^2 / (sum(pca_result$sdev^2)))
plot(pca_result$x[,1],pca_result$x[,2])
plot(pca_result$x[,2],pca_result$x[,3])
plot(pca_result$x[,3],pca_result$x[,4])
pvar <- round(pca_result$sdev^2 * 100 / (sum(pca_result$sdev^2)),3)
df_clean <- df[df$ID %in% X_clean$sample_id, ]
scores <- as.data.frame(pca_result$x)
scores$sample_id <- X_clean$sample_id
df_map <- df[, c("ID", "pheno_Topt_site_p50")]
scores <- merge(scores, df_map, by.x="sample_id", by.y="ID")
scores$pheno <- df_clean$pheno_Topt_site_p50

colnames(df)[grep("pheno",colnames(df))]
merged_clean <- merge(df[, .(ID, pheno_Topt_site_p50,n_occur)],
                      X_clean, by.x="ID", by.y="sample_id")

y <- merged_clean$geno_genomeLength
cors <- sapply(merged_clean[, !c("ID","geno_genomeLength"), with=FALSE],
               function(x) cor(x, y, use="pairwise.complete.obs"))

hist(cors, main="correlation of embeddings with genomic length")
best <- names(cors)[which.max(cors)]
plot(merged_clean[[best]], y,
     xlab=best, ylab="genomic GC",
     main=paste("Best correlated embedding:", best))

df$rbcL_CDS[1:2]


pheno_cols <- colnames(df)[grep("pheno", colnames(df))]
res <- rbindlist(lapply(pheno_cols, function(p) {
  m <- merge(df[, c("ID", p), with=FALSE], X_clean, by.x="ID", by.y="sample_id")
  y <- m[[p]]
  cors <- sapply(m[, !c("ID", p), with=FALSE], function(x) cor(x, y, use="pairwise.complete.obs"))
  data.table(pheno=p, embedding=names(cors), cor=cors)
}))

hist(res$cor, main="Correlations of rbcL embeddings with all wclim variables")
length(unique(res$embedding))

hist(res$cor[grep("occu", res$pheno)], main="Correlation between embeddings and n_occurrences")
library(Biostrings)

keep <- !is.na(df$rbcL_CDS) & !df$rbcL_outlier
seqs <- df$rbcL_CDS[keep]
seq_lengths <- nchar(seqs)
hist(seq_lengths, main="Sequence lengths after removing outliers")

gc_counts <- letterFrequency(DNAStringSet(seqs), letters=c("G","C"))
gc_content <- rowSums(gc_counts) / seq_lengths
hist(gc_content, main="GC content after removing outliers")

df_gc <- data.table(ID=df$ID[keep], gc_seq=gc_content)
merged_gc <- merge(df_gc, X_clean, by.x="ID", by.y="sample_id")

y <- merged_gc$gc_seq
cors <- sapply(merged_gc[, !c("ID","gc_seq"), with=FALSE],
               \(x) cor(x, y, use="pairwise.complete.obs"))
hist(cors, main="rbcL embedding correlations with rbcl_CDS GC")
best <- names(cors)[which.max(cors)]
plot(merged_gc[[best]], y,
     xlab=best, ylab="rbcL GC content",
     main=paste("Best correlated embedding:", best))

best <- names(cors)[which.max(cors)]
plot(merged_gc[[best]], y, xlab=best, ylab="rbcL GC content")

library(glmnet)
library(doParallel)

cl <- makeCluster(10)
registerDoParallel(cl)

merged_clean <- merge(df[, .(ID, pheno_Topt_site_p50, Taxonomy)], X_clean, by.x="ID", by.y="sample_id")

# Hold out random 20% for test
set.seed(42)
test_idx <- sample(seq_len(nrow(merged_clean)), size = floor(0.2*nrow(merged_clean)))
train_idx <- setdiff(seq_len(nrow(merged_clean)), test_idx)

X_train <- as.matrix(merged_clean[train_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_train <- merged_clean$pheno_Topt_site_p50[train_idx]

X_test <- as.matrix(merged_clean[test_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_test <- merged_clean$pheno_Topt_site_p50[test_idx]

fit <- cv.glmnet(X_train, y_train, alpha=1, standardize=TRUE, parallel=TRUE)

pred_train <- predict(fit, X_train, s="lambda.min")
pred_test <- predict(fit, X_test, s="lambda.min")

par(mfrow=c(1,2))
plot(pred_train, y_train, main="Train", xlab="Predicted", ylab="Observed")
abline(0,1,col="red")
plot(pred_test, y_test, main="Test", xlab="Predicted", ylab="Observed")
abline(0,1,col="red")
cor(pred_train, y_train)
cor(pred_test, y_test)

train_idx <- which(!grepl("Poales", merged_clean$Taxonomy))
test_idx  <- which(grepl("Poales", merged_clean$Taxonomy))

X_train <- as.matrix(merged_clean[train_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_train <- merged_clean$pheno_Topt_site_p50[train_idx]

X_test <- as.matrix(merged_clean[test_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_test <- merged_clean$pheno_Topt_site_p50[test_idx]

fit <- cv.glmnet(X_train, y_train, alpha=1, standardize=TRUE, parallel=TRUE)
plot(fit)
pred_train <- predict(fit, X_train, s="lambda.min")
pred_test  <- predict(fit, X_test, s="lambda.min")

par(mfrow=c(1,2))
plot(pred_train, y_train, main="Train (non-Poales)", xlab="Predicted", ylab="Observed")
abline(0,1,col="red")
plot(pred_test, y_test, main="Poales Holdout", xlab="Predicted", ylab="Observed")
abline(0,1,col="red")

cat("Train Pearson:", cor(pred_train, y_train), "\n")
cat("Poales Holdout Pearson:", cor(pred_test, y_test), "\n")

pca <- readRDS("data/tmp/pca.rds")
dim(pca$x)

merged_clean <- merge(df[, .(ID, pheno_Topt_site_p50, Taxonomy)], X_clean, by.x="ID", by.y="sample_id")
train_idx <- which(!grepl("Poales", merged_clean$Taxonomy))
test_idx  <- which(grepl("Poales", merged_clean$Taxonomy))

# ---- Model 1: embeddings ----
X_train_emb <- as.matrix(merged_clean[train_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_train <- merged_clean$pheno_Topt_site_p50[train_idx]
X_test_emb <- as.matrix(merged_clean[test_idx, !c("ID","pheno_Topt_site_p50","Taxonomy"), with=FALSE])
y_test <- merged_clean$pheno_Topt_site_p50[test_idx]

fit_emb <- cv.glmnet(X_train_emb, y_train, alpha=1, standardize=TRUE, parallel=TRUE)
pred_train_emb <- predict(fit_emb, X_train_emb, s="lambda.min")
pred_test_emb  <- predict(fit_emb, X_test_emb, s="lambda.min")

# ---- Model 2: first 10 PCs ----
pca <- readRDS("data/tmp/pca.rds")
pca_df <- data.table(ID = X_clean$sample_id, pca$x[,1:10])
merged_pcs <- merge(df[, .(ID, pheno_Topt_site_p50, Taxonomy)], pca_df, by="ID")
X_train_pc <- as.matrix(merged_pcs[train_idx, paste0("PC",1:10), with=FALSE])
X_test_pc  <- as.matrix(merged_pcs[test_idx, paste0("PC",1:10), with=FALSE])
y_train_pc <- merged_pcs$pheno_Topt_site_p50[train_idx]
y_test_pc  <- merged_pcs$pheno_Topt_site_p50[test_idx]

fit_pc <- cv.glmnet(X_train_pc, y_train_pc, alpha=1, standardize=TRUE, parallel=TRUE)
pred_train_pc <- predict(fit_pc, X_train_pc, s="lambda.min")
pred_test_pc  <- predict(fit_pc, X_test_pc, s="lambda.min")

# ---- Compare ----
cat("Embeddings: Train r =", cor(pred_train_emb, y_train), " Test r =", cor(pred_test_emb, y_test), "\n")
cat("PCs: Train r =", cor(pred_train_pc, y_train_pc), " Test r =", cor(pred_test_pc, y_test_pc), "\n")







library(pheatmap)
mat <- dcast(res, embedding ~ pheno, value.var="cor")
embeddings <- mat$embedding
mat <- as.matrix(mat[,-1])
rownames(mat) <- embeddings

bio_names <- c(
  "bio_1"="Annual Mean Temperature",
  "bio_2"="Mean Diurnal Range",
  "bio_3"="Isothermality (BIO2/BIO7 ×100)",
  "bio_4"="Temperature Seasonality (SD×100)",
  "bio_5"="Max Temp of Warmest Month",
  "bio_6"="Min Temp of Coldest Month",
  "bio_7"="Temperature Annual Range (BIO5-BIO6)",
  "bio_8"="Mean Temp of Wettest Quarter",
  "bio_9"="Mean Temp of Driest Quarter",
  "bio_10"="Mean Temp of Warmest Quarter",
  "bio_11"="Mean Temp of Coldest Quarter",
  "bio_12"="Annual Precipitation",
  "bio_13"="Precipitation of Wettest Month",
  "bio_14"="Precipitation of Driest Month",
  "bio_15"="Precipitation Seasonality (CV)",
  "bio_16"="Precipitation of Wettest Quarter",
  "bio_17"="Precipitation of Driest Quarter",
  "bio_18"="Precipitation of Warmest Quarter",
  "bio_19"="Precipitation of Coldest Quarter"
)

labels <- sub(".*(bio_\\d+).*","\\1",colnames(df))
labels <- ifelse(grepl("bio_",labels), bio_names[labels], labels)

col_labels <- sub(".*(bio_\\d+)_p(\\d+).*","\\1_p\\2",colnames(mat))
col_labels <- ifelse(grepl("bio_",col_labels),
                     paste0(bio_names[sub("(bio_\\d+)_.*","\\1",col_labels)],
                            " (p",sub(".*_p(\\d+)","\\1",col_labels),")"),
                     col_labels)

pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, labels_col=col_labels)





par(mfrow=c(1,1))
hist(cors, main="correlations of cleaned embeddings with T_opt_site")
high_cor <- cors[abs(cors) > 0.7]
high_cor



ggplot(scores, aes(x = PC1, y = PC2, color = pheno)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(
    title = "PCA of rbcL Embeddings Colored by Topt", 
    x=paste0("PC1 ",pvar[1] ,"%"),
    y=paste0("PC2 ",pvar[2] ,"%"),color = "Topt"
  ) +
  theme_minimal()

# Print outlier coordinates for inspection
outlier_indices <- which(is_outlier)
if(length(outlier_indices) > 0) {
  cat("Outlier positions in MDS space:\n")
  for(i in outlier_indices) {
    cat(sprintf("%s: (%.3f, %.3f)\n", 
                sample_ids[i], 
                mds_result[i,1], 
                mds_result[i,2]))
  }
}
