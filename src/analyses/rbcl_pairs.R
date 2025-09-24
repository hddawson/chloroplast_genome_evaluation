library(Biostrings)
library(data.table)
library(ggplot2)
library(stringr)
library(pROC)
library(ggplot2)
library(glmnet)
library(doParallel)
# Read sequences
seqs <- readAAStringSet("data/tmp/alignedGenes/rbcL_AA_aligned.fasta")
seq_names <- names(seqs)
n <- length(seqs)

# Compute pairwise distances (proportion of mismatches)
dist_mat <- as.matrix(stringDist(seqs, method="hamming") / width(seqs[1]))
rownames(dist_mat) <- colnames(dist_mat) <- seq_names

hist(dist_mat)
# Convert to long format for plotting
pairs_idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
dist_df <- data.frame(
  seq1 = seq_names[pairs_idx[,1]],
  seq2 = seq_names[pairs_idx[,2]],
  AA_distance = dist_mat[pairs_idx],
  stringsAsFactors = FALSE
)
dist_df$seq1 <- str_split_i(dist_df$seq1, "\\|", 1)
dist_df$seq2 <- str_split_i(dist_df$seq2, "\\|", 1)

setDT(dist_df)
setDT(df)

# Merge phenotype data
df <- as.data.table(read_parquet("data/processed_data.parquet"))
df <- df[!ID %in% c("MZ475300.1","MZ160994.1","MT876494.1","JN032131.1","MK645903.1","MW817015.1","MW829763.1","MH325132.1")]

merged <- dist_df[
  df, on = c(seq1 = "ID"), nomatch = 0
][
  df, on = c(seq2 = "ID"), nomatch = 0,
  .(seq1, seq2, AA_distance, 
    Organism1 = i.Organism, Topt1 = i.pheno_Topt_site_p50,
    Organism2 = x.Organism, Topt2 = x.pheno_Topt_site_p50)
]

merged <- merged[AA_distance != 0]

# Compute phenotype distance
merged[, Phenotype_distance := Topt1 - Topt2]
hist(merged$Phenotype_distance)
hist(merged$AA_distance, main="rbcL hamming distance")
summary(merged$Phenotype_distance)

# Optional: focus on minimal AA distance pairs with reasonable phenotype difference
aa_cutoff <- quantile(merged$AA_distance, 0.01)
y_cutoff  <- quantile(abs(merged$Phenotype_distance), 0.5)

candidate_pairs <- merged[AA_distance < aa_cutoff & abs(Phenotype_distance) > y_cutoff]

plot(candidate_pairs$AA_distance,candidate_pairs$Phenotype_distance)
# Plot all pairs

saveRDS(candidate_pairs, "data/tmp/candidatePairs.rds")

library(Biostrings)
library(data.table)
library(ggplot2)
library(stringr)
library(pROC)
library(ggplot2)
library(glmnet)
library(doParallel)
candidate_pairs <- readRDS("data/tmp/candidatePairs.rds")

embeds <- as.data.table(read_parquet("data/rbcL_embeddings.parquet"))
embeds[1:3,1:3]

setDT(embeds)
setkey(embeds, sample_id)

# Join embeddings for seq1
data <- embeds[candidate_pairs, on = c(sample_id = "seq1")]
setnames(data, old = names(embeds)[-1], new = paste0("E1_", names(embeds)[-1]))

# Join embeddings for seq2
data <- embeds[data, on = c(sample_id = "seq2")]
setnames(data, old = names(embeds)[-1], new = paste0("E2_", names(embeds)[-1]))

# Subtract embeddings
embed_cols <- grep("^E1_", names(data), value = TRUE)
diff_cols <- paste0("diff_", seq_along(embed_cols))
data[, (diff_cols) := lapply(seq_along(embed_cols), function(i) get(embed_cols[i]) - get(sub("E1_", "E2_", embed_cols[i])))]

# Create logistic response: sign of phenotype difference
data[, y_sign := as.integer(Phenotype_distance > 0)]

table(data$y_sign)

set.seed(42)  # for reproducibility
n <- nrow(data)
train_idx <- sample(n, size = 0.8 * n)

# Training set
X_train <- as.matrix(data[train_idx, ..diff_cols])
y_train <- data$y_sign[train_idx]

# Test set
X_test <- as.matrix(data[-train_idx, ..diff_cols])
y_test <- data$y_sign[-train_idx]

# Check balance
table(y_train)
table(y_test)

fit_simple <- glmnet(X_train, y_train, family="binomial", alpha=1, lambda=0.01,trace.it=1)

# Predictions
pred_train <- predict(fit_simple, X_train, type="response")
pred_test  <- predict(fit_simple, X_test, type="response")

cor(pred_train,y_train)
cor(pred_test,y_test)

fit <- cv.glmnet(X_train, y_train, family="binomial", alpha=0,
                 parallel=F,trace.it = 1)

saveRDS(fit, "data/tmp/fit.rds")
cat("done!!!")
plot(fit)
title("Lasso Logistic Regression: CV Performance", line = 2.5)


coef_min <- coef(fit, s = "lambda.min")
nonzero_coef <- coef_min[coef_min != 0]
print(nonzero_coef)

# Predictions on training and test
pred_train <- predict(fit, X_train, s = "lambda.min", type = "response")
pred_test  <- predict(fit, X_test, s = "lambda.min", type = "response")

# Convert to class labels using 0.5 threshold
pred_class_train <- as.integer(pred_train > 0.5)
pred_class_test  <- as.integer(pred_test  > 0.5)

# Confusion matrix
table(train = pred_class_train, truth = y_train)
table(test  = pred_class_test,  truth = y_test)

# ROC and AUC
roc_train <- roc(y_train, pred_train)
roc_test  <- roc(y_test,  pred_test)

plot(roc_train, col="blue", lwd=2, main="ROC Curves (20% Holdout)")
lines(roc_test, col="red", lwd=2)
abline(a=0,b=1,lty=2,col="gray")

legend("bottomright", legend=c(
  paste0("Train (AUC=", round(auc(roc_train),3),")"),
  paste0("Test (AUC=", round(auc(roc_test),3),")")),
  col=c("blue","red"), lwd=2, bg="white")

# Calibration plot (predicted vs observed)
calib_df <- data.frame(pred = pred_test, obs = y_test)

ggplot(calib_df, aes(x = pred, y = obs)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  labs(x = "Predicted probability", y = "Observed outcome", title = "Calibration Plot: Test Data") +
  theme_minimal()