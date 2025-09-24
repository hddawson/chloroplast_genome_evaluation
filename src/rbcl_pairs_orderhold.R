library(Biostrings)
library(data.table)
library(ggplot2)
library(stringr)
library(pROC)
library(glmnet)
library(doParallel)

# Read sequences
seqs <- readAAStringSet("data/tmp/alignedGenes/rbcL_AA_aligned.fasta")
seq_names <- names(seqs)
n <- length(seqs)

# Compute pairwise distances (proportion of mismatches)
dist_mat <- as.matrix(stringDist(seqs, method="hamming") / width(seqs[1]))
rownames(dist_mat) <- colnames(dist_mat) <- seq_names

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

# Read phenotype data
df <- as.data.table(read_parquet("data/processed_data.parquet"))
df <- df[!ID %in% c("MZ475300.1","MZ160994.1","MT876494.1","JN032131.1","MK645903.1","MW817015.1","MW829763.1","MH325132.1")]

# Identify Poales samples
poales_ids <- df[grep("Poales", Taxonomy), ID]
cat("Found", length(poales_ids), "Poales samples\n")

# Merge phenotype data
merged <- dist_df[
  df, on = c(seq1 = "ID"), nomatch = 0
][
  df, on = c(seq2 = "ID"), nomatch = 0,
  .(seq1, seq2, AA_distance, 
    Organism1 = i.Organism, Topt1 = i.pheno_Topt_site_p50,
    Organism2 = x.Organism, Topt2 = x.pheno_Topt_site_p50,
    Tax1 = i.Taxonomy, Tax2 = x.Taxonomy)
]

merged <- merged[AA_distance != 0]

# Compute phenotype distance
merged[, Phenotype_distance := Topt1 - Topt2]

# Filter candidate pairs
aa_cutoff <- quantile(merged$AA_distance, 0.01)
y_cutoff  <- quantile(abs(merged$Phenotype_distance), 0.5)

candidate_pairs <- merged[AA_distance < aa_cutoff & abs(Phenotype_distance) > y_cutoff]

# Read embeddings
embeds <- as.data.table(read_parquet("data/rbcL_embeddings.parquet"))
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

# Create logistic response
data[, y_sign := as.integer(Phenotype_distance > 0)]
table(data$y_sign)
# Phylogenetic cross-validation: train on non-Poales, test on Poales pairs
data[, is_poales_pair := (sample_id %in% poales_ids) | (i.sample_id %in% poales_ids)]
table(data$is_poales_pair)
train_data <- data[is_poales_pair == FALSE]
test_data <- data[is_poales_pair == TRUE]

cat("Training samples:", nrow(train_data), "\n")
cat("Test samples:", nrow(test_data), "\n")

# Training matrices
X_train <- as.matrix(train_data[, ..diff_cols])
y_train <- train_data$y_sign

# Test matrices  
X_test <- as.matrix(test_data[, ..diff_cols])
y_test <- test_data$y_sign

# Check balance
cat("Training balance:", table(y_train), "\n")
cat("Test balance:", table(y_test), "\n")

# Fit model with cross-validation
set.seed(42)
fit <- glmnet(X_train, y_train, family="binomial", alpha=1, lambda=0.01,trace.it=1)

# Extract coefficients
coef_min <- coef(fit, s = "lambda.min")
nonzero_coef <- coef_min[coef_min != 0]
cat("Non-zero coefficients:", length(nonzero_coef) - 1, "\n")

# Predictions
pred_train <- predict(fit, X_train, s = "lambda.min", type = "response")
pred_test  <- predict(fit, X_test, s = "lambda.min", type = "response")

# Performance metrics
pred_class_train <- as.integer(pred_train > 0.5)
pred_class_test  <- as.integer(pred_test > 0.5)

# Confusion matrices
cat("\nTraining confusion matrix:\n")
print(table(predicted = pred_class_train, truth = y_train))

cat("\nTest confusion matrix:\n") 
print(table(predicted = pred_class_test, truth = y_test))

# ROC analysis
roc_train <- roc(y_train, pred_train)
roc_test  <- roc(y_test, pred_test)

# Plot ROC curves
plot(roc_train, col="blue", lwd=2, main="ROC Curves: Phylogenetic Cross-Validation")
lines(roc_test, col="red", lwd=2)
abline(a=0, b=1, lty=2, col="gray")

legend("bottomright", legend=c(
  paste0("Train (non-Poales, AUC=", round(auc(roc_train), 3), ")"),
  paste0("Test (Poales, AUC=", round(auc(roc_test), 3), ")")),
  col=c("blue", "red"), lwd=2, bg="white")

# Calibration plot
calib_df <- data.frame(pred = as.vector(pred_test), obs = y_test)

p_calib <- ggplot(calib_df, aes(x = pred, y = obs)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Predicted probability", y = "Observed outcome", 
       title = "Calibration Plot: Poales Test Data") +
  theme_minimal()

print(p_calib)

# Summary statistics
cat("\nSummary:\n")
cat("Train AUC:", round(auc(roc_train), 3), "\n")
cat("Test AUC:", round(auc(roc_test), 3), "\n")
cat("Train accuracy:", round(mean(pred_class_train == y_train), 3), "\n")
cat("Test accuracy:", round(mean(pred_class_test == y_test), 3), "\n")

