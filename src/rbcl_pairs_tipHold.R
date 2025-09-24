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

# Merge phenotype data
merged <- dist_df[
  df, on = c(seq1 = "ID"), nomatch = 0
][
  df, on = c(seq2 = "ID"), nomatch = 0,
  .(seq1, seq2, AA_distance, 
    Organism1 = i.Organism, Topt1 = i.pheno_Topt_site_p50,
    Organism2 = x.Organism, Topt2 = x.pheno_Topt_site_p50,
    Order1 = i.Order, Order2 = x.Order)
]

merged <- merged[AA_distance != 0]

# Compute phenotype distance
merged[, Phenotype_distance := Topt1 - Topt2]

# Filter pairs within each order
within_order_pairs <- merged[Order1 == Order2]
pheno_cutoff <- quantile(abs(within_order_pairs$Phenotype_distance), 0.5)

# Select top 3 pairs per order based on criteria
selected_pairs <- within_order_pairs[abs(Phenotype_distance) > pheno_cutoff][
  order(AA_distance), .SD[1:min(1, .N)], by = Order1]

# Get all species involved in selected pairs
holdout_species <- unique(c(selected_pairs$seq1, selected_pairs$seq2))

cat("Selected", nrow(selected_pairs), "pairs from", length(unique(selected_pairs$Order1)), "orders\n")
cat("Holdout species:", length(holdout_species), "\n")

# Filter candidate pairs for modeling (broader criteria)
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

# Define pair types for phylogenetic cross-validation
data[, is_holdout_seq1 := sample_id %in% holdout_species]
data[, is_holdout_seq2 := i.sample_id %in% holdout_species]

data[, pair_type := ifelse(is_holdout_seq1 & is_holdout_seq2, "holdout_holdout",
                           ifelse(is_holdout_seq1 | is_holdout_seq2, "cross_species", "training"))]
table(data$pair_type)
# Split data
train_data <- data[pair_type == "training"]
test_cross_species <- data[pair_type == "cross_species"] 
test_holdout <- data[pair_type == "holdout_holdout"]

cat("Training pairs (no holdout species):", nrow(train_data), "\n")
cat("Cross-species test pairs (training-holdout):", nrow(test_cross_species), "\n")
cat("Holdout test pairs (holdout-holdout):", nrow(test_holdout), "\n")

# Training matrices
X_train <- as.matrix(train_data[, ..diff_cols])
y_train <- train_data$y_sign

# Test matrices
X_test_cross <- as.matrix(test_cross_species[, ..diff_cols])
y_test_cross <- test_cross_species$y_sign

# Check balance
cat("Training balance:", table(y_train), "\n")
cat("Cross-species test balance:", table(y_test_cross), "\n")

# Fit model with cross-validation
set.seed(42)
fit <- glmnet(X_train, y_train, family="binomial", alpha=1, lambda=0.01,trace.it=1)

# Extract coefficients
coef_min <- coef(fit, s = "lambda.min")
nonzero_coef <- coef_min[coef_min != 0]
cat("Non-zero coefficients:", length(nonzero_coef) - 1, "\n")

# Predictions
pred_train <- predict(fit, X_train, s = "lambda.min", type = "response")
pred_test_cross <- predict(fit, X_test_cross, s = "lambda.min", type = "response")

# Performance metrics
pred_class_train <- as.integer(pred_train > 0.5)
pred_class_cross <- as.integer(pred_test_cross > 0.5)

# Confusion matrices
cat("\nTraining confusion matrix:\n")
print(table(predicted = pred_class_train, truth = y_train))

cat("\nCross-species test confusion matrix:\n") 
print(table(predicted = pred_class_cross, truth = y_test_cross))

# ROC analysis
roc_train <- roc(y_train, pred_train)
roc_cross <- roc(y_test_cross, pred_test_cross)

# Plot ROC curves
plot(roc_train, col="blue", lwd=2, main="ROC Curves: Sister-Pair Validation")
lines(roc_cross, col="red", lwd=2)
abline(a=0, b=1, lty=2, col="gray")

legend("bottomright", legend=c(
  paste0("Train (no holdout species, AUC=", round(auc(roc_train), 3), ")"),
  paste0("Test (cross-species, AUC=", round(auc(roc_cross), 3), ")")),
  col=c("blue", "red"), lwd=2, bg="white")

# Test on holdout-holdout pairs if available
if(nrow(test_holdout) > 0) {
  X_test_holdout <- as.matrix(test_holdout[, ..diff_cols])
  y_test_holdout <- test_holdout$y_sign
  
  pred_test_holdout <- predict(fit, X_test_holdout, s = "lambda.min", type = "response")
  pred_class_holdout <- as.integer(pred_test_holdout > 0.5)
  
  cat("\nHoldout-holdout test balance:", table(y_test_holdout), "\n")
  cat("Holdout-holdout confusion matrix:\n")
  print(table(predicted = pred_class_holdout, truth = y_test_holdout))
  
  roc_holdout <- roc(y_test_holdout, pred_test_holdout)
  lines(roc_holdout, col="green", lwd=2)
  
  legend("bottomright", legend=c(
    paste0("Train (no holdout species, AUC=", round(auc(roc_train), 3), ")"),
    paste0("Test (cross-species, AUC=", round(auc(roc_cross), 3), ")"),
    paste0("Test (holdout-holdout, AUC=", round(auc(roc_holdout), 3), ")")),
    col=c("blue", "red", "green"), lwd=2, bg="white")
  
  cat("Holdout-holdout AUC:", round(auc(roc_holdout), 3), "\n")
  cat("Holdout-holdout accuracy:", round(mean(pred_class_holdout == y_test_holdout), 3), "\n")
}

# Show selected pairs for validation
cat("\nSelected validation pairs by order:\n")
print(selected_pairs[, .(Order1, seq1, seq2, AA_distance, Phenotype_distance)])

# Calibration plot for cross-species pairs
calib_df <- data.frame(pred = as.vector(pred_test_cross), obs = y_test_cross)

p_calib <- ggplot(calib_df, aes(x = pred, y = obs)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Predicted probability", y = "Observed outcome", 
       title = "Calibration Plot: Cross-Species Test Data") +
  theme_minimal()

print(p_calib)

# Summary statistics
cat("\nSummary:\n")
cat("Train AUC:", round(auc(roc_train), 3), "\n")
cat("Cross-species AUC:", round(auc(roc_cross), 3), "\n")
cat("Train accuracy:", round(mean(pred_class_train == y_train), 3), "\n")
cat("Cross-species accuracy:", round(mean(pred_class_cross == y_test_cross), 3), "\n