# Combined Analysis: Phylogenetic sister pairs + embedding differentials
library(data.table)
library(arrow)
library(ape)
library(ggplot2)
library(glmnet)
library(pROC)
# 1. Load and clean data
df <- as.data.table(read_parquet("/workdir/hdd29/chloroplast_genome_evaluation/data/processed_data.parquet"))
X <- as.data.table(read_parquet("/workdir/hdd29/chloroplast_genome_evaluation/data/rbcL_embeddings.parquet"))

# 2. Remove outliers from embeddings
mds_result <- as.data.frame(readRDS("/workdir/hdd29/chloroplast_genome_evaluation/results/mds_result.rds"))
mds_df <- data.frame(
  MDS1 = mds_result[,1],
  MDS2 = mds_result[,2],
  sample_id = X$sample_id
)
d <- sqrt(rowSums(scale(mds_df[,c("MDS1","MDS2")], center=TRUE, scale=FALSE)^2))
z <- scale(d)
mds_df$is_outlier <- abs(z) > 4

# Clean embeddings
embedding_cols <- names(X)[names(X) != "sample_id"]
outlier_samples <- mds_df[mds_df$is_outlier, "sample_id"]
keep_idx <- !X$sample_id %in% outlier_samples
X_clean <- X[keep_idx, ]
sum(is.na(X_clean))

# 3. Load phylogenetic data
tree <- read.tree("tree/T2.raxml.rba.raxml.lastTree.TMP")
spermatophyta_data <- df[grep("Spermatophyta", df$Taxonomy),]

spermatophyta_ids <- spermatophyta_data$ID
ids_in_tree <- intersect(spermatophyta_ids, tree$tip.label)

spermatophyta_tree <- keep.tip(tree, ids_in_tree)
spermatophyta_data <- spermatophyta_data[match(spermatophyta_tree$tip.label, spermatophyta_data$ID), ]
complete_cases <- complete.cases(spermatophyta_data$pheno_Topt_site_p50)
spermatophyta_data <- spermatophyta_data[complete_cases, ]
spermatophyta_tree <- keep.tip(spermatophyta_tree, spermatophyta_data$ID)

# Fix zero branch lengths
zero_branches <- spermatophyta_tree$edge.length == 0
very_small <- spermatophyta_tree$edge.length < 1e-6
spermatophyta_tree$edge.length[zero_branches | very_small] <- 1e-6
cat("Fixed", sum(zero_branches), "zero branches and", sum(very_small), "very small branches\n")

# 4. Find sister pairs (adjacent tips with no overlaps)
find_sister_pairs <- function(tree) {
  pairs <- list()
  used_tips <- character(0)
  
  for(i in 1:tree$Nnode) {
    node <- i + length(tree$tip.label)
    children <- which(tree$edge[,1] == node)
    
    # Check if this node has exactly 2 tip children
    child_nodes <- tree$edge[children, 2]
    tip_children <- child_nodes[child_nodes <= length(tree$tip.label)]
    
    if(length(tip_children) == 2) {
      tip_names <- tree$tip.label[tip_children]
      # Only add if neither tip is already used
      if(!any(tip_names %in% used_tips)) {
        pairs[[length(pairs) + 1]] <- tip_names
        used_tips <- c(used_tips, tip_names)
      }
    }
  }
  
  return(do.call(rbind, lapply(pairs, function(p) data.frame(seq1=p[1], seq2=p[2]))))
}

sister_pairs <- find_sister_pairs(spermatophyta_tree)
cat("Found", nrow(sister_pairs), "non-overlapping sister pairs\n")

# 5. Calculate phenotype distances for sister pairs
sister_pairs_dt <- as.data.table(sister_pairs)
pheno_data <- spermatophyta_data[, .(ID, pheno_Topt_site_p50)]
setnames(pheno_data, "ID", "sample_id")

# Add phenotypes for both sequences
sister_pairs_dt <- pheno_data[sister_pairs_dt, on = c(sample_id = "seq1")]
setnames(sister_pairs_dt, "pheno_Topt_site_p50", "pheno1")
sister_pairs_dt <- pheno_data[sister_pairs_dt, on = c(sample_id = "seq2")]  
setnames(sister_pairs_dt, "pheno_Topt_site_p50", "pheno2")

# Calculate phenotype distance
sister_pairs_dt[, Phenotype_distance := pheno1 - pheno2]
sister_pairs_dt <- sister_pairs_dt[complete.cases(sister_pairs_dt)]
hist(sister_pairs_dt$Phenotype_distance * 0.01, main="Adjacent tip pair temp difference (Topt)")
sister_pairs_dt <- sister_pairs_dt[
  which(abs(sister_pairs_dt$Phenotype_distance* 0.01) > 3), 
  ]
#n=1710 if > 2
#about 1227 species left if > 3
#
# 6. Add embedding differences
setkey(X_clean, sample_id)


# 1 C 938
# 2 C 1571
# 3 C 2056
nrow(sister_pairs_dt[
  which(abs(sister_pairs_dt$Phenotype_distance* 0.01) > 3), 
])


# Join embeddings for seq1
data <- X_clean[sister_pairs_dt, on = c(sample_id = "sample_id")]
embed_cols_orig <- names(X_clean)[names(X_clean) != "sample_id"]
setnames(data, old = embed_cols_orig, new = paste0("E1_", embed_cols_orig))

# Join embeddings for seq2  
data <- X_clean[data, on = c(sample_id = "i.sample_id")]
setnames(data, old = embed_cols_orig, new = paste0("E2_", embed_cols_orig))

# Calculate embedding differences
embed_cols <- grep("^E1_", names(data), value = TRUE)
diff_cols <- paste0("diff_", seq_along(embed_cols))
data[, (diff_cols) := lapply(seq_along(embed_cols), function(i) {
  get(embed_cols[i]) - get(sub("E1_", "E2_", embed_cols[i]))
})]

sum(is.na(data[,..diff_cols]))

cat("Pairs with missing embeddings:\n")
missing_E1 <- rowSums(is.na(data[, paste0("E1_", embed_cols_orig), with=FALSE])) > 0
missing_E2 <- rowSums(is.na(data[, paste0("E2_", embed_cols_orig), with=FALSE])) > 0
cat("Missing E1:", sum(missing_E1), "\n")
cat("Missing E2:", sum(missing_E2), "\n")

data_clean <- data[!missing_E1 & !missing_E2]
cat("Clean pairs remaining:", nrow(data_clean), "\n")

# Recalculate embedding differences on clean data
diff_cols <- paste0("diff_", seq_along(embed_cols_orig))
data_clean[, (diff_cols) := lapply(seq_along(embed_cols_orig), function(i) {
  get(paste0("E1_", embed_cols_orig[i])) - get(paste0("E2_", embed_cols_orig[i]))
})]

# Check for NAs
cat("NAs in differences:", sum(is.na(data_clean[,..diff_cols])), "\n")

# Continue with clean data
data <- data_clean

pca <- prcomp(data[,..diff_cols])

pvar <- pca$sdev^2 / sum(pca$sdev^2)
barplot(pvar)
plot(pca$x[,1],pca$x[,2])

# 7. Create binary response variable
data[, y_sign := as.integer(Phenotype_distance > 0)]
cat("Response distribution:\n")
print(table(data$y_sign))

# 8. Fit logistic regression on embedding differences
diff_matrix <- as.matrix(data[, ..diff_cols])
y <- data$y_sign

# Fit model (use subset of embeddings to avoid overfitting)
set.seed(42)

logistic_model <- glm(y ~ diff_matrix, family = binomial())
summary(logistic_model)

# 9. Model diagnostics
pred_prob <- predict(logistic_model, type = "response")
pred_class <- as.integer(pred_prob > 0.5)
accuracy <- mean(pred_class == y)

cat("\nLogistic Regression Results:\n")
cat("Accuracy:", round(accuracy, 3), "\n")
cat("AUC:", round(pROC::auc(y, pred_prob), 3), "\n")

# 10. Visualization
data[, pred_prob := pred_prob]

# Plot phenotype distance vs prediction probability
ggplot(data, aes(x = Phenotype_distance, y = pred_prob)) +
  geom_point(aes(color = factor(y_sign)), alpha = 0.6) +
  geom_smooth(method = "loess") +
  labs(x = "Phenotype Distance (Topt)", 
       y = "Predicted Probability", 
       color = "Actual Sign",
       title = "Sister Pair Analysis: Phenotype vs Embedding Predictions") +
  theme_minimal()

# Summary statistics
cat("\nSister Pairs Summary:\n")
cat("Number of pairs:", nrow(data), "\n")
cat("Mean |phenotype difference|:", round(mean(abs(data$Phenotype_distance)), 3), "\n")
print(head(data[, .(sample_id, i.sample_id, Phenotype_distance, y_sign, pred_prob)]))

set.seed(123)
n_pairs <- nrow(data)
test_idx_random <- sample(n_pairs, round(0.1 * n_pairs))
train_idx_random <- setdiff(1:n_pairs, test_idx_random)

order_info <- df[, .(ID, Order)]

# Join Order for seq1
setnames(order_info, c("ID", "Order"), c("sample_id", "Order1"))
data <- order_info[data, on = c(sample_id = "sample_id")]

# Join Order for seq2  
setnames(order_info, c("sample_id", "Order1"), c("sample_id", "Order2"))
data <- order_info[data, on = c(sample_id = "i.sample_id")]

train_data_random <- data[train_idx_random]
test_data_random <- data[test_idx_random]

# Fit Lasso on random split
X_train_random <- as.matrix(train_data_random[, ..diff_cols])
y_train_random <- train_data_random$y_sign
X_test_random <- as.matrix(test_data_random[, ..diff_cols])
y_test_random <- test_data_random$y_sign

fold_ids <- as.integer(as.factor(train_data_random$Order1))
sum(fold_ids == as.integer(as.factor(train_data_random$Order2)))
nrow(train_data_random) #7 misses .. okay
nfolds <- length(unique(fold_ids))


cv_lasso_random <- cv.glmnet(X_train_random, y_train_random, family="binomial", alpha=1,trace.it=1,
                             nfolds=nfolds,foldid = fold_ids)
saveRDS(cv_lasso_random, "data/tmp/cv_lasso_random.rds")
pred_random <- predict(cv_lasso_random, X_test_random, s="lambda.min", type="response")
acc_random <- mean((pred_random > 0.5) == y_test_random)
auc_random <- pROC::auc(y_test_random, as.vector(pred_random))

data[, poales_pairs := Order1 == "Poales" & Order2 == "Poales"]
data[is.na(poales_pairs), poales_pairs := FALSE]

sum(is.na(data$poales_pairs))
table(data$poales_pairs)
# Split data based on Poales pairs
train_data_poales <- data[poales_pairs == FALSE]
test_data_poales <- data[poales_pairs == TRUE]

# Fit Lasso excluding Poales
X_train_poales <- as.matrix(train_data_poales[, ..diff_cols])
y_train_poales <- train_data_poales$y_sign
X_test_poales <- as.matrix(test_data_poales[, ..diff_cols])
y_test_poales <- test_data_poales$y_sign

fold_ids <- as.integer(as.factor(train_data_poales$Order1))
sum(fold_ids == as.integer(as.factor(train_data_poales$Order2)))
nrow(train_data_poales) #7 misses .. okay
nfolds <- length(unique(fold_ids))

sums <- rowSums(train_data_poales[, ..diff_cols])
a <- as.data.frame(table(sums))
plot(as.numeric(a$sums), a$Freq)
summary(sums)
sum(sums==0)

cv_lasso_poales <- cv.glmnet(X_train_poales, y_train_poales, family="binomial", alpha=1,trace.it=1,
                             nfolds=nfolds,foldid = fold_ids)
pred_poales <- predict(cv_lasso_poales, X_test_poales, s="lambda.min", type="response")
acc_poales <- mean((pred_poales > 0.5) == y_test_poales)
auc_poales <- pROC::auc(y_test_poales, as.vector(pred_poales))

cat("\nCross-validation Results:\n")
cat("Random 10% holdout:\n")
cat("  Train size:", nrow(train_data_random), "Test size:", nrow(test_data_random), "\n")
cat("  Accuracy:", round(acc_random, 3), "AUC:", round(auc_random, 3), "\n")

cat("Poales holdout:\n") 
cat("  Train size:", nrow(train_data_poales), "Test size:", nrow(test_data_poales), "\n")
cat("  Accuracy:", round(acc_poales, 3), "AUC:", round(auc_poales, 3), "\n")

# Compare feature selection
coefs_random <- coef(cv_lasso_random, s="lambda.min")
coefs_poales <- coef(cv_lasso_poales, s="lambda.min")


cat("Non-zero coefficients:\n")
cat("  Random model:", sum(coefs_random != 0) - 1, "\n")  # -1 for intercept
cat("  Poales model:", sum(coefs_poales != 0) - 1, "\n")

# 12. Save important results for later analysis
save_objects <- list(
  # Core data
  sister_pairs_final = data,
  embedding_differences = as.matrix(data[, ..diff_cols]),
  
  # Models
  cv_lasso_random = cv_lasso_random,
  cv_lasso_poales = if(exists("cv_lasso_poales")) cv_lasso_poales else NULL,
  
  # Results
  results_summary = list(
    n_pairs = nrow(data),
    mean_pheno_diff = mean(abs(data$Phenotype_distance)),
    random_accuracy = acc_random,
    random_auc = auc_random,
    poales_accuracy = if(exists("acc_poales")) acc_poales else NA,
    poales_auc = if(exists("auc_poales")) auc_poales else NA
  ),
  
  # Metadata
  embedding_cols = diff_cols,
  tree_info = list(
    n_tips = length(spermatophyta_tree$tip.label),
    n_pairs_found = nrow(sister_pairs)
  )
)
cat("TENFOLD!!")

saveRDS(save_objects, "data/tmp/sister_pairs_analysis_results.rds")
cat("Saved results to: sister_pairs_analysis_results.rds\n")
a <- readRDS("data/tmp/sister_pairs_analysis_results.rds")

zero_diffs <- rowSums(a$embedding_differences==0)

sum(a$embedding_differences==0)
dim(a$embedding_differences)
table(zero_diffs)
plot(as.vector(a$embedding_differences))

par(mfrow=c(1,3))
hist(a$embedding_differences, main="All emb differences")
hist(a$embedding_differences[a$embedding_differences!=0], main="All non-zero")

res <- a$results_summary

zero_map <- (a$embedding_differences == 0)

image(t(zero_map[nrow(zero_map):1, ]), 
      col=c("white","black"), 
      main="Zeros (black) vs Nonzeros (white)")

# Count zeros per column
col_zero_counts <- colSums(a$embedding_differences == 0)

# Which columns are *all* zero?
all_zero_cols <- which(col_zero_counts == nrow(a$embedding_differences))

length(all_zero_cols)   # how many
all_zero_cols[1:20] 

summary(colSums(a$embedding_differences == 0))
summary(rowSums(a$embedding_differences == 0))

