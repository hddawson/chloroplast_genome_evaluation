library(data.table)
library(ggplot2)
library(pROC)
library(ggplot2)
library(glmnet)
library(doParallel)
library(arrow)
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

cl <- makeCluster(10)
registerDoParallel(cl)

fit <- cv.glmnet(X_train, y_train, family="binomial", alpha=1,
                 parallel=FALSE,trace.it = 1)

saveRDS(fit, "data/tmp/fit.rds")
fit <- readRDS("data/tmp/fit.rds")
cat("done!!!")
plot(fit)
title("Lasso Logistic Regression: CV Performance", line = 2.5)
fit
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

plot(roc_train, col = "blue", main = "ROC Curves", lwd = 2)
lines(roc_test, col = "red", lwd = 2)
legend("bottomright", legend = c("Train", "Test"), col = c("blue","red"), lwd = 2)

auc(roc_train)
auc(roc_test)

# Calibration plot (predicted vs observed)
calib_df <- data.frame(pred = pred_test, obs = y_test)

ggplot(calib_df, aes(x = pred, y = obs)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  labs(x = "Predicted probability", y = "Observed outcome", title = "Calibration Plot: Test Data") +
  theme_minimal()