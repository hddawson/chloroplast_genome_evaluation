library(data.table)
library(arrow)
library(ape)

rmse <- function(a,b) sqrt(mean((a - b)^2))

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
table(embeds$ManualOutlier)
dim(embeds)
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

df <- read_parquet("data/processed_data.parquet")
setDT(df)

### --- PREVENT CONTAMINATION: Split data first ---
poaceae_ids <- df[grepl("Poaceae", Taxonomy), ID]

# Create training set (exclude Poales for variable selection)
train_df <- df[!ID %in% poaceae_ids]
train_embeds <- clean_embeds[!ID %in% poaceae_ids]

# Keep test set separate
test_df <- df[ID %in% poaceae_ids]

### --- Variable selection on TRAINING data only ---
setDT(train_embeds)
pheno_col <- "pheno_wc2.1_2.5m_bio_13_p50"
n_top <- 1  # Fixed at 1
genes <- unique(train_embeds$Gene)

# Pre-merge TRAINING data only
merged_data <- merge(train_embeds[, c("ID", "Gene", embed_cols), with = FALSE], 
                     train_df[, .(ID, Order, pheno = get(pheno_col))], 
                     by = "ID")

# Pre-filter genes with sufficient samples
gene_counts <- merged_data[, .N, by = Gene]
valid_genes <- gene_counts[N >= 5]$Gene
summary(gene_counts$N)

compute_gene_cors <- function(gene_data, embed_cols) {
  embed_matrix <- as.matrix(gene_data[, embed_cols, with = FALSE])
  cors <- cor(embed_matrix, gene_data$pheno, use = "complete.obs")[, 1]
  names(cors) <- embed_cols
  return(cors)
}

# Vectorized correlation computation (unchanged)
compute_gene_cors <- function(gene_data, embed_cols) {
  order_groups <- split(gene_data, gene_data$Order)
  
  cors_by_order <- lapply(order_groups, function(group) {
    embed_matrix <- as.matrix(group[, embed_cols, with = FALSE])
    cors <- cor(embed_matrix, group$pheno, use = "complete.obs")[, 1]
    return(cors)
  })
  
  cors_matrix <- do.call(rbind, cors_by_order)
  mean_cors <- colMeans(cors_matrix, na.rm = TRUE)
  names(mean_cors) <- embed_cols
  
  return(mean_cors)
}

# Variable selection on training data
sel_list <- vector("list", length(valid_genes))
names(sel_list) <- valid_genes

for (i in seq_along(valid_genes)) {
  gene <- valid_genes[i]
  cat(gene, "\n")
  
  gdt <- merged_data[Gene == gene]
  cors_by_order <- compute_gene_cors(gdt, embed_cols)
  
  # Select top 1 dimension
  available <- names(cors_by_order)[!is.na(cors_by_order)]
  k <- min(n_top, length(available))
  top_dims <- names(sort(abs(cors_by_order[available]), decreasing = TRUE))[1:k]
  
  sel <- gdt[, c("ID", top_dims), with = FALSE]
  sel <- unique(sel)
  setnames(sel, old = top_dims, new = paste0(gene, "__", top_dims))
  sel_list[[gene]] <- sel
}

# Remove NULL entries
sel_list <- sel_list[!sapply(sel_list, is.null)]

# Merge selected features across all training IDs
combined_train <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)
combined_train <- merge(combined_train, train_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)

sum(is.na(combined_train))
# Median impute missing predictors
pred_cols <- setdiff(names(combined_train), c("ID", "pheno"))
for (col in pred_cols) {
  med <- median(combined_train[[col]], na.rm = TRUE)
  combined_train[is.na(get(col)), (col) := med]
}

cat("Training samples:", nrow(combined_train), "\n")
cat("Number of predictors:", length(pred_cols), "\n")

### --- Model fitting and validation ---
pred_cols <- grep("embedding", colnames(combined_train), value=TRUE)
formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))

# Fit on training data
fit <- lm(as.formula(formula_str), data=combined_train, na.action = na.exclude)
summary(fit)
plot(fit$fitted.values, combined_train$pheno)

test_combined <- merge(clean_embeds[ID %in% poaceae_ids, c("ID", "Gene", embed_cols), with = FALSE],
                       test_df[, .(ID, pheno = get(pheno_col))], by = "ID")

# Reshape test data to match training format (one row per ID)
test_wide_list <- list()
for (gene in unique(test_combined$Gene)) {
  gene_data <- test_combined[Gene == gene, c("ID", embed_cols), with = FALSE]
  if (nrow(gene_data) > 0) {
    # Rename columns to match training format
    setnames(gene_data, old = embed_cols, new = paste0(gene, "__", embed_cols))
    test_wide_list[[gene]] <- gene_data
  }
}

# Merge all genes for test data
combined_test <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), test_wide_list)
combined_test <- merge(combined_test, test_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)
sum(is.na(combined_test[,..pred_cols]))
dim(combined_test[,..pred_cols])
# Median impute missing predictors (using test medians!)
test_pred_cols <- setdiff(names(combined_test), c("ID", "pheno"))
for (col in test_pred_cols) {
  if (col %in% pred_cols) {  # Only impute columns that exist in training
    med <- median(combined_test[[col]], na.rm = TRUE)
    combined_test[is.na(get(col)), (col) := med]
  }
}

### --- Model fitting and validation ---
pred_cols <- grep("embedding", colnames(combined_train), value=TRUE)
formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))

# Fit on training data (selected features only)
fit <- lm(as.formula(formula_str), data=combined_train)
summary(fit)

# Predict on test data (model will use available features, ignore missing ones)
pred_test <- predict(fit, newdata = combined_test)

# Validation plot
plot(pred_test, combined_test$pheno,
     main="Predicted vs. observed growth temperature across grasses",
     xlab="Predicted",
     ylab="Observed",col="black")
#abline(line(pred_test, combined_test$pheno)) Predicted vs observed temperature tolerance across grasses (Poales excluded from training)
abline(a=0, b=1, col="red")
sum(is.na(pred_test))
corrr <- cor(pred_test[!is.na(pred_test)], combined_test[!is.na(pred_test)]$pheno,method="pearson")
rmseeee <- rmse(pred_test[!is.na(pred_test)], combined_test[!is.na(pred_test)]$pheno)
text(24, 3,
     paste0("Pearson=", round(corrr,3)), col="red")
text(24, 1,
     paste0("RMSE=", round(rmseeee, 3)), col="red")

LSD::heatscatterpoints(pred_test, combined_test$pheno)

cat("Test correlation (Spearman):", cor(pred_test, combined_test$pheno, method="spearman"), "\n")
cat("Test RMSE:", rmse(pred_test, combined_test$pheno), "\n")

library(ggplot2)
coefs <- summary(fit)$coefficients
df <- data.frame(term=rownames(coefs), Estimate=coefs[,1], pval=coefs[,4])
df <- subset(df, pval < 0.05 & term != "(Intercept)")
ci <- confint(fit)[rownames(df),]
df$lower <- ci[,1]; df$upper <- ci[,2]

df$gene <- sub("__embedding.*","",df$term)
df$prefix <- substr(df$gene,1,3)
df$prefix <- ifelse(df$prefix %in% c("psb","psa","rbc","atp","ndh","pet","rpo","rps"), df$prefix, "other")
table(df$prefix)
plot(table(df$prefix), main="frequency of ")
# Option 1: simple barplot
ggplot(df, aes(x=reorder(term, Estimate), y=Estimate, fill=Estimate > 0)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(x="", y="Effect size", title="Significant predictor effects")

# Option 2: forest plot with CI
ci <- confint(fit)[rownames(df),]
df$lower <- ci[,1]; df$upper <- ci[,2]
ggplot(df, aes(x=term, y=Estimate, ymin=lower, ymax=upper, color=Estimate>0)) +
  geom_pointrange() + coord_flip() + theme_minimal() +
  labs(x="", y="Estimate (95% CI)", title="Significant predictors")

# Option 3: group by gene (strip plot)
df$gene <- sub("__embedding.*","",df$term)
ggplot(df, aes(x=gene, y=Estimate, color=gene)) +
  geom_jitter(width=0.2, height=0, size=3) + theme_minimal() +
  labs(x="Gene", y="Effect size", title="Significant embeddings per gene")

ggplot(df, aes(x=reorder(term, Estimate), y=Estimate, ymin=lower, ymax=upper, color=prefix)) +
  geom_pointrange() + coord_flip() + theme_minimal() +
  labs(x="", y="Estimate (95% CI)", title="Significant predictors by complex")

ggplot(df, aes(x=Estimate, y=term, color=prefix)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  theme_minimal() + labs(y="", x="Effect size", title="Predictor effects with complexes highlighted")

ggplot(df, aes(x=prefix, y=Estimate, color=prefix)) +
  geom_jitter(width=0.2, height=0, alpha=0.7, size=2) +
  geom_boxplot(outlier.shape=NA, alpha=0.2) +
  theme_minimal() + labs(x="Complex", y="Effect size", title="Distribution of significant effects per complex")


