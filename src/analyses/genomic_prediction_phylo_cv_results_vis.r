# Reload and visualize best model results
library(data.table)
library(ggplot2)

res <- readRDS("results/phylo_cv_clade_results.rds")

# Check grid results
print(res$grid_results[1:10])
print(res$best_config)

# To make the scaled plot for best config, need to re-run that config
# Extract best params
best <- res$best_config

# Quick re-run of best config to get predictions
library(ape)
library(arrow)

tree <- read.tree("results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree")
data <- as.data.table(read_parquet("data/processed_data.parquet"))
embeds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds)
embed_cols <- grep("embedding", colnames(embeds), value = TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
pheno <- setNames(data[[pheno_col]], data$ID)
pheno <- pheno[tree$tip.label]
pheno <- pheno[!is.na(pheno)]
tree <- keep.tip(tree, names(pheno))

clean_embeds <- clean_embeds[ID %in% tree$tip.label]
data <- data[ID %in% tree$tip.label]

# Use same seed as main script
set.seed(42)

# Recreate holdout clades (same function as main script)
select_holdout_clades <- function(tree, target_frac = 0.15, min_clade = 10, max_clade = 30) {
  n_tips <- Ntip(tree)
  target_n <- floor(n_tips * target_frac)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  candidates <- list()
  for (node in internal_nodes) {
    tips <- extract.clade(tree, node)$tip.label
    if (length(tips) >= min_clade && length(tips) <= max_clade) {
      candidates[[as.character(node)]] <- tips
    }
  }
  selected <- list()
  used_tips <- character(0)
  candidate_order <- sample(names(candidates))
  for (node in candidate_order) {
    tips <- candidates[[node]]
    if (!any(tips %in% used_tips)) {
      selected[[node]] <- tips
      used_tips <- c(used_tips, tips)
      if (length(used_tips) >= target_n) break
    }
  }
  return(selected)
}

holdout_clades <- select_holdout_clades(tree)
holdout_ids <- unlist(holdout_clades)
train_ids <- setdiff(tree$tip.label, holdout_ids)

# Source the helper functions or redefine them inline
compute_gene_cors <- function(gene_data, embed_cols) {
  order_groups <- split(gene_data, gene_data$Order)
  cors_by_order <- lapply(order_groups, function(group) {
    if (nrow(group) < 3) return(rep(NA, length(embed_cols)))
    embed_matrix <- as.matrix(group[, embed_cols, with = FALSE])
    cors <- cor(embed_matrix, group$pheno, use = "complete.obs")[, 1]
    return(cors)
  })
  cors_matrix <- do.call(rbind, cors_by_order)
  mean_cors <- colMeans(cors_matrix, na.rm = TRUE)
  names(mean_cors) <- embed_cols
  return(mean_cors)
}

# Simplified feature selection for best config
train_df <- data[ID %in% train_ids]
train_embeds <- clean_embeds[ID %in% train_ids]
test_df <- data[ID %in% holdout_ids]
test_embeds <- clean_embeds[ID %in% holdout_ids]

merged <- merge(
  train_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
  train_df[, .(ID, Order, pheno = get(pheno_col))],
  by = "ID"
)

gene_counts <- merged[, .N, by = Gene]
valid_genes <- gene_counts[N >= 5]$Gene

# Get gene correlations
gene_best_cors <- list()
for (gene in valid_genes) {
  gdt <- merged[Gene == gene]
  cors <- compute_gene_cors(gdt, embed_cols)
  available <- names(cors)[!is.na(cors)]
  if (length(available) > 0) gene_best_cors[[gene]] <- max(abs(cors[available]), na.rm = TRUE)
}

gene_ranks <- sort(unlist(gene_best_cors), decreasing = TRUE)
top_gene_names <- names(gene_ranks)[1:min(best$n_top_genes, length(gene_ranks))]

# Select features with best config params
sel_list <- list()
for (gene in valid_genes) {
  gdt <- merged[Gene == gene]
  cors <- compute_gene_cors(gdt, embed_cols)
  available <- names(cors)[!is.na(cors)]
  if (length(available) == 0) next
  
  n_select <- if (gene %in% top_gene_names) min(best$n_dims_per_gene, length(available)) else 1
  
  embed_matrix <- as.matrix(gdt[, embed_cols, with = FALSE])
  pheno_vec <- gdt$pheno
  selected_dims <- character(0)
  residual_pheno <- pheno_vec
  
  for (d in 1:n_select) {
    if (d == 1) {
      current_cors <- cors[available]
    } else {
      current_cors <- cor(embed_matrix[, available, drop = FALSE], residual_pheno, use = "complete.obs")[, 1]
      names(current_cors) <- available
    }
    best_dim <- names(sort(abs(current_cors), decreasing = TRUE))[1]
    selected_dims <- c(selected_dims, best_dim)
    if (d < n_select) {
      fit_resid <- lm(residual_pheno ~ embed_matrix[, best_dim])
      residual_pheno <- residuals(fit_resid)
      available <- setdiff(available, best_dim)
      if (length(available) == 0) break
    }
  }
  
  sel <- unique(gdt[, c("ID", selected_dims), with = FALSE])
  setnames(sel, old = selected_dims, new = paste0(gene, "__", selected_dims))
  sel_list[[gene]] <- sel
}

combined_train <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)
combined_train <- merge(combined_train, train_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)

pred_cols_all <- grep("embedding", names(combined_train), value = TRUE)

# Median impute
for (col in pred_cols_all) {
  med <- median(combined_train[[col]], na.rm = TRUE)
  combined_train[is.na(get(col)), (col) := med]
}

# Filter by best threshold
# Need to recalc feature cors
feature_cors <- data.table(
  gene = sub("__embedding.*", "", pred_cols_all),
  predictor = pred_cols_all
)
feature_cors[, max_cor := gene_best_cors[gene], by = gene]

keep_genes <- names(gene_best_cors)[gene_best_cors >= best$cor_threshold]
keep_cols <- pred_cols_all[sub("__embedding.*", "", pred_cols_all) %in% keep_genes]

message("Best config: ", length(keep_cols), " predictors")

# Fit
formula_best <- paste("pheno ~", paste(keep_cols, collapse = " + "))
fit_best <- lm(as.formula(formula_best), data = combined_train)

# Prepare test
test_merged <- merge(
  test_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
  test_df[, .(ID, pheno = get(pheno_col))],
  by = "ID"
)

test_wide_list <- list()
for (gene in unique(test_merged$Gene)) {
  gene_data <- test_merged[Gene == gene, c("ID", embed_cols), with = FALSE]
  if (nrow(gene_data) > 0) {
    setnames(gene_data, old = embed_cols, new = paste0(gene, "__", embed_cols))
    test_wide_list[[gene]] <- gene_data
  }
}

combined_test <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), test_wide_list)
combined_test <- merge(combined_test, test_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)

for (col in names(combined_test)) {
  if (col %in% pred_cols_all && any(is.na(combined_test[[col]]))) {
    med <- median(combined_test[[col]], na.rm = TRUE)
    combined_test[is.na(get(col)), (col) := med]
  }
}

# Predict
combined_test$pred <- predict(fit_best, newdata = combined_test)

# Add clade info
id_to_clade <- data.table(
  ID = unlist(holdout_clades),
  clade = rep(names(holdout_clades), sapply(holdout_clades, length))
)
combined_test <- merge(combined_test, id_to_clade, by = "ID")

# Scale within clade
combined_test[, pred_scaled := scale(pred), by = clade]
combined_test[, pheno_scaled := scale(pheno), by = clade]

# Calc within-clade cors
clade_results <- combined_test[, .(
  n = .N,
  spearman = cor(pred, pheno, method = "spearman", use = "complete.obs")
), by = clade]

message("Mean within-clade Spearman: ", round(mean(clade_results$spearman, na.rm = TRUE), 3))

# THE PLOT
plot_data <- combined_test[!is.na(pred_scaled) & !is.na(pheno_scaled)]

p <- ggplot(plot_data, aes(x = pred_scaled, y = pheno_scaled)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Within-Clade Prediction: Best Model",
    subtitle = paste0("n_top_genes=", best$n_top_genes, ", n_dims=", best$n_dims_per_gene,
                      ", threshold=", best$cor_threshold, 
                      "\nMean within-clade r = ", round(mean(clade_results$spearman, na.rm = TRUE), 3)),
    x = "Predicted (z-scored within clade)",
    y = "Observed (z-scored within clade)"
  ) +
  theme_minimal() +
  coord_fixed()

print(p)
ggsave("results/phylo_cv_best_model_scaled.png", p, width = 7, height = 7)

png("results/phylo_cv_best_model_scaled_heatscatter.png", width = 7, height = 7, units = "in", res = 300)
LSD::heatscatter(plot_data$pred_scaled, plot_data$pheno_scaled, xlab = "Predicted (z-scored within clade)", ylab = "Observed (z-scored within clade)", main = paste0("Within-Clade Prediction: Best Model\nMean within-clade r = ", round(mean(clade_results$spearman, na.rm = TRUE), 3)))
abline(a = 0, b = 1, col = "red", lty = "dashed")
dev.off()
# Also show clade results
print(clade_results[order(-spearman)])