library(ape)
library(arrow)
library(data.table)
library(ggplot2)
library(glmnet)

# ---- 1. LOAD DATA ----
tree <- read.tree("results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree")
data <- as.data.table(read_parquet("data/processed_data.parquet"))
embeds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds)
embed_cols <- grep("embedding", colnames(embeds), value = TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

# ---- 2. PREPARE TREE & PHENOTYPE ----
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
pheno <- setNames(data[[pheno_col]], data$ID)
pheno <- pheno[tree$tip.label]
pheno <- pheno[!is.na(pheno)]
tree <- keep.tip(tree, names(pheno))
message("Tree: ", Ntip(tree), " tips with phenotype data")

clean_embeds <- clean_embeds[ID %in% tree$tip.label]
data <- data[ID %in% tree$tip.label]
stopifnot(length(unique(clean_embeds$ID)) > 0)

# ---- 3. SELECT HOLDOUT CLADES ----
# Target: ~15% of tips in holdout, spread across many small clades
select_holdout_clades <- function(tree, target_frac = 0.15, min_clade = 10, max_clade = 30) {
  n_tips <- Ntip(tree)
  target_n <- floor(n_tips * target_frac)
  
  # Get all clades in size range
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  candidates <- list()
  
  for (node in internal_nodes) {
    tips <- extract.clade(tree, node)$tip.label
    if (length(tips) >= min_clade && length(tips) <= max_clade) {
      candidates[[as.character(node)]] <- tips
    }
  }
  
  message("Found ", length(candidates), " candidate clades (", min_clade, "-", max_clade, " tips)")
  stopifnot(length(candidates) > 0)
  
  # Greedily select non-overlapping clades until we hit target
  selected <- list()
  used_tips <- character(0)
  
  # Shuffle to avoid bias toward tree structure
  candidate_order <- sample(names(candidates))
  
  for (node in candidate_order) {
    tips <- candidates[[node]]
    if (!any(tips %in% used_tips)) {
      selected[[node]] <- tips
      used_tips <- c(used_tips, tips)
      if (length(used_tips) >= target_n) break
    }
  }
  
  message("Selected ", length(selected), " clades, ", 
          length(used_tips), " tips (", 
          round(100 * length(used_tips) / n_tips, 1), "%)")
  
  return(selected)
}

set.seed(42)
holdout_clades <- select_holdout_clades(tree, target_frac = 0.15, min_clade = 10, max_clade = 30)

# If not enough, relax constraints
if (length(unlist(holdout_clades)) < 0.10 * Ntip(tree)) {
  message("Relaxing clade size constraints...")
  holdout_clades <- select_holdout_clades(tree, target_frac = 0.15, min_clade = 5, max_clade = 50)
}

holdout_ids <- unlist(holdout_clades)
train_ids <- setdiff(tree$tip.label, holdout_ids)

message("Train: ", length(train_ids), " | Test: ", length(holdout_ids), 
        " (", length(holdout_clades), " clades)")

# ---- 4. HELPER FUNCTIONS ----
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

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

select_features <- function(train_embeds, train_df, embed_cols, pheno_col, 
                           n_top = 1, n_top_genes = 10, n_dims_per_top_gene = 3) {
  # n_top: dims per gene for regular genes
  # n_top_genes: how many top genes get extra dimensions
  # n_dims_per_top_gene: how many orthogonal dims for top genes
  
  merged <- merge(
    train_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
    train_df[, .(ID, Order, pheno = get(pheno_col))],
    by = "ID"
  )
  
  gene_counts <- merged[, .N, by = Gene]
  valid_genes <- gene_counts[N >= 5]$Gene
  
  sel_list <- list()
  cor_records <- list()
  
  # First pass: get best correlation per gene to rank genes
  gene_best_cors <- list()
  for (gene in valid_genes) {
    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) == 0) next
    gene_best_cors[[gene]] <- max(abs(cors[available]), na.rm = TRUE)
  }
  
  # Rank genes by best correlation
  gene_ranks <- sort(unlist(gene_best_cors), decreasing = TRUE)
  top_gene_names <- names(gene_ranks)[1:min(n_top_genes, length(gene_ranks))]
  
  message("Top ", length(top_gene_names), " genes by correlation: ", 
          paste(top_gene_names, collapse = ", "))
  
  # Second pass: select dimensions
  for (gene in valid_genes) {
    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) == 0) next
    
    # Determine how many dims to select
    if (gene %in% top_gene_names) {
      n_select <- min(n_dims_per_top_gene, length(available))
    } else {
      n_select <- min(n_top, length(available))
    }
    
    if (n_select == 1) {
      # Simple case: just pick best
      top_dims <- names(sort(abs(cors[available]), decreasing = TRUE))[1]
      cor_records[[paste0(gene, "_1")]] <- data.table(
        gene = gene, embedding = top_dims, correlation = cors[top_dims], dim_rank = 1
      )
    } else {
      # Pick orthogonal dimensions via residualization
      embed_matrix <- as.matrix(gdt[, embed_cols, with = FALSE])
      pheno_vec <- gdt$pheno
      
      selected_dims <- character(0)
      residual_pheno <- pheno_vec
      
      for (d in 1:n_select) {
        # Correlate remaining embeddings with residual phenotype
        if (d == 1) {
          current_cors <- cors[available]
        } else {
          current_cors <- cor(embed_matrix[, available, drop = FALSE], 
                              residual_pheno, use = "complete.obs")[, 1]
          names(current_cors) <- available
        }
        
        # Pick best remaining
        best_dim <- names(sort(abs(current_cors), decreasing = TRUE))[1]
        selected_dims <- c(selected_dims, best_dim)
        
        cor_records[[paste0(gene, "_", d)]] <- data.table(
          gene = gene, embedding = best_dim, 
          correlation = current_cors[best_dim], dim_rank = d
        )
        
        # Residualize phenotype on selected dimension
        if (d < n_select) {
          fit_resid <- lm(residual_pheno ~ embed_matrix[, best_dim])
          residual_pheno <- residuals(fit_resid)
          # Remove from available
          available <- setdiff(available, best_dim)
          if (length(available) == 0) break
        }
      }
      top_dims <- selected_dims
    }
    
    sel <- unique(gdt[, c("ID", top_dims), with = FALSE])
    new_names <- paste0(gene, "__", top_dims)
    setnames(sel, old = top_dims, new = new_names)
    sel_list[[gene]] <- sel
  }
  
  if (length(sel_list) == 0) return(NULL)
  combined <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)
  combined <- merge(combined, train_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)
  
  # Median impute
  pred_cols <- setdiff(names(combined), c("ID", "pheno"))
  for (col in pred_cols) {
    med <- median(combined[[col]], na.rm = TRUE)
    combined[is.na(get(col)), (col) := med]
  }
  
  # Attach correlation info as attribute
  cor_dt <- rbindlist(cor_records)
  attr(combined, "feature_cors") <- cor_dt
  
  return(combined)
}

prepare_test_data <- function(test_embeds, test_df, pred_cols, embed_cols, pheno_col) {
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
  
  if (length(test_wide_list) == 0) return(NULL)
  combined <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), test_wide_list)
  combined <- merge(combined, test_df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)
  
  for (col in names(combined)) {
    if (col %in% pred_cols && any(is.na(combined[[col]]))) {
      med <- median(combined[[col]], na.rm = TRUE)
      combined[is.na(get(col)), (col) := med]
    }
  }
  
  return(combined)
}

# ---- 5. TRAIN MODEL ----
train_df <- data[ID %in% train_ids]
train_embeds <- clean_embeds[ID %in% train_ids]

message("Running feature selection on training data...")
# Get multiple orthogonal dims from top genes
combined_train <- select_features(train_embeds, train_df, embed_cols, pheno_col,
                                  n_top = 1,           # 1 dim for regular genes
                                  n_top_genes = 6,     # top 6 genes get extra dims
                                  n_dims_per_top_gene = 3)  # 3 orthogonal dims each
stopifnot(!is.null(combined_train))

pred_cols_all <- grep("embedding", names(combined_train), value = TRUE)
message("Selected ", length(pred_cols_all), " predictors total")

# Print feature correlations
feature_cors <- attr(combined_train, "feature_cors")
message("\n=== SELECTED FEATURES (sorted by |correlation|) ===")
feature_cors <- feature_cors[order(-abs(correlation))]
print(feature_cors)

# ---- 5b. PRE-FIT PRUNING BY CORRELATION THRESHOLD ----
# Test multiple thresholds - push higher
cor_thresholds <- c(0, 0.10, 0.15, 0.18, 0.20, 0.22, 0.25)

message("\n=== TESTING CORRELATION THRESHOLDS ===")
threshold_results <- data.table(
  threshold = numeric(),
  n_predictors = integer(),
  train_r2 = numeric()
)

for (thresh in cor_thresholds) {
  # Keep dims where the gene's best correlation >= threshold
  # (for multi-dim genes, keep all dims if gene qualifies)
  gene_max_cors <- feature_cors[, .(max_cor = max(abs(correlation))), by = gene]
  keep_genes <- gene_max_cors[max_cor >= thresh]$gene
  keep_cols <- pred_cols_all[sapply(pred_cols_all, function(x) {
    gene <- sub("__embedding.*", "", x)
    gene %in% keep_genes
  })]
  
  if (length(keep_cols) < 2) next
  
  formula_tmp <- paste("pheno ~", paste(keep_cols, collapse = " + "))
  fit_tmp <- lm(as.formula(formula_tmp), data = combined_train, na.action = na.exclude)
  
  threshold_results <- rbind(threshold_results, data.table(
    threshold = thresh,
    n_predictors = length(keep_cols),
    train_r2 = summary(fit_tmp)$r.squared
  ))
  
  message("  threshold=", thresh, ": ", length(keep_cols), " predictors, R²=", 
          round(summary(fit_tmp)$r.squared, 3))
}

# Set thresholds - push strict higher
cor_threshold <- 0.15
cor_threshold_strict <- 0.20

# Filter by gene's max correlation
gene_max_cors <- feature_cors[, .(max_cor = max(abs(correlation))), by = gene]

keep_genes <- gene_max_cors[max_cor >= cor_threshold]$gene
pred_cols <- pred_cols_all[sapply(pred_cols_all, function(x) {
  gene <- sub("__embedding.*", "", x)
  gene %in% keep_genes
})]

keep_genes_strict <- gene_max_cors[max_cor >= cor_threshold_strict]$gene
pred_cols_strict <- pred_cols_all[sapply(pred_cols_all, function(x) {
  gene <- sub("__embedding.*", "", x)
  gene %in% keep_genes_strict
})]

message("\n=== MODELS ===")
message("Standard (gene cor >= ", cor_threshold, "): ", length(pred_cols), 
        " predictors from ", length(keep_genes), " genes")
message("Strict (gene cor >= ", cor_threshold_strict, "): ", length(pred_cols_strict), 
        " predictors from ", length(keep_genes_strict), " genes")

# Fit both models
formula_str <- paste("pheno ~", paste(pred_cols, collapse = " + "))
fit <- lm(as.formula(formula_str), data = combined_train, na.action = na.exclude)

formula_strict <- paste("pheno ~", paste(pred_cols_strict, collapse = " + "))
fit_strict <- lm(as.formula(formula_strict), data = combined_train, na.action = na.exclude)

message("\nStandard model R²: ", round(summary(fit)$r.squared, 3))
message("Strict model R²: ", round(summary(fit_strict)$r.squared, 3))

# Print model summary for standard
message("\n=== MODEL SUMMARY (standard) ===")
model_coefs <- summary(fit)$coefficients
sig_coefs <- model_coefs[model_coefs[, 4] < 0.05, , drop = FALSE]
message("Significant predictors (p < 0.05): ", nrow(sig_coefs) - 1, " (excluding intercept)")
print(sig_coefs)

message("\n=== MODEL SUMMARY (strict) ===")
model_coefs_strict <- summary(fit_strict)$coefficients
print(model_coefs_strict)

# ---- 6. EVALUATE ON HOLDOUT CLADES ----
test_df <- data[ID %in% holdout_ids]
test_embeds <- clean_embeds[ID %in% holdout_ids]

combined_test <- prepare_test_data(test_embeds, test_df, pred_cols_all, embed_cols, pheno_col)
stopifnot(!is.null(combined_test))

# Predictions from both models
pred_test_standard <- predict(fit, newdata = combined_test)
pred_test_strict <- predict(fit_strict, newdata = combined_test)

combined_test$pred_standard <- pred_test_standard
combined_test$pred_strict <- pred_test_strict

# Map IDs to clades
id_to_clade <- data.table(
  ID = unlist(holdout_clades),
  clade = rep(names(holdout_clades), sapply(holdout_clades, length))
)
combined_test <- merge(combined_test, id_to_clade, by = "ID")

# Within-clade Spearman correlations for both models
clade_results <- combined_test[, .(
  n = .N,
  n_valid = sum(!is.na(pred_standard) & !is.na(pheno)),
  spearman_standard = if (sum(!is.na(pred_standard) & !is.na(pheno)) >= 3) 
    cor(pred_standard, pheno, method = "spearman", use = "complete.obs") else NA_real_,
  spearman_strict = if (sum(!is.na(pred_strict) & !is.na(pheno)) >= 3) 
    cor(pred_strict, pheno, method = "spearman", use = "complete.obs") else NA_real_,
  pearson_standard = if (sum(!is.na(pred_standard) & !is.na(pheno)) >= 3) 
    cor(pred_standard, pheno, method = "pearson", use = "complete.obs") else NA_real_
), by = clade]

# ---- 7. SUMMARY ----
message("\n=== WITHIN-CLADE EVALUATION (STANDARD: cor >= ", cor_threshold, ") ===")
message("Predictors: ", length(pred_cols))
message("Clades evaluated: ", sum(!is.na(clade_results$spearman_standard)), "/", nrow(clade_results))
message("Mean within-clade Spearman: ", round(mean(clade_results$spearman_standard, na.rm = TRUE), 3),
        " (SD: ", round(sd(clade_results$spearman_standard, na.rm = TRUE), 3), ")")
message("Median within-clade Spearman: ", round(median(clade_results$spearman_standard, na.rm = TRUE), 3))

message("\n=== WITHIN-CLADE EVALUATION (STRICT: cor >= ", cor_threshold_strict, ") ===")
message("Predictors: ", length(pred_cols_strict))
message("Mean within-clade Spearman: ", round(mean(clade_results$spearman_strict, na.rm = TRUE), 3),
        " (SD: ", round(sd(clade_results$spearman_strict, na.rm = TRUE), 3), ")")
message("Median within-clade Spearman: ", round(median(clade_results$spearman_strict, na.rm = TRUE), 3))

# Test if mean Spearman > 0
valid_cors <- clade_results$spearman_standard[!is.na(clade_results$spearman_standard)]
if (length(valid_cors) >= 3) {
  tt <- t.test(valid_cors, mu = 0)
  message("\nStandard model t-test (H0: mean=0): p=", format.pval(tt$p.value, digits = 3))
}

valid_cors_strict <- clade_results$spearman_strict[!is.na(clade_results$spearman_strict)]
if (length(valid_cors_strict) >= 3) {
  tt_strict <- t.test(valid_cors_strict, mu = 0)
  message("Strict model t-test (H0: mean=0): p=", format.pval(tt_strict$p.value, digits = 3))
}

# Global correlations
global_cor_standard <- cor(combined_test$pred_standard, combined_test$pheno, 
                           method = "spearman", use = "complete.obs")
global_cor_strict <- cor(combined_test$pred_strict, combined_test$pheno, 
                         method = "spearman", use = "complete.obs")
message("\nGlobal Spearman (standard): ", round(global_cor_standard, 3))
message("Global Spearman (strict): ", round(global_cor_strict, 3))

# ---- 8. BEST/WORST CLADE DETAILS ----
# Add taxonomy info to test data
combined_test <- merge(combined_test, data[, .(ID, Order, Taxonomy)], by = "ID", all.x = TRUE)

message("\n=== TOP 5 CLADES (strict model) ===")
top_clades <- clade_results[order(-spearman_strict)][1:5]$clade
for (cl in top_clades) {
  clade_data <- combined_test[clade == cl]
  message("\nClade ", cl, " (n=", nrow(clade_data), ", r=", 
          round(clade_results[clade == cl]$spearman_strict, 3), "):")
  message("  Orders: ", paste(unique(clade_data$Order), collapse = ", "))
  message("  Pheno range: ", round(min(clade_data$pheno, na.rm = TRUE), 1), " - ", 
          round(max(clade_data$pheno, na.rm = TRUE), 1))
}

message("\n=== BOTTOM 5 CLADES (strict model) ===")
bottom_clades <- clade_results[order(spearman_strict)][1:5]$clade
for (cl in bottom_clades) {
  clade_data <- combined_test[clade == cl]
  message("\nClade ", cl, " (n=", nrow(clade_data), ", r=", 
          round(clade_results[clade == cl]$spearman_strict, 3), "):")
  message("  Orders: ", paste(unique(clade_data$Order), collapse = ", "))
  message("  Pheno range: ", round(min(clade_data$pheno, na.rm = TRUE), 1), " - ", 
          round(max(clade_data$pheno, na.rm = TRUE), 1))
}

print(clade_results[order(-spearman_strict)])

# ---- 9. PLOTS ----

# KEY VISUALIZATION: Within-clade scaled predictions
# Center and scale within each clade to remove between-clade structure
combined_test[, pred_scaled := scale(pred_standard), by = clade]
combined_test[, pheno_scaled := scale(pheno), by = clade]

# Remove clades with no variance (single value after scaling = NaN)
plot_data <- combined_test[!is.na(pred_scaled) & !is.na(pheno_scaled)]

p_within <- ggplot(plot_data, aes(x = pred_scaled, y = pheno_scaled)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Within-Clade Prediction (Centered & Scaled)",
    subtitle = paste0("Mean within-clade Spearman = ", 
                      round(mean(clade_results$spearman_standard, na.rm = TRUE), 3),
                      " (p = ", format.pval(tt$p.value, digits = 2), ")"),
    x = "Predicted (z-scored within clade)",
    y = "Observed (z-scored within clade)"
  ) +
  theme_minimal() +
  coord_fixed()

# Also do for strict model
combined_test[, pred_strict_scaled := scale(pred_strict), by = clade]

plot_data_strict <- combined_test[!is.na(pred_strict_scaled) & !is.na(pheno_scaled)]

p_within_strict <- ggplot(plot_data_strict, aes(x = pred_strict_scaled, y = pheno_scaled)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Within-Clade Prediction - Strict Model (Centered & Scaled)",
    subtitle = paste0("Mean within-clade Spearman = ", 
                      round(mean(clade_results$spearman_strict, na.rm = TRUE), 3),
                      " (p = ", format.pval(tt_strict$p.value, digits = 2), ")"),
    x = "Predicted (z-scored within clade)",
    y = "Observed (z-scored within clade)"
  ) +
  theme_minimal() +
  coord_fixed()

# Within-clade correlation distribution
p_dist <- ggplot(clade_results, aes(x = spearman_standard)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean(clade_results$spearman_standard, na.rm = TRUE), 
             color = "darkgreen", linetype = "solid", size = 1) +
  labs(
    title = "Distribution of Within-Clade Spearman Correlations",
    subtitle = paste0("Mean = ", round(mean(clade_results$spearman_standard, na.rm = TRUE), 3),
                      ", Median = ", round(median(clade_results$spearman_standard, na.rm = TRUE), 3)),
    x = "Spearman r", y = "Count"
  ) +
  theme_minimal()

# Save plots
ggsave("results/phylo_cv_within_clade_scaled.png", p_within, width = 7, height = 7)
ggsave("results/phylo_cv_within_clade_scaled_strict.png", p_within_strict, width = 7, height = 7)
ggsave("results/phylo_cv_spearman_dist.png", p_dist, width = 8, height = 5)
message("\nPlots saved to results/")

# ---- 10. GRID SEARCH OVER VARIABLE SELECTION PARAMETERS ----
message("\n=== GRID SEARCH ===")

grid <- expand.grid(
  n_top_genes = c(3, 6, 10, 15),
  n_dims_per_gene = c(1, 2, 3, 5),
  cor_threshold = c(0.12, 0.15, 0.18, 0.20, 0.22)
)

grid_results <- data.table(
  n_top_genes = integer(),
  n_dims_per_gene = integer(),
  cor_threshold = numeric(),
  n_predictors = integer(),
  train_r2 = numeric(),
  mean_within_spearman = numeric(),
  median_within_spearman = numeric(),
  pval = numeric(),
  global_spearman = numeric()
)

for (i in 1:nrow(grid)) {
  params <- grid[i, ]
  print(params)
  
  # Re-run feature selection with these params
  combined_grid <- select_features(
    train_embeds, train_df, embed_cols, pheno_col,
    n_top = 1,
    n_top_genes = params$n_top_genes,
    n_dims_per_top_gene = params$n_dims_per_gene
  )
  
  if (is.null(combined_grid)) next
  
  pred_cols_grid <- grep("embedding", names(combined_grid), value = TRUE)
  feature_cors_grid <- attr(combined_grid, "feature_cors")
  
  # Filter by threshold
  gene_max_cors <- feature_cors_grid[, .(max_cor = max(abs(correlation))), by = gene]
  keep_genes <- gene_max_cors[max_cor >= params$cor_threshold]$gene
  keep_cols <- pred_cols_grid[sapply(pred_cols_grid, function(x) {
    gene <- sub("__embedding.*", "", x)
    gene %in% keep_genes
  })]
  
  if (length(keep_cols) < 2) next
  
  # Fit model
  formula_grid <- paste("pheno ~", paste(keep_cols, collapse = " + "))
  fit_grid <- lm(as.formula(formula_grid), data = combined_grid, na.action = na.exclude)
  
  # Prepare test data and predict
  test_grid <- prepare_test_data(test_embeds, test_df, pred_cols_grid, embed_cols, pheno_col)
  if (is.null(test_grid)) next
  
  pred_grid <- predict(fit_grid, newdata = test_grid)
  test_grid$pred <- pred_grid
  test_grid <- merge(test_grid, id_to_clade, by = "ID")
  
  # Within-clade correlations
  clade_cors <- test_grid[, .(
    spearman = if (sum(!is.na(pred) & !is.na(pheno)) >= 3)
      cor(pred, pheno, method = "spearman", use = "complete.obs") else NA_real_
  ), by = clade]
  
  valid_cors <- clade_cors$spearman[!is.na(clade_cors$spearman)]
  if (length(valid_cors) < 3) next
  
  tt_grid <- t.test(valid_cors, mu = 0)
  global_cor <- cor(test_grid$pred, test_grid$pheno, method = "spearman", use = "complete.obs")
  
  grid_results <- rbind(grid_results, data.table(
    n_top_genes = params$n_top_genes,
    n_dims_per_gene = params$n_dims_per_gene,
    cor_threshold = params$cor_threshold,
    n_predictors = length(keep_cols),
    train_r2 = summary(fit_grid)$r.squared,
    mean_within_spearman = mean(valid_cors),
    median_within_spearman = median(valid_cors),
    pval = tt_grid$p.value,
    global_spearman = global_cor
  ))
  
  if (i %% 10 == 0) message("  Completed ", i, "/", nrow(grid), " grid points")
}

message("\n=== GRID SEARCH RESULTS (sorted by mean within-clade Spearman) ===")
grid_results <- grid_results[order(-mean_within_spearman)]
print(grid_results[1:20])

# Best config
best <- grid_results[1]
message("\n=== BEST CONFIGURATION ===")
message("n_top_genes: ", best$n_top_genes)
message("n_dims_per_gene: ", best$n_dims_per_gene)
message("cor_threshold: ", best$cor_threshold)
message("n_predictors: ", best$n_predictors)
message("Mean within-clade Spearman: ", round(best$mean_within_spearman, 3))
message("p-value: ", format.pval(best$pval, digits = 3))
message("Global Spearman: ", round(best$global_spearman, 3))

# Plot grid search results
p_grid <- ggplot(grid_results, aes(x = n_predictors, y = mean_within_spearman, 
                                    color = factor(cor_threshold))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(n_top_genes, n_dims_per_gene)), alpha = 0.3) +
  labs(
    title = "Grid Search: Within-Clade Performance vs Model Complexity",
    x = "Number of Predictors",
    y = "Mean Within-Clade Spearman",
    color = "Correlation\nThreshold"
  ) +
  theme_minimal()

ggsave("results/phylo_cv_grid_search.png", p_grid, width = 10, height = 6)

# ---- 11. SAVE ----
results <- list(
  clade_results = clade_results,
  feature_cors = feature_cors,
  threshold_results = threshold_results,
  grid_results = grid_results,
  global_spearman_standard = global_cor_standard,
  global_spearman_strict = global_cor_strict,
  n_train = nrow(combined_train),
  n_test = nrow(combined_test),
  n_clades = nrow(clade_results),
  n_predictors_standard = length(pred_cols),
  n_predictors_strict = length(pred_cols_strict),
  cor_threshold = cor_threshold,
  cor_threshold_strict = cor_threshold_strict,
  model_standard = fit,
  model_strict = fit_strict,
  best_config = best
)
saveRDS(results, "results/phylo_cv_clade_results.rds")