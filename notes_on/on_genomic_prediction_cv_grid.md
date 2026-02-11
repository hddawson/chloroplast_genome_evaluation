    ## Run: 2026-02-10T20:26:56

    **Script**: `src/analyses/genomic_prediction_cv_grid.r`
    **Interpreter**: `Rscript`
    **Hash**: `88037e3a1056`

    ### Summary
    This script loads cross-validation results from a phylogenetic machine learning analysis and recreates the best-performing model configuration to generate predictions and visualizations. It uses protein embeddings as features to predict a climate-related phenotype (bio_8_p50), applying the optimal hyperparameters (number of top genes, dimensions per gene, and correlation threshold) identified from grid search. The script performs within-clade scaling of both predictions and observed values, then creates a scatter plot showing the relationship between predicted and observed phenotypes scaled within phylogenetic clades. It calculates and displays within-clade Spearman correlations to evaluate model performance while controlling for phylogenetic relationships.

    ### Script
    ```
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

# Also show clade results
print(clade_results[order(-spearman)])
    ```

    ### Output
    ```
        n_top_genes n_dims_per_gene cor_threshold n_predictors  train_r2
          <num>           <num>         <num>        <int>     <num>
 1:          15               5          0.18           55 0.3924839
 2:          10               5          0.18           40 0.3746781
 3:          10               5          0.15           50 0.3845307
 4:          10               5          0.20           30 0.3666204
 5:          15               5          0.20           40 0.3823032
 6:           6               5          0.15           33 0.3650585
 7:           6               5          0.18           30 0.3604315
 8:          15               5          0.15           70 0.4013329
 9:          15               5          0.12           85 0.4103527
10:          10               5          0.12           65 0.4000511
    mean_within_spearman median_within_spearman         pval global_spearman
                   <num>                  <num>        <num>           <num>
 1:            0.2292242              0.3090909 1.448731e-06       0.5553862
 2:            0.2270809              0.2484848 1.623835e-07       0.5508764
 3:            0.2221841              0.2456140 2.489701e-08       0.5693558
 4:            0.2199926              0.2363636 9.412614e-08       0.5506997
 5:            0.2189764              0.2647059 5.380559e-07       0.5419249
 6:            0.2183719              0.2311983 3.355142e-07       0.5652679
 7:            0.2183437              0.2235294 9.560154e-07       0.5604816
 8:            0.2182237              0.2500000 1.403115e-06       0.5728868
 9:            0.2175120              0.2466460 9.184730e-06       0.5837983
10:            0.2107036              0.2382353 3.950978e-06       0.5814215
   n_top_genes n_dims_per_gene cor_threshold n_predictors  train_r2
         <num>           <num>         <num>        <int>     <num>
1:          15               5          0.18           55 0.3924839
   mean_within_spearman median_within_spearman         pval global_spearman
                  <num>                  <num>        <num>           <num>
1:            0.2292242              0.3090909 1.448731e-06       0.5553862
     clade     n      spearman
    <char> <int>         <num>
 1:   8230    11  0.8009085750
 2:   7195    19  0.6263157895
 3:  10183    28  0.6212370005
 4:   9660    18  0.5990708091
 5:   7167    18  0.5761487603
 6:   8604    13  0.5219780220
 7:  10271    11  0.5068545991
 8:   7026    12  0.5034965035
 9:   8427    18  0.4994840052
10:  11811    11  0.4862589963
11:   7569    10  0.4817162723
12:   6114    15  0.4704837582
13:   9548    11  0.4692495090
14:   8725    25  0.4362377461
15:   9812    16  0.4317678503
16:   8130    20  0.3699248120
17:   8574    14  0.3626373626
18:   7835    29  0.3291045720
19:   7596    29  0.3247505458
20:   7384    16  0.3235294118
21:   9513    13  0.3034508738
22:   8478    13  0.2971608194
23:  11944    17  0.2858947222
24:   9429    20  0.2527266068
25:   7725    14  0.2497275683
26:   7995    19  0.1826164005
27:   9327    25  0.1800346254
28:  11215    22  0.1724627112
29:  11706    17  0.1523346121
30:  12138    11  0.1454545455
31:   6541    16  0.1411764706
32:  12002    11  0.1369877295
33:  12112    25  0.1350259690
34:   7484    28  0.1212868001
35:   6924    13  0.0990372326
36:  12100    12  0.0595447498
37:   6395    15  0.0535714286
38:  10385    16  0.0264705882
39:  10718    15  0.0232350406
40:   6942    25  0.0165416429
41:   8198    11  0.0045558205
42:  10034    15  0.0017905168
43:   8100    24 -0.0004349718
44:   7132    12 -0.0069930070
45:   8842    18 -0.0092879257
46:  10927    21 -0.0207927269
47:  11882    10 -0.0547114989
48:   8075    26 -0.0925690811
49:  10353    10 -0.1043023303
50:   6268    12 -0.1048951049
51:   8860    23 -0.2643933859
52:  10127    10 -0.2969696970
53:  11690    11 -0.4181818182
54:  11970    10 -0.6242424242
55:   8457    10 -0.6319494127
     clade     n      spearman

--- stderr ---

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

There were 27 warnings (use warnings() to see them)
There were 27 warnings (use warnings() to see them)
Best config: 15 predictors
Mean within-clade Spearman: 0.185
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'

    ```

    ---

