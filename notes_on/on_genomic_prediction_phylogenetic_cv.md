    ## Run: 2026-02-10T18:20:00

    **Script**: `src/analyses/genomic_prediction_phylogenetic_cv.r`
    **Interpreter**: `Rscript`
    **Hash**: `8f73f7f7068c`

    ### Summary
    This script performs phylogenetic cross-validation to evaluate protein embedding-based prediction of climate adaptation traits in species. It uses a phylogenetic tree and protein embeddings to predict a climate variable (bio_8_p50, likely related to temperature of wettest quarter) by selecting holdout clades (10-30 species each, ~15% of data) as test sets while training on the remaining phylogenetically distant species. The analysis compares two linear regression models with different feature selection stringency based on correlation thresholds (0.12 vs 0.16), evaluating performance using within-clade Spearman correlations to assess how well the models generalize across phylogenetic distances. The script outputs detailed results including per-clade correlations, global performance metrics, and visualizations comparing predicted vs observed values.

    ### Script
    ```
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

select_features <- function(train_embeds, train_df, embed_cols, pheno_col, n_top = 1) {
  merged <- merge(
    train_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
    train_df[, .(ID, Order, pheno = get(pheno_col))],
    by = "ID"
  )

  gene_counts <- merged[, .N, by = Gene]
  valid_genes <- gene_counts[N >= 5]$Gene

  sel_list <- list()
  cor_records <- list()  # Track correlations

  for (gene in valid_genes) {
    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) == 0) next
    k <- min(n_top, length(available))
    top_dims <- names(sort(abs(cors[available]), decreasing = TRUE))[1:k]

    # Record correlation for selected dimension
    cor_records[[gene]] <- data.table(
      gene = gene,
      embedding = top_dims,
      correlation = cors[top_dims]
    )

    sel <- unique(gdt[, c("ID", top_dims), with = FALSE])
    setnames(sel, old = top_dims, new = paste0(gene, "__", top_dims))
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
combined_train <- select_features(train_embeds, train_df, embed_cols, pheno_col)
stopifnot(!is.null(combined_train))

pred_cols_all <- grep("embedding", names(combined_train), value = TRUE)
message("Selected ", length(pred_cols_all), " predictors")

# Print feature correlations
feature_cors <- attr(combined_train, "feature_cors")
message("\n=== SELECTED FEATURES (sorted by |correlation|) ===")
feature_cors <- feature_cors[order(-abs(correlation))]
print(feature_cors)

# ---- 5b. PRE-FIT PRUNING BY CORRELATION THRESHOLD ----
# Test multiple thresholds
cor_thresholds <- c(0, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20)

message("\n=== TESTING CORRELATION THRESHOLDS ===")
threshold_results <- data.table(
  threshold = numeric(),
  n_predictors = integer(),
  train_r2 = numeric()
)

for (thresh in cor_thresholds) {
  keep_genes <- feature_cors[abs(correlation) >= thresh]$gene
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

# Use 0.12 as default (keeps ~25 predictors based on your output)
# But also test a stricter one
cor_threshold <- 0.12
cor_threshold_strict <- 0.16

keep_genes <- feature_cors[abs(correlation) >= cor_threshold]$gene
pred_cols <- pred_cols_all[sapply(pred_cols_all, function(x) {
  gene <- sub("__embedding.*", "", x)
  gene %in% keep_genes
})]

keep_genes_strict <- feature_cors[abs(correlation) >= cor_threshold_strict]$gene
pred_cols_strict <- pred_cols_all[sapply(pred_cols_all, function(x) {
  gene <- sub("__embedding.*", "", x)
  gene %in% keep_genes_strict
})]

message("\n=== MODELS ===")
message("Standard (cor >= ", cor_threshold, "): ", length(pred_cols), " predictors")
message("Strict (cor >= ", cor_threshold_strict, "): ", length(pred_cols_strict), " predictors")

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
combined_test <- merge(combined_test, data[, .(ID, Order, Family, Taxonomy)], by = "ID", all.x = TRUE)

message("\n=== TOP 5 CLADES (strict model) ===")
top_clades <- clade_results[order(-spearman_strict)][1:5]$clade
for (cl in top_clades) {
  clade_data <- combined_test[clade == cl]
  message("\nClade ", cl, " (n=", nrow(clade_data), ", r=", 
          round(clade_results[clade == cl]$spearman_strict, 3), "):")
  message("  Orders: ", paste(unique(clade_data$Order), collapse = ", "))
  message("  Families: ", paste(unique(clade_data$Family), collapse = ", "))
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
  message("  Families: ", paste(unique(clade_data$Family), collapse = ", "))
  message("  Pheno range: ", round(min(clade_data$pheno, na.rm = TRUE), 1), " - ", 
          round(max(clade_data$pheno, na.rm = TRUE), 1))
}

print(clade_results[order(-spearman_strict)])

# ---- 9. PLOTS ----
# Global prediction plot
p1 <- ggplot(combined_test, aes(x = pred_strict, y = pheno)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = paste0("Predicted vs Observed (Strict: cor >= ", cor_threshold_strict, ")"),
    subtitle = paste0("Global Spearman = ", round(global_cor_strict, 3)),
    x = "Predicted", y = "Observed"
  ) +
  theme_minimal()

# Faceted by top/bottom clades
highlight_clades <- c(top_clades[1:3], bottom_clades[1:3])
p2 <- ggplot(combined_test[clade %in% highlight_clades], 
             aes(x = pred_strict, y = pheno)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~clade, scales = "free") +
  labs(
    title = "Top 3 vs Bottom 3 Clades (strict model)",
    x = "Predicted", y = "Observed"
  ) +
  theme_minimal()

# Within-clade correlation distribution
p3 <- ggplot(clade_results, aes(x = spearman_strict)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean(clade_results$spearman_strict, na.rm = TRUE), 
             color = "darkgreen", linetype = "solid") +
  labs(
    title = "Distribution of Within-Clade Spearman Correlations (strict)",
    subtitle = paste0("Mean = ", round(mean(clade_results$spearman_strict, na.rm = TRUE), 3)),
    x = "Spearman r", y = "Count"
  ) +
  theme_minimal()

# Compare standard vs strict
p4 <- ggplot(clade_results, aes(x = spearman_standard, y = spearman_strict)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Standard vs Strict Model (Within-Clade Spearman)",
    x = paste0("Standard (cor >= ", cor_threshold, ")"), 
    y = paste0("Strict (cor >= ", cor_threshold_strict, ")")
  ) +
  theme_minimal()

# Save plots
ggsave("results/phylo_cv_global_pred.png", p1, width = 8, height = 6)
ggsave("results/phylo_cv_clade_facets.png", p2, width = 10, height = 6)
ggsave("results/phylo_cv_spearman_dist.png", p3, width = 8, height = 5)
ggsave("results/phylo_cv_standard_vs_strict.png", p4, width = 6, height = 6)
message("\nPlots saved to results/")

# ---- 10. SAVE ----
results <- list(
  clade_results = clade_results,
  feature_cors = feature_cors,
  threshold_results = threshold_results,
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
  model_strict = fit_strict
)
saveRDS(results, "results/phylo_cv_clade_results.rds")
    ```

    ### Output
    ```
          gene     embedding correlation
    <char>        <char>       <num>
 1:   rbcL embedding_916 -0.25250013
 2:   psbJ  embedding_72  0.22293014
 3:   psbB embedding_295  0.19134844
 4:   psaB embedding_199  0.17085064
 5:   psbK  embedding_57  0.16038676
 6:   psbM embedding_815  0.16022707
 7:   psbT embedding_803 -0.15698894
 8:   petA embedding_906  0.15426775
 9:   psbD embedding_897  0.15215471
10:   atpA embedding_171  0.13902136
11:   rpoA embedding_805 -0.13603285
12:   psbC  embedding_36  0.13510181
13:   psbF embedding_546  0.13503082
14:   ycf3 embedding_381 -0.13452911
15:  rpl36 embedding_716  0.13243696
16:   ndhC embedding_228 -0.12960872
17:  rpl14 embedding_388  0.12888506
18:  rps18 embedding_109 -0.12774389
19:   psbH embedding_857  0.12749911
20:   atpB   embedding_3  0.12612358
21:   matK embedding_655  0.12590059
22:   ndhH embedding_551 -0.12490483
23:   psbA  embedding_49  0.12428956
24:   psaC embedding_749  0.12285694
25:  rpl16 embedding_776 -0.12050352
26:   atpF embedding_954 -0.11730631
27:   petG embedding_646 -0.11638646
28:   psaJ embedding_795 -0.11632859
29:   atpH embedding_708  0.11486318
30:   ndhJ embedding_373 -0.11410104
31:   ndhE embedding_527 -0.11056010
32:  rps11 embedding_397 -0.11048576
33:  rpoC1 embedding_892 -0.11033718
34:   psaA embedding_814  0.10928489
35:   rps3 embedding_812 -0.10875713
36:   atpE embedding_162  0.10726716
37:   ndhG embedding_126 -0.10724990
38:   ndhK embedding_160 -0.10692423
39:   cemA embedding_705  0.10648624
40:  rps14 embedding_702  0.10638836
41:   rps8 embedding_302  0.10524383
42:  rpl23 embedding_533 -0.10513011
43:  rps19  embedding_72 -0.10374081
44:   psbE  embedding_17  0.10259995
45:   rps4 embedding_697 -0.10215691
46:   ndhI embedding_805 -0.10192681
47:   ndhB embedding_548  0.10027748
48:  rpl20 embedding_151 -0.09996291
49:   ycf4 embedding_577 -0.09948303
50:   petN embedding_274 -0.09823504
51:  rpl33 embedding_666  0.09718504
52:   psbN embedding_779 -0.09711985
53:   psbZ embedding_272  0.09588992
54:   ccsA  embedding_91 -0.09554074
55:   rpoB embedding_428  0.09224667
56:   atpI embedding_521 -0.09201985
57:   psbI embedding_669 -0.08516818
58:   ndhA embedding_677  0.08468371
59:   rps7 embedding_855 -0.08223373
60:   ndhD embedding_323 -0.07845143
61:   rpl2 embedding_226  0.07628314
      gene     embedding correlation
                        Estimate Std. Error    t value     Pr(>|t|)
(Intercept)            -32.01903   5.817323  -5.504083 3.890441e-08
petA__embedding_906    358.11997  60.848755   5.885412 4.221440e-09
psaB__embedding_199    773.28216 132.215728   5.848640 5.261375e-09
psaC__embedding_749    576.15047 158.649073   3.631603 2.844032e-04
psbA__embedding_49     928.07031 208.123381   4.459231 8.400147e-06
psbB__embedding_295    945.81181 122.405400   7.726880 1.315950e-14
psbC__embedding_36    1650.30827 134.960111  12.228119 6.466110e-34
psbD__embedding_897    435.63687 208.856350   2.085821 3.704394e-02
psbF__embedding_546     89.88701  10.157892   8.848982 1.194517e-18
psbH__embedding_857    147.40796  46.790558   3.150378 1.639925e-03
psbJ__embedding_72    1222.04552  79.291836  15.411997 1.978669e-52
psbT__embedding_803   -192.32117  47.399567  -4.057446 5.034604e-05
rbcL__embedding_916  -1363.24654 117.343949 -11.617527 8.104530e-31
rpl36__embedding_716    85.88021  40.748567   2.107564 3.511689e-02
ycf3__embedding_381   -232.51183  49.590194  -4.688665 2.820860e-06
matK__embedding_655    166.88642  45.925338   3.633864 2.819280e-04
ndhH__embedding_551   -401.24569  96.453087  -4.160009 3.234126e-05
psbM__embedding_815    319.33678  52.616335   6.069157 1.377722e-09
                       Estimate Std. Error    t value      Pr(>|t|)
(Intercept)            13.19018   1.757609   7.504615  7.212267e-14
psaB__embedding_199   392.31867 124.386105   3.154039  1.619497e-03
psbB__embedding_295  1297.69993 113.388357  11.444737  5.717790e-30
psbJ__embedding_72   1646.35237  72.226317  22.794356 1.125488e-109
psbK__embedding_57    322.55450  42.253311   7.633828  2.696963e-14
rbcL__embedding_916 -1705.13262 119.403629 -14.280409  2.095057e-45
psbM__embedding_815   451.77776  48.463558   9.322010  1.651371e-20

--- stderr ---

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: Matrix
Loaded glmnet 4.1-10
Tree: 6094 tips with phenotype data
Found 915 candidate clades (10-30 tips)
Selected 55 clades, 914 tips (15%)
Train: 5180 | Test: 914 (55 clades)
Running feature selection on training data...
There were 27 warnings (use warnings() to see them)
Selected 61 predictors

=== SELECTED FEATURES (sorted by |correlation|) ===

=== TESTING CORRELATION THRESHOLDS ===
  threshold=0: 61 predictors, R²=0.359
  threshold=0.1: 47 predictors, R²=0.351
  threshold=0.12: 25 predictors, R²=0.321
  threshold=0.14: 9 predictors, R²=0.272
  threshold=0.16: 6 predictors, R²=0.248
  threshold=0.18: 3 predictors, R²=0.221
  threshold=0.2: 2 predictors, R²=0.201

=== MODELS ===
Standard (cor >= 0.12): 25 predictors
Strict (cor >= 0.16): 6 predictors

Standard model R²: 0.321
Strict model R²: 0.248

=== MODEL SUMMARY (standard) ===
Significant predictors (p < 0.05): 17 (excluding intercept)

=== MODEL SUMMARY (strict) ===

=== WITHIN-CLADE EVALUATION (STANDARD: cor >= 0.12) ===
Predictors: 25
Clades evaluated: 55/55
Mean within-clade Spearman: 0.129 (SD: 0.351)
Median within-clade Spearman: 0.183

=== WITHIN-CLADE EVALUATION (STRICT: cor >= 0.16) ===
Predictors: 6
Mean within-clade Spearman: 0.13 (SD: 0.305)
Median within-clade Spearman: 0.11

Standard model t-test (H0: mean=0): p=0.00872
Strict model t-test (H0: mean=0): p=0.00256

Global Spearman (standard): 0.515
Global Spearman (strict): 0.454
Error in eval(jsub, SDenv, parent.frame()) : object 'Family' not found
Calls: merge ... merge.data.table -> is.data.table -> [ -> [.data.table -> eval -> eval
Execution halted

    ```

    ---

    ## Run: 2026-02-10T18:33:34

    **Script**: `src/analyses/genomic_prediction_phylogenetic_cv.r`
    **Interpreter**: `Rscript`
    **Hash**: `62370197c174`

    ### Summary
    This script performs phylogenetic cross-validation to predict a climate phenotype (bio_8_p50, mean temperature of wettest quarter) using protein embeddings from multiple genes. It trains linear regression models on ~85% of phylogenetic tree tips and tests on held-out clades (10-30 tips each, ~15% total), using feature selection that identifies the most correlated embedding dimensions per gene and selects orthogonal dimensions for top-performing genes. The script compares two model variants with different correlation thresholds (0.15 vs 0.20) and evaluates performance using within-clade Spearman correlations to assess the model's ability to predict phenotypic variation within phylogenetically related groups.

    ### Script
    ```
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
# Global prediction plot
p1 <- ggplot(combined_test, aes(x = pred_strict, y = pheno)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = paste0("Predicted vs Observed (Strict: cor >= ", cor_threshold_strict, ")"),
    subtitle = paste0("Global Spearman = ", round(global_cor_strict, 3)),
    x = "Predicted", y = "Observed"
  ) +
  theme_minimal()

# Faceted by top/bottom clades
highlight_clades <- c(top_clades[1:3], bottom_clades[1:3])
p2 <- ggplot(combined_test[clade %in% highlight_clades], 
             aes(x = pred_strict, y = pheno)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~clade, scales = "free") +
  labs(
    title = "Top 3 vs Bottom 3 Clades (strict model)",
    x = "Predicted", y = "Observed"
  ) +
  theme_minimal()

# Within-clade correlation distribution
p3 <- ggplot(clade_results, aes(x = spearman_strict)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = mean(clade_results$spearman_strict, na.rm = TRUE), 
             color = "darkgreen", linetype = "solid") +
  labs(
    title = "Distribution of Within-Clade Spearman Correlations (strict)",
    subtitle = paste0("Mean = ", round(mean(clade_results$spearman_strict, na.rm = TRUE), 3)),
    x = "Spearman r", y = "Count"
  ) +
  theme_minimal()

# Compare standard vs strict
p4 <- ggplot(clade_results, aes(x = spearman_standard, y = spearman_strict)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Standard vs Strict Model (Within-Clade Spearman)",
    x = paste0("Standard (cor >= ", cor_threshold, ")"), 
    y = paste0("Strict (cor >= ", cor_threshold_strict, ")")
  ) +
  theme_minimal()

# Save plots
ggsave("results/phylo_cv_global_pred.png", p1, width = 8, height = 6)
ggsave("results/phylo_cv_clade_facets.png", p2, width = 10, height = 6)
ggsave("results/phylo_cv_spearman_dist.png", p3, width = 8, height = 5)
ggsave("results/phylo_cv_standard_vs_strict.png", p4, width = 6, height = 6)
message("\nPlots saved to results/")

# ---- 10. SAVE ----
results <- list(
  clade_results = clade_results,
  feature_cors = feature_cors,
  threshold_results = threshold_results,
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
  model_strict = fit_strict
)
saveRDS(results, "results/phylo_cv_clade_results.rds")
    ```

    ### Output
    ```
          gene     embedding correlation dim_rank
    <char>        <char>       <num>    <num>
 1:   rbcL embedding_875 -0.26580203        2
 2:   psbB embedding_868  0.25514263        2
 3:   rbcL embedding_916 -0.25250013        1
 4:   psaB embedding_949  0.22577962        2
 5:   psbJ  embedding_72  0.22293014        1
 6:   psaB embedding_707 -0.22253588        3
 7:   rbcL embedding_212  0.21992603        3
 8:   psbK  embedding_53  0.20429816        2
 9:   psbJ embedding_304 -0.20112298        2
10:   psbB  embedding_62 -0.19401704        3
11:   psbB embedding_295  0.19134844        1
12:   psbM embedding_161 -0.18862276        2
13:   psbK embedding_503  0.17795665        3
14:   psaB embedding_199  0.17085064        1
15:   psbK  embedding_57  0.16038676        1
16:   psbM embedding_815  0.16022707        1
17:   psbT embedding_803 -0.15698894        1
18:   petA embedding_906  0.15426775        1
19:   psbD embedding_897  0.15215471        1
20:   psbJ embedding_942  0.15024990        3
21:   atpA embedding_171  0.13902136        1
22:   rpoA embedding_805 -0.13603285        1
23:   psbC  embedding_36  0.13510181        1
24:   psbF embedding_546  0.13503082        1
25:   ycf3 embedding_381 -0.13452911        1
26:  rpl36 embedding_716  0.13243696        1
27:   ndhC embedding_228 -0.12960872        1
28:  rpl14 embedding_388  0.12888506        1
29:  rps18 embedding_109 -0.12774389        1
30:   psbH embedding_857  0.12749911        1
31:   atpB   embedding_3  0.12612358        1
32:   matK embedding_655  0.12590059        1
33:   ndhH embedding_551 -0.12490483        1
34:   psbA  embedding_49  0.12428956        1
35:   psaC embedding_749  0.12285694        1
36:  rpl16 embedding_776 -0.12050352        1
37:   atpF embedding_954 -0.11730631        1
38:   petG embedding_646 -0.11638646        1
39:   psaJ embedding_795 -0.11632859        1
40:   atpH embedding_708  0.11486318        1
41:   ndhJ embedding_373 -0.11410104        1
42:   ndhE embedding_527 -0.11056010        1
43:  rps11 embedding_397 -0.11048576        1
44:  rpoC1 embedding_892 -0.11033718        1
45:   psaA embedding_814  0.10928489        1
46:   rps3 embedding_812 -0.10875713        1
47:   atpE embedding_162  0.10726716        1
48:   ndhG embedding_126 -0.10724990        1
49:   ndhK embedding_160 -0.10692423        1
50:   cemA embedding_705  0.10648624        1
51:  rps14 embedding_702  0.10638836        1
52:   rps8 embedding_302  0.10524383        1
53:  rpl23 embedding_533 -0.10513011        1
54:  rps19  embedding_72 -0.10374081        1
55:   psbE  embedding_17  0.10259995        1
56:   rps4 embedding_697 -0.10215691        1
57:   ndhI embedding_805 -0.10192681        1
58:   psbM embedding_571  0.10096880        3
59:   ndhB embedding_548  0.10027748        1
60:  rpl20 embedding_151 -0.09996291        1
61:   ycf4 embedding_577 -0.09948303        1
62:   petN embedding_274 -0.09823504        1
63:  rpl33 embedding_666  0.09718504        1
64:   psbN embedding_779 -0.09711985        1
65:   psbZ embedding_272  0.09588992        1
66:   ccsA  embedding_91 -0.09554074        1
67:   rpoB embedding_428  0.09224667        1
68:   atpI embedding_521 -0.09201985        1
69:   psbI embedding_669 -0.08516818        1
70:   ndhA embedding_677  0.08468371        1
71:   rps7 embedding_855 -0.08223373        1
72:   ndhD embedding_323 -0.07845143        1
73:   rpl2 embedding_226  0.07628314        1
      gene     embedding correlation dim_rank
                       Estimate Std. Error    t value     Pr(>|t|)
(Intercept)            41.15747   4.615801   8.916648 6.558566e-19
petA__embedding_906   429.77342  56.668715   7.583963 3.951487e-14
psaB__embedding_199   772.48402 130.112054   5.937067 3.091311e-09
psaB__embedding_949   411.65038 139.693507   2.946811 3.225069e-03
psbB__embedding_295  1164.39295 117.116204   9.942202 4.394578e-23
psbB__embedding_868   895.79304 113.982228   7.859059 4.677108e-15
psbB__embedding_62   -560.62249 133.559913  -4.197536 2.743751e-05
psbD__embedding_897   428.77811 206.741896   2.073978 3.813087e-02
psbJ__embedding_72   1027.28681  85.202062  12.057065 4.933214e-33
psbJ__embedding_304  -708.99379 121.369291  -5.841624 5.486079e-09
psbJ__embedding_942   198.07778  66.948471   2.958660 3.103834e-03
psbK__embedding_53    310.69075  55.467559   5.601306 2.237653e-08
psbK__embedding_503   187.50733  48.926197   3.832453 1.283815e-04
rbcL__embedding_916 -1650.02662 126.484179 -13.045320 2.712291e-38
rbcL__embedding_875 -1542.08367 139.265505 -11.072977 3.517266e-28
rbcL__embedding_212  1104.17572 134.529002   8.207715 2.821820e-16
psbM__embedding_815   270.35536  55.904972   4.835981 1.363423e-06
psbM__embedding_161  -186.27687  49.082299  -3.795194 1.492193e-04
                       Estimate Std. Error     t value     Pr(>|t|)
(Intercept)            45.84918   3.718249  12.3308535 1.879050e-34
psaB__embedding_199   637.97243 128.057733   4.9819126 6.501185e-07
psaB__embedding_949   647.03687 134.141962   4.8235232 1.450987e-06
psaB__embedding_707  -210.44044 166.729953  -1.2621634 2.069470e-01
psbB__embedding_295  1245.83461 116.148036  10.7262649 1.456860e-26
psbB__embedding_868  1035.03182 108.470357   9.5420707 2.099966e-21
psbB__embedding_62   -835.32588 130.881696  -6.3822972 1.896991e-10
psbJ__embedding_72   1076.19583  84.250548  12.7737546 8.247699e-37
psbJ__embedding_304  -695.67569 117.987300  -5.8961913 3.956089e-09
psbJ__embedding_942   162.67506  66.621499   2.4417802 1.464837e-02
psbK__embedding_57     22.38940  46.638508   0.4800625 6.312033e-01
psbK__embedding_53    315.07653  53.049261   5.9393199 3.049198e-09
psbK__embedding_503   271.44765  48.142111   5.6384659 1.806867e-08
rbcL__embedding_916 -1784.06248 126.186995 -14.1382437 1.482817e-44
rbcL__embedding_875 -1591.35615 139.755992 -11.3866757 1.098338e-29
rbcL__embedding_212  1204.24837 130.712620   9.2129465 4.516627e-20
     clade     n n_valid spearman_standard spearman_strict pearson_standard
    <char> <int>   <int>             <num>           <num>            <num>
 1:   7195    19      19        0.70000000      0.70701754       0.67722414
 2:   8230    11      11        0.56363636      0.64835456       0.54121355
 3:  11811    11      11        0.50909091      0.59681248       0.34037456
 4:   9660    18      18        0.45740842      0.53175020       0.46241719
 5:   9812    16      16        0.31470588      0.47275456       0.24840373
 6:   7167    18      18        0.54723807      0.44192055       0.53887260
 7:   9327    25      25        0.61088672      0.43469899       0.51282256
 8:  10271    11      11        0.39635638      0.43280295       0.43086104
 9:   9548    11      11        0.46363636      0.41818182       0.30176012
10:  10183    28      28        0.47892720      0.40339354       0.23607572
11:   9513    13      13        0.46556650      0.39724478       0.26723964
12:  12002    11      11        0.29612833      0.39635638       0.36915235
13:   6114    15      15        0.37500000      0.38571429       0.21726798
14:   8604    13      13        0.30219780      0.37912088       0.35084100
15:  11215    22      22        0.26180405      0.35001449       0.14984078
16:   8427    18      18        0.28600933      0.33970061       0.29391461
17:   7725    14      14        0.41538462      0.33186813       0.34358481
18:   7596    29      29        0.26748768      0.32787289       0.28045203
19:   8725    25      25        0.38846154      0.31769231       0.36706845
20:  11944    17      17        0.29328855      0.31547004       0.21471109
21:   7995    19      19        0.21052632      0.28596491       0.41728909
22:   7569    10      10        0.28484848      0.28484848       0.13337141
23:   8130    20      20        0.20601504      0.22255639       0.06445477
24:   7835    29      29        0.28866995      0.21579012       0.34187077
25:   7484    28      28        0.20087575      0.21182266       0.14035535
26:   8574    14      14        0.23516484      0.20439560       0.14435664
27:   8478    13      13        0.27235239      0.20385752      -0.00545877
28:  12112    25      25        0.19119062      0.17734180       0.26537543
29:   8198    11      11        0.17272727      0.17272727       0.17907822
30:  11706    17      17        0.22549020      0.17177927       0.20800535
31:   6395    15      15        0.25714286      0.15357143       0.33451937
32:  12138    11      11        0.06363636      0.14545455       0.13095543
33:   6541    16      16        0.13235294      0.13823529       0.03191594
34:   9429    20      20        0.09924812      0.11729323       0.16647956
35:   7026    12      12       -0.24475524      0.10489510      -0.15146666
36:   8842    18      18        0.17027864      0.10423117       0.12146301
37:   8100    24      24        0.07610350      0.09175908       0.07886741
38:   7384    16      16        0.13823529      0.08529412       0.19867976
39:  12100    12      12        0.02797203      0.06293706       0.04733110
40:  10385    16      16       -0.01764706      0.01470588      -0.13192236
41:   6924    13      13       -0.13204964     -0.05502068      -0.17373966
42:   6942    25      25       -0.18000000     -0.07461538      -0.17111471
43:  10353    10      10        0.20668788     -0.09203147      -0.01621120
44:  10927    21      21       -0.05456317     -0.09223774       0.00529930
45:  10034    15      15       -0.14503186     -0.15577496      -0.17683262
46:   8457    10      10        0.03030303     -0.18787879      -0.32249816
47:  11882    10      10       -0.11550205     -0.20060883      -0.27512843
48:   8075    26      26       -0.20184746     -0.23263776      -0.17541006
49:   6268    12      12       -0.25174825     -0.27272727      -0.19531304
50:   8860    23      23       -0.27032370     -0.29898691      -0.30103856
51:  11690    11      11       -0.12727273     -0.34545455      -0.29553376
52:  10718    15      15       -0.32500000     -0.35924947      -0.42496101
53:  10127    10      10       -0.28484848     -0.47878788      -0.29035241
54:   7132    12      12       -0.59440559     -0.61538462      -0.60446381
55:  11970    10      10       -0.11515152     -0.73333333       0.07185236
     clade     n n_valid spearman_standard spearman_strict pearson_standard

--- stderr ---

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: Matrix
Loaded glmnet 4.1-10
Tree: 6094 tips with phenotype data
Found 915 candidate clades (10-30 tips)
Selected 55 clades, 914 tips (15%)
Train: 5180 | Test: 914 (55 clades)
Running feature selection on training data...
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
There were 50 or more warnings (use warnings() to see the first 50)
Selected 73 predictors total

=== SELECTED FEATURES (sorted by |correlation|) ===

=== TESTING CORRELATION THRESHOLDS ===
  threshold=0: 73 predictors, R²=0.391
  threshold=0.1: 59 predictors, R²=0.385
  threshold=0.15: 21 predictors, R²=0.342
  threshold=0.18: 18 predictors, R²=0.334
  threshold=0.2: 15 predictors, R²=0.328
  threshold=0.22: 12 predictors, R²=0.321
  threshold=0.25: 6 predictors, R²=0.267

=== MODELS ===
Standard (gene cor >= 0.15): 21 predictors from 9 genes
Strict (gene cor >= 0.2): 15 predictors from 5 genes

Standard model R²: 0.342
Strict model R²: 0.328

=== MODEL SUMMARY (standard) ===
Significant predictors (p < 0.05): 17 (excluding intercept)

=== MODEL SUMMARY (strict) ===

=== WITHIN-CLADE EVALUATION (STANDARD: cor >= 0.15) ===
Predictors: 21
Clades evaluated: 55/55
Mean within-clade Spearman: 0.16 (SD: 0.273)
Median within-clade Spearman: 0.207

=== WITHIN-CLADE EVALUATION (STRICT: cor >= 0.2) ===
Predictors: 15
Mean within-clade Spearman: 0.138 (SD: 0.311)
Median within-clade Spearman: 0.177

Standard model t-test (H0: mean=0): p=5.98e-05
Strict model t-test (H0: mean=0): p=0.00176

Global Spearman (standard): 0.532
Global Spearman (strict): 0.517

=== TOP 5 CLADES (strict model) ===

Clade 7195 (n=19, r=0.707):
  Orders: Gentianales
  Pheno range: 8.9 - 27.1

Clade 8230 (n=11, r=0.648):
  Orders: Asterales
  Pheno range: 5.9 - 23.5

Clade 11811 (n=11, r=0.597):
  Orders: Poales
  Pheno range: 26.2 - 29.4

Clade 9660 (n=18, r=0.532):
  Orders: Rosales
  Pheno range: 20 - 27.7

Clade 9812 (n=16, r=0.473):
  Orders: Fagales
  Pheno range: 15.4 - 26.8

=== BOTTOM 5 CLADES (strict model) ===

Clade 11970 (n=10, r=-0.733):
  Orders: Poales
  Pheno range: 9.7 - 28.1

Clade 7132 (n=12, r=-0.615):
  Orders: Gentianales
  Pheno range: 6.6 - 18.5

Clade 10127 (n=10, r=-0.479):
  Orders: Malpighiales
  Pheno range: 18.6 - 27.8

Clade 10718 (n=15, r=-0.359):
  Orders: Sapindales
  Pheno range: 14.6 - 27.5

Clade 11690 (n=11, r=-0.345):
  Orders: Poales
  Pheno range: 20.6 - 28
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'

Plots saved to results/

    ```

    ---

    ## Run: 2026-02-10T19:46:11

    **Script**: `src/analyses/genomic_prediction_phylogenetic_cv.r`
    **Interpreter**: `Rscript`
    **Hash**: `9560d3a94099`

    ### Summary
    This script performs phylogenetically-informed cross-validation to evaluate protein embedding-based prediction of a climate phenotype (mean temperature of the wettest quarter). The analysis uses a phylogenetic tree to select holdout clades (10-30 tips each, ~15% of data) as test sets, ensuring that closely related species are grouped together during evaluation. The method trains linear regression models on gene-specific protein embeddings, using correlation-based feature selection to identify the most predictive embedding dimensions from each gene, with top-performing genes getting multiple orthogonal dimensions. Performance is evaluated by computing within-clade Spearman correlations to assess the model's ability to predict phenotypic variation among closely related species, and includes a grid search over feature selection parameters to optimize the approach.

    ### Script
    ```
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
    ```

    ### Output
    ```
          gene     embedding correlation dim_rank
    <char>        <char>       <num>    <num>
 1:   rbcL embedding_875 -0.26580203        2
 2:   psbB embedding_868  0.25514263        2
 3:   rbcL embedding_916 -0.25250013        1
 4:   psaB embedding_949  0.22577962        2
 5:   psbJ  embedding_72  0.22293014        1
 6:   psaB embedding_707 -0.22253588        3
 7:   rbcL embedding_212  0.21992603        3
 8:   psbK  embedding_53  0.20429816        2
 9:   psbJ embedding_304 -0.20112298        2
10:   psbB  embedding_62 -0.19401704        3
11:   psbB embedding_295  0.19134844        1
12:   psbM embedding_161 -0.18862276        2
13:   psbK embedding_503  0.17795665        3
14:   psaB embedding_199  0.17085064        1
15:   psbK  embedding_57  0.16038676        1
16:   psbM embedding_815  0.16022707        1
17:   psbT embedding_803 -0.15698894        1
18:   petA embedding_906  0.15426775        1
19:   psbD embedding_897  0.15215471        1
20:   psbJ embedding_942  0.15024990        3
21:   atpA embedding_171  0.13902136        1
22:   rpoA embedding_805 -0.13603285        1
23:   psbC  embedding_36  0.13510181        1
24:   psbF embedding_546  0.13503082        1
25:   ycf3 embedding_381 -0.13452911        1
26:  rpl36 embedding_716  0.13243696        1
27:   ndhC embedding_228 -0.12960872        1
28:  rpl14 embedding_388  0.12888506        1
29:  rps18 embedding_109 -0.12774389        1
30:   psbH embedding_857  0.12749911        1
31:   atpB   embedding_3  0.12612358        1
32:   matK embedding_655  0.12590059        1
33:   ndhH embedding_551 -0.12490483        1
34:   psbA  embedding_49  0.12428956        1
35:   psaC embedding_749  0.12285694        1
36:  rpl16 embedding_776 -0.12050352        1
37:   atpF embedding_954 -0.11730631        1
38:   petG embedding_646 -0.11638646        1
39:   psaJ embedding_795 -0.11632859        1
40:   atpH embedding_708  0.11486318        1
41:   ndhJ embedding_373 -0.11410104        1
42:   ndhE embedding_527 -0.11056010        1
43:  rps11 embedding_397 -0.11048576        1
44:  rpoC1 embedding_892 -0.11033718        1
45:   psaA embedding_814  0.10928489        1
46:   rps3 embedding_812 -0.10875713        1
47:   atpE embedding_162  0.10726716        1
48:   ndhG embedding_126 -0.10724990        1
49:   ndhK embedding_160 -0.10692423        1
50:   cemA embedding_705  0.10648624        1
51:  rps14 embedding_702  0.10638836        1
52:   rps8 embedding_302  0.10524383        1
53:  rpl23 embedding_533 -0.10513011        1
54:  rps19  embedding_72 -0.10374081        1
55:   psbE  embedding_17  0.10259995        1
56:   rps4 embedding_697 -0.10215691        1
57:   ndhI embedding_805 -0.10192681        1
58:   psbM embedding_571  0.10096880        3
59:   ndhB embedding_548  0.10027748        1
60:  rpl20 embedding_151 -0.09996291        1
61:   ycf4 embedding_577 -0.09948303        1
62:   petN embedding_274 -0.09823504        1
63:  rpl33 embedding_666  0.09718504        1
64:   psbN embedding_779 -0.09711985        1
65:   psbZ embedding_272  0.09588992        1
66:   ccsA  embedding_91 -0.09554074        1
67:   rpoB embedding_428  0.09224667        1
68:   atpI embedding_521 -0.09201985        1
69:   psbI embedding_669 -0.08516818        1
70:   ndhA embedding_677  0.08468371        1
71:   rps7 embedding_855 -0.08223373        1
72:   ndhD embedding_323 -0.07845143        1
73:   rpl2 embedding_226  0.07628314        1
      gene     embedding correlation dim_rank
                       Estimate Std. Error    t value     Pr(>|t|)
(Intercept)            41.15747   4.615801   8.916648 6.558566e-19
petA__embedding_906   429.77342  56.668715   7.583963 3.951487e-14
psaB__embedding_199   772.48402 130.112054   5.937067 3.091311e-09
psaB__embedding_949   411.65038 139.693507   2.946811 3.225069e-03
psbB__embedding_295  1164.39295 117.116204   9.942202 4.394578e-23
psbB__embedding_868   895.79304 113.982228   7.859059 4.677108e-15
psbB__embedding_62   -560.62249 133.559913  -4.197536 2.743751e-05
psbD__embedding_897   428.77811 206.741896   2.073978 3.813087e-02
psbJ__embedding_72   1027.28681  85.202062  12.057065 4.933214e-33
psbJ__embedding_304  -708.99379 121.369291  -5.841624 5.486079e-09
psbJ__embedding_942   198.07778  66.948471   2.958660 3.103834e-03
psbK__embedding_53    310.69075  55.467559   5.601306 2.237653e-08
psbK__embedding_503   187.50733  48.926197   3.832453 1.283815e-04
rbcL__embedding_916 -1650.02662 126.484179 -13.045320 2.712291e-38
rbcL__embedding_875 -1542.08367 139.265505 -11.072977 3.517266e-28
rbcL__embedding_212  1104.17572 134.529002   8.207715 2.821820e-16
psbM__embedding_815   270.35536  55.904972   4.835981 1.363423e-06
psbM__embedding_161  -186.27687  49.082299  -3.795194 1.492193e-04
                       Estimate Std. Error     t value     Pr(>|t|)
(Intercept)            45.84918   3.718249  12.3308535 1.879050e-34
psaB__embedding_199   637.97243 128.057733   4.9819126 6.501185e-07
psaB__embedding_949   647.03687 134.141962   4.8235232 1.450987e-06
psaB__embedding_707  -210.44044 166.729953  -1.2621634 2.069470e-01
psbB__embedding_295  1245.83461 116.148036  10.7262649 1.456860e-26
psbB__embedding_868  1035.03182 108.470357   9.5420707 2.099966e-21
psbB__embedding_62   -835.32588 130.881696  -6.3822972 1.896991e-10
psbJ__embedding_72   1076.19583  84.250548  12.7737546 8.247699e-37
psbJ__embedding_304  -695.67569 117.987300  -5.8961913 3.956089e-09
psbJ__embedding_942   162.67506  66.621499   2.4417802 1.464837e-02
psbK__embedding_57     22.38940  46.638508   0.4800625 6.312033e-01
psbK__embedding_53    315.07653  53.049261   5.9393199 3.049198e-09
psbK__embedding_503   271.44765  48.142111   5.6384659 1.806867e-08
rbcL__embedding_916 -1784.06248 126.186995 -14.1382437 1.482817e-44
rbcL__embedding_875 -1591.35615 139.755992 -11.3866757 1.098338e-29
rbcL__embedding_212  1204.24837 130.712620   9.2129465 4.516627e-20
     clade     n n_valid spearman_standard spearman_strict pearson_standard
    <char> <int>   <int>             <num>           <num>            <num>
 1:   7195    19      19        0.70000000      0.70701754       0.67722414
 2:   8230    11      11        0.56363636      0.64835456       0.54121355
 3:  11811    11      11        0.50909091      0.59681248       0.34037456
 4:   9660    18      18        0.45740842      0.53175020       0.46241719
 5:   9812    16      16        0.31470588      0.47275456       0.24840373
 6:   7167    18      18        0.54723807      0.44192055       0.53887260
 7:   9327    25      25        0.61088672      0.43469899       0.51282256
 8:  10271    11      11        0.39635638      0.43280295       0.43086104
 9:   9548    11      11        0.46363636      0.41818182       0.30176012
10:  10183    28      28        0.47892720      0.40339354       0.23607572
11:   9513    13      13        0.46556650      0.39724478       0.26723964
12:  12002    11      11        0.29612833      0.39635638       0.36915235
13:   6114    15      15        0.37500000      0.38571429       0.21726798
14:   8604    13      13        0.30219780      0.37912088       0.35084100
15:  11215    22      22        0.26180405      0.35001449       0.14984078
16:   8427    18      18        0.28600933      0.33970061       0.29391461
17:   7725    14      14        0.41538462      0.33186813       0.34358481
18:   7596    29      29        0.26748768      0.32787289       0.28045203
19:   8725    25      25        0.38846154      0.31769231       0.36706845
20:  11944    17      17        0.29328855      0.31547004       0.21471109
21:   7995    19      19        0.21052632      0.28596491       0.41728909
22:   7569    10      10        0.28484848      0.28484848       0.13337141
23:   8130    20      20        0.20601504      0.22255639       0.06445477
24:   7835    29      29        0.28866995      0.21579012       0.34187077
25:   7484    28      28        0.20087575      0.21182266       0.14035535
26:   8574    14      14        0.23516484      0.20439560       0.14435664
27:   8478    13      13        0.27235239      0.20385752      -0.00545877
28:  12112    25      25        0.19119062      0.17734180       0.26537543
29:   8198    11      11        0.17272727      0.17272727       0.17907822
30:  11706    17      17        0.22549020      0.17177927       0.20800535
31:   6395    15      15        0.25714286      0.15357143       0.33451937
32:  12138    11      11        0.06363636      0.14545455       0.13095543
33:   6541    16      16        0.13235294      0.13823529       0.03191594
34:   9429    20      20        0.09924812      0.11729323       0.16647956
35:   7026    12      12       -0.24475524      0.10489510      -0.15146666
36:   8842    18      18        0.17027864      0.10423117       0.12146301
37:   8100    24      24        0.07610350      0.09175908       0.07886741
38:   7384    16      16        0.13823529      0.08529412       0.19867976
39:  12100    12      12        0.02797203      0.06293706       0.04733110
40:  10385    16      16       -0.01764706      0.01470588      -0.13192236
41:   6924    13      13       -0.13204964     -0.05502068      -0.17373966
42:   6942    25      25       -0.18000000     -0.07461538      -0.17111471
43:  10353    10      10        0.20668788     -0.09203147      -0.01621120
44:  10927    21      21       -0.05456317     -0.09223774       0.00529930
45:  10034    15      15       -0.14503186     -0.15577496      -0.17683262
46:   8457    10      10        0.03030303     -0.18787879      -0.32249816
47:  11882    10      10       -0.11550205     -0.20060883      -0.27512843
48:   8075    26      26       -0.20184746     -0.23263776      -0.17541006
49:   6268    12      12       -0.25174825     -0.27272727      -0.19531304
50:   8860    23      23       -0.27032370     -0.29898691      -0.30103856
51:  11690    11      11       -0.12727273     -0.34545455      -0.29553376
52:  10718    15      15       -0.32500000     -0.35924947      -0.42496101
53:  10127    10      10       -0.28484848     -0.47878788      -0.29035241
54:   7132    12      12       -0.59440559     -0.61538462      -0.60446381
55:  11970    10      10       -0.11515152     -0.73333333       0.07185236
     clade     n n_valid spearman_standard spearman_strict pearson_standard
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
11:           6               5          0.20           25 0.3547907
12:          15               5          0.22           25 0.3635672
13:           6               5          0.22           20 0.3453777
14:          10               5          0.22           20 0.3453777
15:           6               5          0.12           49 0.3835891
16:           3               5          0.18           15 0.3380522
17:           3               5          0.20           15 0.3380522
18:           3               5          0.22           15 0.3380522
19:          10               3          0.18           24 0.3501233
20:          15               2          0.12           40 0.3625900
    n_top_genes n_dims_per_gene cor_threshold n_predictors  train_r2
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
11:            0.2064017              0.2269667 2.115408e-06       0.5551149
12:            0.2014256              0.2141236 1.124896e-05       0.5364635
13:            0.1957699              0.1776770 7.070346e-06       0.5369148
14:            0.1957699              0.1776770 7.070346e-06       0.5369148
15:            0.1925685              0.2132353 6.211877e-05       0.5714245
16:            0.1845913              0.1724627 3.219289e-05       0.5010324
17:            0.1845913              0.1724627 3.219289e-05       0.5010324
18:            0.1845913              0.1724627 3.219289e-05       0.5010324
19:            0.1831558              0.2035714 5.024579e-06       0.5297128
20:            0.1810192              0.1970588 4.753290e-05       0.5568569
    mean_within_spearman median_within_spearman         pval global_spearman

--- stderr ---

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: Matrix
Loaded glmnet 4.1-10
Tree: 6094 tips with phenotype data
Found 915 candidate clades (10-30 tips)
Selected 55 clades, 914 tips (15%)
Train: 5180 | Test: 914 (55 clades)
Running feature selection on training data...
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
There were 50 or more warnings (use warnings() to see the first 50)
Selected 73 predictors total

=== SELECTED FEATURES (sorted by |correlation|) ===

=== TESTING CORRELATION THRESHOLDS ===
  threshold=0: 73 predictors, R²=0.391
  threshold=0.1: 59 predictors, R²=0.385
  threshold=0.15: 21 predictors, R²=0.342
  threshold=0.18: 18 predictors, R²=0.334
  threshold=0.2: 15 predictors, R²=0.328
  threshold=0.22: 12 predictors, R²=0.321
  threshold=0.25: 6 predictors, R²=0.267

=== MODELS ===
Standard (gene cor >= 0.15): 21 predictors from 9 genes
Strict (gene cor >= 0.2): 15 predictors from 5 genes

Standard model R²: 0.342
Strict model R²: 0.328

=== MODEL SUMMARY (standard) ===
Significant predictors (p < 0.05): 17 (excluding intercept)

=== MODEL SUMMARY (strict) ===

=== WITHIN-CLADE EVALUATION (STANDARD: cor >= 0.15) ===
Predictors: 21
Clades evaluated: 55/55
Mean within-clade Spearman: 0.16 (SD: 0.273)
Median within-clade Spearman: 0.207

=== WITHIN-CLADE EVALUATION (STRICT: cor >= 0.2) ===
Predictors: 15
Mean within-clade Spearman: 0.138 (SD: 0.311)
Median within-clade Spearman: 0.177

Standard model t-test (H0: mean=0): p=5.98e-05
Strict model t-test (H0: mean=0): p=0.00176

Global Spearman (standard): 0.532
Global Spearman (strict): 0.517

=== TOP 5 CLADES (strict model) ===

Clade 7195 (n=19, r=0.707):
  Orders: Gentianales
  Pheno range: 8.9 - 27.1

Clade 8230 (n=11, r=0.648):
  Orders: Asterales
  Pheno range: 5.9 - 23.5

Clade 11811 (n=11, r=0.597):
  Orders: Poales
  Pheno range: 26.2 - 29.4

Clade 9660 (n=18, r=0.532):
  Orders: Rosales
  Pheno range: 20 - 27.7

Clade 9812 (n=16, r=0.473):
  Orders: Fagales
  Pheno range: 15.4 - 26.8

=== BOTTOM 5 CLADES (strict model) ===

Clade 11970 (n=10, r=-0.733):
  Orders: Poales
  Pheno range: 9.7 - 28.1

Clade 7132 (n=12, r=-0.615):
  Orders: Gentianales
  Pheno range: 6.6 - 18.5

Clade 10127 (n=10, r=-0.479):
  Orders: Malpighiales
  Pheno range: 18.6 - 27.8

Clade 10718 (n=15, r=-0.359):
  Orders: Sapindales
  Pheno range: 14.6 - 27.5

Clade 11690 (n=11, r=-0.345):
  Orders: Poales
  Pheno range: 20.6 - 28
Warning message:
Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
ℹ Please use `linewidth` instead. 
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'

Plots saved to results/

=== GRID SEARCH ===
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
  Completed 10/80 grid points
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
  Completed 20/80 grid points
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
  Completed 30/80 grid points
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
  Completed 40/80 grid points
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
  Completed 50/80 grid points
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
  Completed 60/80 grid points
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
  Completed 70/80 grid points
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
Top 3 genes by correlation: rbcL, psbJ, psbB
Top 6 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM
Top 10 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA
Top 15 genes by correlation: rbcL, psbJ, psbB, psaB, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36
  Completed 80/80 grid points
There were 50 or more warnings (use warnings() to see the first 50)

=== GRID SEARCH RESULTS (sorted by mean within-clade Spearman) ===

=== BEST CONFIGURATION ===
n_top_genes: 15
n_dims_per_gene: 5
cor_threshold: 0.18
n_predictors: 55
Mean within-clade Spearman: 0.229
p-value: 1.45e-06
Global Spearman: 0.555

    ```

    ---

