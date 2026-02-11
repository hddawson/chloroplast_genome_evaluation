    ## Run: 2026-02-11T16:50:44

    **Script**: `src/analyses/genomic_prediction_with_embed_and_sites.r`
    **Interpreter**: `Rscript`
    **Hash**: `bebceba1b5b2`

    ### Summary
    This script combines protein embedding features with GWAS-identified genetic variants to predict a climate-related phenotype (bio_8_p50) in phylogenetically structured data. It loads orthogonal GWAS sites previously identified as significant but uncorrelated with embeddings, creates binary amino acid features from these sites using aligned protein sequences, and merges them with embedding features selected through correlation analysis. The script fits three linear models (embeddings-only, GWAS-only, and combined) and evaluates their performance using within-clade cross-validation on holdout phylogenetic clades, comparing Spearman correlations to assess whether adding GWAS features improves prediction accuracy beyond embeddings alone.

    ### Script
    ```
    library(data.table)
library(ape)
library(arrow)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD PREVIOUS RESULTS ----
message("Loading previous results...")

# GWAS-embedding correlation results
gwas_embed_cor <- readRDS("results/gwas_embedding_correlation_results.rds")
orthogonal_sites <- gwas_embed_cor$orthogonal_sites
message("Orthogonal sites available: ", nrow(orthogonal_sites))
print(table(orthogonal_sites$class))

# CV results for comparison
cv_results <- readRDS("results/phylo_cv_clade_results.rds")
feature_cors <- cv_results$feature_cors
best_config <- cv_results$best_config

message("\nBaseline model performance:")
message("  Mean within-clade Spearman: ", round(mean(cv_results$clade_results$spearman_standard, na.rm = TRUE), 3))

# ---- 2. LOAD DATA ----
message("\nLoading data...")

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

message("Tree: ", Ntip(tree), " tips")

# ---- 3. RECREATE HOLDOUT CLADES ----
set.seed(42)

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

id_to_clade <- data.table(
  ID = unlist(holdout_clades),
  clade = rep(names(holdout_clades), sapply(holdout_clades, length))
)

message("Train: ", length(train_ids), " | Test: ", length(holdout_ids), " (", length(holdout_clades), " clades)")

# ---- 4. CREATE GWAS SITE FEATURES ----
message("\nCreating GWAS site features...")

# Load alignments for orthogonal sites
gwas_genes <- unique(orthogonal_sites$gene)
aln_list <- list()
for (gene in gwas_genes) {
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  if (file.exists(aln_file)) {
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (!is.null(aln)) {
      names(aln) <- sub("\\|.*", "", names(aln))
      aln_list[[gene]] <- as.matrix(aln)
    }
  }
}

# Create binary features for orthogonal sites
# Focus on sig_control (the 90 orthogonal ones)
ortho_control <- orthogonal_sites[class == "sig_control"]
message("Using ", nrow(ortho_control), " orthogonal sig_control sites")

create_site_feature <- function(gene, pos, aln_list, all_ids) {
  if (!gene %in% names(aln_list)) return(NULL)
  aln_mat <- aln_list[[gene]]
  if (pos > ncol(aln_mat)) return(NULL)

  # Get residues
  available_ids <- intersect(all_ids, rownames(aln_mat))
  residues <- aln_mat[available_ids, pos]
  residues[residues == "-"] <- NA

  # Find major allele
  res_table <- table(residues, useNA = "no")
  if (length(res_table) < 2) return(NULL)

  major_allele <- names(res_table)[which.max(res_table)]

  # Binary: 1 = major, 0 = minor
  binary <- as.integer(residues == major_allele)
  binary[is.na(residues)] <- NA

  dt <- data.table(ID = available_ids, value = binary)
  return(dt)
}

# Create features for all orthogonal sites
site_features <- list()
for (i in 1:nrow(ortho_control)) {
  site_id <- ortho_control$site[i]
  gene <- ortho_control$gene[i]
  pos <- as.integer(sub(".*_", "", site_id))

  feat <- create_site_feature(gene, pos, aln_list, tree$tip.label)
  if (!is.null(feat) && sum(!is.na(feat$value)) > 100) {
    setnames(feat, "value", paste0("GWAS_", site_id))
    site_features[[site_id]] <- feat
  }
}

message("Created ", length(site_features), " GWAS site features")

# Merge into wide format
if (length(site_features) > 0) {
  gwas_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), site_features)
  gwas_cols <- setdiff(names(gwas_wide), "ID")
  message("GWAS feature matrix: ", nrow(gwas_wide), " x ", length(gwas_cols))
} else {
  stop("No GWAS features created")
}

# ---- 5. PREPARE EMBEDDING FEATURES ----
message("\nPreparing embedding features...")

# Helper functions from main CV script
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

# Use best config params
n_top_genes <- best_config$n_top_genes
n_dims_per_gene <- best_config$n_dims_per_gene
cor_threshold <- best_config$cor_threshold

train_df <- data[ID %in% train_ids]
train_embeds <- clean_embeds[ID %in% train_ids]
test_df <- data[ID %in% holdout_ids]
test_embeds <- clean_embeds[ID %in% holdout_ids]

# Merge for feature selection
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
top_gene_names <- names(gene_ranks)[1:min(n_top_genes, length(gene_ranks))]

# Select embedding features with orthogonal dims
sel_list <- list()
for (gene in valid_genes) {
  gdt <- merged[Gene == gene]
  cors <- compute_gene_cors(gdt, embed_cols)
  available <- names(cors)[!is.na(cors)]
  if (length(available) == 0) next

  n_select <- if (gene %in% top_gene_names) min(n_dims_per_gene, length(available)) else 1

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
    if (d < n_select && length(available) > 1) {
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

# Combine embedding features
embed_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)
embed_cols_all <- setdiff(names(embed_wide), "ID")

# Filter by threshold
gene_max_cors_dt <- data.table(gene = names(gene_best_cors), max_cor = unlist(gene_best_cors))
keep_genes <- gene_max_cors_dt[max_cor >= cor_threshold]$gene
embed_cols_use <- embed_cols_all[sapply(embed_cols_all, function(x) {
  gene <- sub("__.*", "", x)
  gene %in% keep_genes
})]

message("Embedding features: ", length(embed_cols_use))

# ---- 6. COMBINE FEATURES AND FIT MODELS ----
message("\nFitting models...")

# Merge all features for training
train_features <- merge(embed_wide[ID %in% train_ids], gwas_wide[ID %in% train_ids], by = "ID", all = TRUE)
train_features <- merge(train_features, train_df[, .(ID, pheno = get(pheno_col))], by = "ID")

# Median impute
all_pred_cols <- c(embed_cols_use, gwas_cols)
for (col in all_pred_cols) {
  if (col %in% names(train_features)) {
    med <- median(train_features[[col]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    train_features[is.na(get(col)), (col) := med]
  }
}

# Model 1: Embeddings only (baseline)
formula_embed <- paste("pheno ~", paste(embed_cols_use, collapse = " + "))
fit_embed <- lm(as.formula(formula_embed), data = train_features, na.action = na.exclude)

# Model 2: GWAS only
formula_gwas <- paste("pheno ~", paste(gwas_cols, collapse = " + "))
fit_gwas <- lm(as.formula(formula_gwas), data = train_features, na.action = na.exclude)

# Model 3: Combined
formula_combined <- paste("pheno ~", paste(c(embed_cols_use, gwas_cols), collapse = " + "))
fit_combined <- lm(as.formula(formula_combined), data = train_features, na.action = na.exclude)

message("\nTraining R²:")
message("  Embeddings only: ", round(summary(fit_embed)$r.squared, 3))
message("  GWAS only: ", round(summary(fit_gwas)$r.squared, 3))
message("  Combined: ", round(summary(fit_combined)$r.squared, 3))

# ---- 7. PREPARE TEST DATA ----
message("\nPreparing test data...")

# Merge test embeddings
test_merged <- merge(
  test_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
  test_df[, .(ID, pheno = get(pheno_col))],
  by = "ID"
)

test_embed_list <- list()
for (gene in unique(test_merged$Gene)) {
  gene_data <- test_merged[Gene == gene, c("ID", embed_cols), with = FALSE]
  if (nrow(gene_data) > 0) {
    setnames(gene_data, old = embed_cols, new = paste0(gene, "__", embed_cols))
    test_embed_list[[gene]] <- gene_data
  }
}

test_embed_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), test_embed_list)

# Merge with GWAS features
test_features <- merge(test_embed_wide, gwas_wide[ID %in% holdout_ids], by = "ID", all.x = TRUE)
test_features <- merge(test_features, test_df[, .(ID, pheno = get(pheno_col))], by = "ID")

# Median impute test
for (col in all_pred_cols) {
  if (col %in% names(test_features)) {
    med <- median(test_features[[col]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    test_features[is.na(get(col)), (col) := med]
  }
}

# ---- 8. PREDICT AND EVALUATE ----
message("\nEvaluating on holdout clades...")

test_features$pred_embed <- predict(fit_embed, newdata = test_features)
test_features$pred_gwas <- predict(fit_gwas, newdata = test_features)
test_features$pred_combined <- predict(fit_combined, newdata = test_features)

# Add clade info
test_features <- merge(test_features, id_to_clade, by = "ID")

# Within-clade correlations
clade_results <- test_features[, .(
  n = .N,
  spearman_embed = if (sum(!is.na(pred_embed) & !is.na(pheno)) >= 3)
    cor(pred_embed, pheno, method = "spearman", use = "complete.obs") else NA_real_,
  spearman_gwas = if (sum(!is.na(pred_gwas) & !is.na(pheno)) >= 3)
    cor(pred_gwas, pheno, method = "spearman", use = "complete.obs") else NA_real_,
  spearman_combined = if (sum(!is.na(pred_combined) & !is.na(pheno)) >= 3)
    cor(pred_combined, pheno, method = "spearman", use = "complete.obs") else NA_real_
), by = clade]

# ---- 9. SUMMARY ----
message("\n=== WITHIN-CLADE EVALUATION ===")

message("\nEmbeddings only:")
message("  Mean Spearman: ", round(mean(clade_results$spearman_embed, na.rm = TRUE), 3))
message("  Median Spearman: ", round(median(clade_results$spearman_embed, na.rm = TRUE), 3))
tt_embed <- t.test(clade_results$spearman_embed, mu = 0)
message("  t-test p-value: ", format.pval(tt_embed$p.value, digits = 3))

message("\nGWAS only:")
message("  Mean Spearman: ", round(mean(clade_results$spearman_gwas, na.rm = TRUE), 3))
message("  Median Spearman: ", round(median(clade_results$spearman_gwas, na.rm = TRUE), 3))
tt_gwas <- t.test(clade_results$spearman_gwas, mu = 0)
message("  t-test p-value: ", format.pval(tt_gwas$p.value, digits = 3))

message("\nCombined (Embeddings + GWAS):")
message("  Mean Spearman: ", round(mean(clade_results$spearman_combined, na.rm = TRUE), 3))
message("  Median Spearman: ", round(median(clade_results$spearman_combined, na.rm = TRUE), 3))
tt_combined <- t.test(clade_results$spearman_combined, mu = 0)
message("  t-test p-value: ", format.pval(tt_combined$p.value, digits = 3))

# Paired comparison
message("\n=== PAIRED COMPARISON ===")
diff_combined_embed <- clade_results$spearman_combined - clade_results$spearman_embed
message("Combined - Embed: mean diff = ", round(mean(diff_combined_embed, na.rm = TRUE), 4))
tt_diff <- t.test(diff_combined_embed, mu = 0)
message("  Paired t-test p-value: ", format.pval(tt_diff$p.value, digits = 3))

# Global correlations
global_embed <- cor(test_features$pred_embed, test_features$pheno, method = "spearman", use = "complete.obs")
global_gwas <- cor(test_features$pred_gwas, test_features$pheno, method = "spearman", use = "complete.obs")
global_combined <- cor(test_features$pred_combined, test_features$pheno, method = "spearman", use = "complete.obs")

message("\n=== GLOBAL CORRELATIONS ===")
message("  Embeddings: ", round(global_embed, 3))
message("  GWAS: ", round(global_gwas, 3))
message("  Combined: ", round(global_combined, 3))

# ---- 10. PLOTS ----
message("\nGenerating plots...")

# Within-clade scaled plot for combined model
test_features[, pred_combined_scaled := scale(pred_combined), by = clade]
test_features[, pheno_scaled := scale(pheno), by = clade]

plot_data <- test_features[!is.na(pred_combined_scaled) & !is.na(pheno_scaled)]

p1 <- ggplot(plot_data, aes(x = pred_combined_scaled, y = pheno_scaled)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Within-Clade Prediction: Combined Model (Embeddings + GWAS)",
    subtitle = paste0("Mean within-clade r = ", round(mean(clade_results$spearman_combined, na.rm = TRUE), 3),
                      " | p = ", format.pval(tt_combined$p.value, digits = 2)),
    x = "Predicted (z-scored within clade)",
    y = "Observed (z-scored within clade)"
  ) +
  theme_minimal() +
  coord_fixed()

# Compare models
clade_long <- melt(clade_results, id.vars = c("clade", "n"),
                   measure.vars = c("spearman_embed", "spearman_gwas", "spearman_combined"),
                   variable.name = "model", value.name = "spearman")
clade_long[, model := factor(model, 
                              levels = c("spearman_embed", "spearman_gwas", "spearman_combined"),
                              labels = c("Embeddings", "GWAS", "Combined"))]

p2 <- ggplot(clade_long, aes(x = model, y = spearman, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Within-Clade Spearman by Model",
    x = "", y = "Within-clade Spearman r"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Scatter: embed vs combined
p3 <- ggplot(clade_results, aes(x = spearman_embed, y = spearman_combined)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = "Embeddings vs Combined Model",
    x = "Within-clade r (Embeddings only)",
    y = "Within-clade r (Combined)",
    size = "Clade size"
  ) +
  theme_minimal()

ggsave("results/combined_model_within_clade.png", p1, width = 7, height = 7)
ggsave("results/model_comparison_boxplot.png", p2, width = 6, height = 5)
ggsave("results/embed_vs_combined_scatter.png", p3, width = 7, height = 6)

message("Plots saved to results/")

# ---- 11. SAVE ----
results <- list(
  clade_results = clade_results,
  n_embed_features = length(embed_cols_use),
  n_gwas_features = length(gwas_cols),
  train_r2_embed = summary(fit_embed)$r.squared,
  train_r2_gwas = summary(fit_gwas)$r.squared,
  train_r2_combined = summary(fit_combined)$r.squared,
  global_cor_embed = global_embed,
  global_cor_gwas = global_gwas,
  global_cor_combined = global_combined,
  within_clade_mean_embed = mean(clade_results$spearman_embed, na.rm = TRUE),
  within_clade_mean_gwas = mean(clade_results$spearman_gwas, na.rm = TRUE),
  within_clade_mean_combined = mean(clade_results$spearman_combined, na.rm = TRUE),
  pval_embed = tt_embed$p.value,
  pval_gwas = tt_gwas$p.value,
  pval_combined = tt_combined$p.value,
  model_embed = fit_embed,
  model_gwas = fit_gwas,
  model_combined = fit_combined
)

saveRDS(results, "results/combined_model_results.rds")
message("\nResults saved to results/combined_model_results.rds")

# Print clade results
message("\n=== CLADE RESULTS ===")
print(clade_results[order(-spearman_combined)])
    ```

    ### Output
    ```

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’

The following object is masked from ‘package:arrow’:

    type

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    table, tapply, union, unique, unsplit, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:data.table’:

    first, second

The following object is masked from ‘package:utils’:

    findMatches

The following objects are masked from ‘package:base’:

    expand.grid, I, unname

Loading required package: IRanges

Attaching package: ‘IRanges’

The following object is masked from ‘package:data.table’:

    shift

Loading required package: XVector
Loading required package: GenomeInfoDb

Attaching package: ‘Biostrings’

The following object is masked from ‘package:ape’:

    complement

The following object is masked from ‘package:base’:

    strsplit

Loading previous results...
Orthogonal sites available: 95

   sig_control sig_no_control 
            90              5 

Baseline model performance:
  Mean within-clade Spearman: 0.16

Loading data...
Tree: 6094 tips
Train: 5180 | Test: 914 (55 clades)

Creating GWAS site features...
Using 90 orthogonal sig_control sites
Created 89 GWAS site features
GWAS feature matrix: 6094 x 89

Preparing embedding features...
There were 27 warnings (use warnings() to see them)
There were 27 warnings (use warnings() to see them)
Embedding features: 15

Fitting models...

Training R²:
  Embeddings only: 0.338
  GWAS only: 0.049
  Combined: 0.369

Preparing test data...

Evaluating on holdout clades...
Warning message:
In predict.lm(fit_gwas, newdata = test_features) :
  prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
Warning message:
In predict.lm(fit_combined, newdata = test_features) :
  prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
There were 17 warnings (use warnings() to see them)

=== WITHIN-CLADE EVALUATION ===

Embeddings only:
  Mean Spearman: 0.185
  Median Spearman: 0.172
  t-test p-value: 3.22e-05

GWAS only:
  Mean Spearman: 0.089
  Median Spearman: 0.082
  t-test p-value: 0.0931

Combined (Embeddings + GWAS):
  Mean Spearman: 0.193
  Median Spearman: 0.248
  t-test p-value: 1.21e-05

=== PAIRED COMPARISON ===
Combined - Embed: mean diff = 0.0084
  Paired t-test p-value: 0.493

=== GLOBAL CORRELATIONS ===
  Embeddings: 0.501
  GWAS: 0.153
  Combined: 0.515

Generating plots...
`geom_smooth()` using formula = 'y ~ x'
Warning message:
Removed 17 rows containing non-finite outside the scale range
(`stat_boxplot()`). 
`geom_smooth()` using formula = 'y ~ x'
Plots saved to results/

Results saved to results/combined_model_results.rds

=== CLADE RESULTS ===
     clade     n spearman_embed spearman_gwas spearman_combined
    <char> <int>          <num>         <num>             <num>
 1:   8230    11   0.8009085750            NA       0.800908575
 2:  10183    28   0.6212370005   -0.20674698       0.645320197
 3:   7195    19   0.6263157895   -0.54718845       0.596491228
 4:   8427    18   0.4994840052    0.23979785       0.576009505
 5:   7167    18   0.5761487603    0.34463281       0.569953612
 6:  11811    11   0.4862589963            NA       0.568831279
 7:  10271    11   0.5068545991    0.01740777       0.506854599
 8:   8604    13   0.5219780220            NA       0.505494505
 9:   9812    16   0.4317678503    0.38348249       0.484819676
10:   7569    10   0.4817162723            NA       0.481716272
11:   7026    12   0.5034965035            NA       0.475524476
12:   9548    11   0.4692495090            NA       0.469249509
13:   7596    29   0.3247505458    0.49795500       0.407046966
14:   9513    13   0.3034508738    0.51282259       0.391727492
15:   6114    15   0.4704837582    0.24908280       0.366398716
16:   7384    16   0.3235294118   -0.36407282       0.347058824
17:   7835    29   0.3291045720    0.31622777       0.343392046
18:   9660    18   0.5990708091    0.45648340       0.325919352
19:   7725    14   0.2497275683   -0.14404037       0.316026923
20:   6395    15   0.0535714286    0.30929479       0.310714286
21:   8725    25   0.4362377461    0.22093497       0.294672058
22:  11944    17   0.2858947222            NA       0.285894722
23:  10034    15   0.0017905168    0.43301270       0.266547833
24:   9327    25   0.1800346254    0.13180210       0.256203121
25:   7995    19   0.1826164005    0.43438918       0.254609405
26:   8130    20   0.3699248120   -0.24678382       0.254135338
27:   8198    11   0.0045558205    0.40000000       0.250570126
28:   8574    14   0.3626373626    0.44721360       0.248351648
29:   8478    13   0.2971608194    0.46291005       0.240558759
30:  12138    11   0.1454545455   -0.37267800       0.190909091
31:  11706    17   0.1523346121            NA       0.164619662
32:  11215    22   0.1724627112   -0.05159394       0.164546390
33:  12112    25   0.1350259690            NA       0.144258514
34:   9429    20   0.2527266068   -0.19194297       0.115833028
35:   6541    16   0.1411764706   -0.28697202       0.105882353
36:   6924    13   0.0990372326            NA       0.099037233
37:   7484    28   0.1212868001   -0.23232639       0.081861745
38:  12002    11   0.1369877295   -0.10000000       0.063927607
39:   8100    24  -0.0004349718    0.10542918       0.063070906
40:  12100    12   0.0595447498            NA       0.059544750
41:  10718    15   0.0232350406   -0.08451543       0.044682770
42:  10385    16   0.0264705882    0.02800560       0.032352941
43:   8842    18  -0.0092879257    0.05815182       0.019607843
44:   6268    12  -0.1048951049   -0.16829046       0.006993007
45:   7132    12  -0.0069930070            NA       0.006993007
46:   6942    25   0.0165416429   -0.08498366       0.005770341
47:  10927    21  -0.0207927269    0.02678999      -0.015594545
48:  11882    10  -0.0547114989   -0.52223297      -0.054711499
49:   8075    26  -0.0925690811    0.22666667      -0.102892770
50:  10353    10  -0.1043023303            NA      -0.104302330
51:   8860    23  -0.2643933859            NA      -0.268346932
52:  11690    11  -0.4181818182            NA      -0.390909091
53:  10127    10  -0.2969696970            NA      -0.393939394
54:   8457    10  -0.6319494127            NA      -0.631949413
55:  11970    10  -0.6242424242    0.68334907      -0.636363636
     clade     n spearman_embed spearman_gwas spearman_combined

    ```

    ---

    ## Run: 2026-02-11T17:48:54

    **Script**: `src/analyses/genomic_prediction_with_embed_and_sites.r`
    **Interpreter**: `Rscript`
    **Hash**: `751319824148`

    ### Summary
    This script performs cross-validation to evaluate a combined genomic prediction model that integrates protein sequence embeddings with GWAS-identified amino acid sites. The script uses phylogenetic cross-validation by holding out entire clades, then trains linear models on embedding features (selected using correlation-based feature selection) combined with binary features from 90 orthogonal GWAS sites to predict a climate-related phenotype. It repeats this process across 10 random train-test splits and evaluates performance using within-clade Spearman correlations, comparing the combined embedding+GWAS model against baseline embedding-only results from previous analyses.

    ### Script
    ```
    library(data.table)
library(ape)
library(arrow)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD PREVIOUS RESULTS ----
message("Loading previous results...")

# GWAS-embedding correlation results
gwas_embed_cor <- readRDS("results/gwas_embedding_correlation_results.rds")
orthogonal_sites <- gwas_embed_cor$orthogonal_sites
message("Orthogonal sites available: ", nrow(orthogonal_sites))
print(table(orthogonal_sites$class))

# CV results for comparison
cv_results <- readRDS("results/phylo_cv_clade_results.rds")
feature_cors <- cv_results$feature_cors
best_config <- cv_results$best_config

message("\nBaseline model performance:")
message("  Mean within-clade Spearman: ", round(mean(cv_results$clade_results$spearman_standard, na.rm = TRUE), 3))

# ---- 2. LOAD DATA ----
message("\nLoading data...")

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

message("Tree: ", Ntip(tree), " tips")

# ---- 3. RECREATE HOLDOUT CLADES ----
set.seed(42)
# ---- 3. PRECOMPUTE 10 HOLDOUT SPLITS ----
set.seed(42)

n_reps <- 10
splits <- vector("list", n_reps)

for (r in seq_len(n_reps)) {
  holdout_clades <- select_holdout_clades(tree)
  holdout_ids <- unlist(holdout_clades)
  train_ids <- setdiff(tree$tip.label, holdout_ids)

  id_to_clade <- data.table(
    ID = unlist(holdout_clades),
    clade = rep(names(holdout_clades), sapply(holdout_clades, length))
  )

  splits[[r]] <- list(
    train_ids = train_ids,
    holdout_ids = holdout_ids,
    id_to_clade = id_to_clade
  )
}

# ---- 4. CREATE GWAS SITE FEATURES ----
message("\nCreating GWAS site features...")

# Load alignments for orthogonal sites
gwas_genes <- unique(orthogonal_sites$gene)
aln_list <- list()
for (gene in gwas_genes) {
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  if (file.exists(aln_file)) {
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (!is.null(aln)) {
      names(aln) <- sub("\\|.*", "", names(aln))
      aln_list[[gene]] <- as.matrix(aln)
    }
  }
}

# Create binary features for orthogonal sites
# Focus on sig_control (the 90 orthogonal ones)
ortho_control <- orthogonal_sites[class == "sig_control"]
message("Using ", nrow(ortho_control), " orthogonal sig_control sites")

create_site_feature <- function(gene, pos, aln_list, all_ids) {
  if (!gene %in% names(aln_list)) return(NULL)
  aln_mat <- aln_list[[gene]]
  if (pos > ncol(aln_mat)) return(NULL)

  # Get residues
  available_ids <- intersect(all_ids, rownames(aln_mat))
  residues <- aln_mat[available_ids, pos]
  residues[residues == "-"] <- NA

  # Find major allele
  res_table <- table(residues, useNA = "no")
  if (length(res_table) < 2) return(NULL)

  major_allele <- names(res_table)[which.max(res_table)]

  # Binary: 1 = major, 0 = minor
  binary <- as.integer(residues == major_allele)
  binary[is.na(residues)] <- NA

  dt <- data.table(ID = available_ids, value = binary)
  return(dt)
}

# Create features for all orthogonal sites
site_features <- list()
for (i in 1:nrow(ortho_control)) {
  site_id <- ortho_control$site[i]
  gene <- ortho_control$gene[i]
  pos <- as.integer(sub(".*_", "", site_id))

  feat <- create_site_feature(gene, pos, aln_list, tree$tip.label)
  if (!is.null(feat) && sum(!is.na(feat$value)) > 100) {
    setnames(feat, "value", paste0("GWAS_", site_id))
    site_features[[site_id]] <- feat
  }
}

message("Created ", length(site_features), " GWAS site features")

# Merge into wide format
if (length(site_features) > 0) {
  gwas_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), site_features)
  gwas_cols <- setdiff(names(gwas_wide), "ID")
  message("GWAS feature matrix: ", nrow(gwas_wide), " x ", length(gwas_cols))
} else {
  stop("No GWAS features created")
}
dir.create("results/embed_cv_", showWarnings = FALSE)

all_rep_results <- list()
all_clade_results <- list()

for (r in seq_len(n_reps)) {

  message("\n=== CV Replicate ", r, " ===")

  train_ids <- splits[[r]]$train_ids
  holdout_ids <- splits[[r]]$holdout_ids
  id_to_clade <- splits[[r]]$id_to_clade

  train_df <- data[ID %in% train_ids]
  test_df  <- data[ID %in% holdout_ids]
  train_embeds <- clean_embeds[ID %in% train_ids]
  test_embeds  <- clean_embeds[ID %in% holdout_ids]

  # --- feature selection identical to original ---
  merged <- merge(
    train_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
    train_df[, .(ID, Order, pheno = get(pheno_col))],
    by = "ID"
  )

  gene_counts <- merged[, .N, by = Gene]
  valid_genes <- gene_counts[N >= 5]$Gene

  gene_best_cors <- list()
  for (gene in valid_genes) {
    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) > 0)
      gene_best_cors[[gene]] <- max(abs(cors[available]), na.rm = TRUE)
  }

  gene_ranks <- sort(unlist(gene_best_cors), decreasing = TRUE)
  top_gene_names <- names(gene_ranks)[1:min(n_top_genes, length(gene_ranks))]

  # --- embedding selection identical logic ---
  sel_list <- list()
  for (gene in valid_genes) {

    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) == 0) next

    n_select <- if (gene %in% top_gene_names)
      min(n_dims_per_gene, length(available)) else 1

    embed_matrix <- as.matrix(gdt[, embed_cols, with = FALSE])
    pheno_vec <- gdt$pheno
    selected_dims <- character(0)
    residual_pheno <- pheno_vec

    for (d in seq_len(n_select)) {

      current_cors <- if (d == 1) cors[available] else
        cor(embed_matrix[, available, drop = FALSE],
            residual_pheno, use = "complete.obs")[,1]

      best_dim <- names(sort(abs(current_cors), decreasing = TRUE))[1]
      selected_dims <- c(selected_dims, best_dim)

      if (d < n_select && length(available) > 1) {
        residual_pheno <- residuals(lm(residual_pheno ~ embed_matrix[, best_dim]))
        available <- setdiff(available, best_dim)
      }
    }

    sel <- unique(gdt[, c("ID", selected_dims), with = FALSE])
    setnames(sel, old = selected_dims,
             new = paste0(gene, "__", selected_dims))
    sel_list[[gene]] <- sel
  }

  embed_wide <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), sel_list)
  embed_cols_all <- setdiff(names(embed_wide), "ID")

  gene_max_cors_dt <- data.table(
    gene = names(gene_best_cors),
    max_cor = unlist(gene_best_cors)
  )
  keep_genes <- gene_max_cors_dt[max_cor >= cor_threshold]$gene

  embed_cols_use <- embed_cols_all[sapply(embed_cols_all, function(x) {
    sub("__.*","",x) %in% keep_genes
  })]

  # --- merge train features ---
  train_features <- merge(embed_wide[ID %in% train_ids],
                          gwas_wide[ID %in% train_ids],
                          by="ID", all=TRUE)
  train_features <- merge(train_features,
                          train_df[,.(ID,pheno=get(pheno_col))],
                          by="ID")

  all_pred_cols <- c(embed_cols_use, gwas_cols)

  for (col in all_pred_cols) {
    if (col %in% names(train_features)) {
      med <- median(train_features[[col]], na.rm=TRUE)
      if (is.na(med)) med <- 0
      train_features[is.na(get(col)), (col):=med]
    }
  }

  fit_combined <- lm(
    as.formula(paste("pheno ~",
                     paste(all_pred_cols, collapse=" + "))),
    data=train_features
  )

  # --- prepare test ---
  test_features <- merge(embed_wide[ID %in% holdout_ids],
                         gwas_wide[ID %in% holdout_ids],
                         by="ID", all=TRUE)
  test_features <- merge(test_features,
                         test_df[,.(ID,pheno=get(pheno_col))],
                         by="ID")

  for (col in all_pred_cols) {
    if (col %in% names(test_features)) {
      med <- median(test_features[[col]], na.rm=TRUE)
      if (is.na(med)) med <- 0
      test_features[is.na(get(col)), (col):=med]
    }
  }

  test_features$pred <- predict(fit_combined, newdata=test_features)
  test_features <- merge(test_features, id_to_clade, by="ID")

  clade_results <- test_features[, .(
    n=.N,
    spearman = if (.N>=3)
      cor(pred, pheno, method="spearman", use="complete.obs")
      else NA_real_
  ), by=clade]

  clade_results[, rep := r]
  all_clade_results[[r]] <- clade_results

  all_rep_results[[r]] <- list(
    rep=r,
    train_r2=summary(fit_combined)$r.squared,
    mean_within_clade=mean(clade_results$spearman, na.rm=TRUE),
    global_cor=cor(test_features$pred, test_features$pheno,
                   method="spearman", use="complete.obs")
  )
}
clade_all <- rbindlist(all_clade_results)
rep_summary <- rbindlist(all_rep_results)

message(rep_summary)

saveRDS(list(
  rep_summary=rep_summary,
  clade_results=clade_all
), "results/embed_cv_/cv_results.rds")

p_cv <- ggplot(rep_summary, aes(x=mean_within_clade)) +
  geom_histogram(bins=10, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=mean(rep_summary$mean_within_clade),
             linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="Distribution of Mean Within-Clade Spearman (10 splits)",
       x="Mean within-clade Spearman", y="Count")

ggsave("results/embed_cv_/within_clade_distribution.png",
       p_cv, width=6, height=5)

    ```

    ### Output
    ```

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’

The following object is masked from ‘package:arrow’:

    type

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    table, tapply, union, unique, unsplit, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:data.table’:

    first, second

The following object is masked from ‘package:utils’:

    findMatches

The following objects are masked from ‘package:base’:

    expand.grid, I, unname

Loading required package: IRanges

Attaching package: ‘IRanges’

The following object is masked from ‘package:data.table’:

    shift

Loading required package: XVector
Loading required package: GenomeInfoDb

Attaching package: ‘Biostrings’

The following object is masked from ‘package:ape’:

    complement

The following object is masked from ‘package:base’:

    strsplit

Loading previous results...
Orthogonal sites available: 95

   sig_control sig_no_control 
            90              5 

Baseline model performance:
  Mean within-clade Spearman: 0.16

Loading data...
Tree: 6094 tips
Error in select_holdout_clades(tree) : 
  could not find function "select_holdout_clades"
Execution halted

    ```

    ---

    ## Run: 2026-02-11T18:15:56

    **Script**: `src/analyses/genomic_prediction_with_embed_and_sites.r`
    **Interpreter**: `Rscript`
    **Hash**: `663fafc99435`

    ### Summary
    This script performs cross-validation to evaluate a combined predictive model using protein sequence embeddings and GWAS-identified amino acid sites to predict phenotypic traits. The script loads previous GWAS-embedding correlation results to identify orthogonal control sites, creates binary features from amino acid alignments at these sites, and combines them with pre-selected protein embeddings in a linear regression model. It conducts 10-fold cross-validation using phylogenetic clade holdouts, training on most clades while testing on held-out clades to assess within-clade prediction performance via Spearman correlation. The approach integrates machine learning embeddings with traditional genetic association data to improve phenotype prediction across phylogenetically related groups.

    ### Script
    ```
    library(data.table)
library(ape)
library(arrow)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD PREVIOUS RESULTS ----
message("Loading previous results...")

# GWAS-embedding correlation results
gwas_embed_cor <- readRDS("results/gwas_embedding_correlation_results.rds")
orthogonal_sites <- gwas_embed_cor$orthogonal_sites
message("Orthogonal sites available: ", nrow(orthogonal_sites))
print(table(orthogonal_sites$class))

# CV results for comparison
cv_results <- readRDS("results/phylo_cv_clade_results.rds")
feature_cors <- cv_results$feature_cors
best_config <- cv_results$best_config

message("\nBaseline model performance:")
message("  Mean within-clade Spearman: ", round(mean(cv_results$clade_results$spearman_standard, na.rm = TRUE), 3))

# ---- 2. LOAD DATA ----
message("\nLoading data...")

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

message("Tree: ", Ntip(tree), " tips")

# ---- 3. RECREATE HOLDOUT CLADES ----
set.seed(42)
# ---- 3. PRECOMPUTE 10 HOLDOUT SPLITS ----
set.seed(42)

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

n_reps <- 10
splits <- vector("list", n_reps)

for (r in seq_len(n_reps)) {
  holdout_clades <- select_holdout_clades(tree)
  holdout_ids <- unlist(holdout_clades)
  train_ids <- setdiff(tree$tip.label, holdout_ids)

  id_to_clade <- data.table(
    ID = unlist(holdout_clades),
    clade = rep(names(holdout_clades), sapply(holdout_clades, length))
  )

  splits[[r]] <- list(
    train_ids = train_ids,
    holdout_ids = holdout_ids,
    id_to_clade = id_to_clade
  )
}

# ---- 4. CREATE GWAS SITE FEATURES ----
message("\nCreating GWAS site features...")

# Load alignments for orthogonal sites
gwas_genes <- unique(orthogonal_sites$gene)
aln_list <- list()
for (gene in gwas_genes) {
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  if (file.exists(aln_file)) {
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (!is.null(aln)) {
      names(aln) <- sub("\\|.*", "", names(aln))
      aln_list[[gene]] <- as.matrix(aln)
    }
  }
}

# Create binary features for orthogonal sites
# Focus on sig_control (the 90 orthogonal ones)
ortho_control <- orthogonal_sites[class == "sig_control"]
message("Using ", nrow(ortho_control), " orthogonal sig_control sites")

create_site_feature <- function(gene, pos, aln_list, all_ids) {
  if (!gene %in% names(aln_list)) return(NULL)
  aln_mat <- aln_list[[gene]]
  if (pos > ncol(aln_mat)) return(NULL)

  # Get residues
  available_ids <- intersect(all_ids, rownames(aln_mat))
  residues <- aln_mat[available_ids, pos]
  residues[residues == "-"] <- NA

  # Find major allele
  res_table <- table(residues, useNA = "no")
  if (length(res_table) < 2) return(NULL)

  major_allele <- names(res_table)[which.max(res_table)]

  # Binary: 1 = major, 0 = minor
  binary <- as.integer(residues == major_allele)
  binary[is.na(residues)] <- NA

  dt <- data.table(ID = available_ids, value = binary)
  return(dt)
}

# Create features for all orthogonal sites
site_features <- list()
for (i in 1:nrow(ortho_control)) {
  site_id <- ortho_control$site[i]
  gene <- ortho_control$gene[i]
  pos <- as.integer(sub(".*_", "", site_id))

  feat <- create_site_feature(gene, pos, aln_list, tree$tip.label)
  if (!is.null(feat) && sum(!is.na(feat$value)) > 100) {
    setnames(feat, "value", paste0("GWAS_", site_id))
    site_features[[site_id]] <- feat
  }
}

message("Created ", length(site_features), " GWAS site features")

# Merge into wide format
if (length(site_features) > 0) {
  gwas_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), site_features)
  gwas_cols <- setdiff(names(gwas_wide), "ID")
  message("GWAS feature matrix: ", nrow(gwas_wide), " x ", length(gwas_cols))
} else {
  stop("No GWAS features created")
}
dir.create("results/embed_cv_", showWarnings = FALSE)

all_rep_results <- list()
all_clade_results <- list()

for (r in seq_len(n_reps)) {

  message("\n=== CV Replicate ", r, " ===")

  train_ids <- splits[[r]]$train_ids
  holdout_ids <- splits[[r]]$holdout_ids
  id_to_clade <- splits[[r]]$id_to_clade

  train_df <- data[ID %in% train_ids]
  test_df  <- data[ID %in% holdout_ids]
  train_embeds <- clean_embeds[ID %in% train_ids]
  test_embeds  <- clean_embeds[ID %in% holdout_ids]

  # --- feature selection identical to original ---
  merged <- merge(
    train_embeds[, c("ID", "Gene", embed_cols), with = FALSE],
    train_df[, .(ID, Order, pheno = get(pheno_col))],
    by = "ID"
  )

  gene_counts <- merged[, .N, by = Gene]
  valid_genes <- gene_counts[N >= 5]$Gene

  gene_best_cors <- list()
  for (gene in valid_genes) {
    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) > 0)
      gene_best_cors[[gene]] <- max(abs(cors[available]), na.rm = TRUE)
  }

  gene_ranks <- sort(unlist(gene_best_cors), decreasing = TRUE)
  top_gene_names <- names(gene_ranks)[1:min(n_top_genes, length(gene_ranks))]

  # --- embedding selection identical logic ---
  sel_list <- list()
  for (gene in valid_genes) {

    gdt <- merged[Gene == gene]
    cors <- compute_gene_cors(gdt, embed_cols)
    available <- names(cors)[!is.na(cors)]
    if (length(available) == 0) next

    n_select <- if (gene %in% top_gene_names)
      min(n_dims_per_gene, length(available)) else 1

    embed_matrix <- as.matrix(gdt[, embed_cols, with = FALSE])
    pheno_vec <- gdt$pheno
    selected_dims <- character(0)
    residual_pheno <- pheno_vec

    for (d in seq_len(n_select)) {

      current_cors <- if (d == 1) cors[available] else
        cor(embed_matrix[, available, drop = FALSE],
            residual_pheno, use = "complete.obs")[,1]

      best_dim <- names(sort(abs(current_cors), decreasing = TRUE))[1]
      selected_dims <- c(selected_dims, best_dim)

      if (d < n_select && length(available) > 1) {
        residual_pheno <- residuals(lm(residual_pheno ~ embed_matrix[, best_dim]))
        available <- setdiff(available, best_dim)
      }
    }

    sel <- unique(gdt[, c("ID", selected_dims), with = FALSE])
    setnames(sel, old = selected_dims,
             new = paste0(gene, "__", selected_dims))
    sel_list[[gene]] <- sel
  }

  embed_wide <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), sel_list)
  embed_cols_all <- setdiff(names(embed_wide), "ID")

  gene_max_cors_dt <- data.table(
    gene = names(gene_best_cors),
    max_cor = unlist(gene_best_cors)
  )
  keep_genes <- gene_max_cors_dt[max_cor >= cor_threshold]$gene

  embed_cols_use <- embed_cols_all[sapply(embed_cols_all, function(x) {
    sub("__.*","",x) %in% keep_genes
  })]

  # --- merge train features ---
  train_features <- merge(embed_wide[ID %in% train_ids],
                          gwas_wide[ID %in% train_ids],
                          by="ID", all=TRUE)
  train_features <- merge(train_features,
                          train_df[,.(ID,pheno=get(pheno_col))],
                          by="ID")

  all_pred_cols <- c(embed_cols_use, gwas_cols)

  for (col in all_pred_cols) {
    if (col %in% names(train_features)) {
      med <- median(train_features[[col]], na.rm=TRUE)
      if (is.na(med)) med <- 0
      train_features[is.na(get(col)), (col):=med]
    }
  }

  fit_combined <- lm(
    as.formula(paste("pheno ~",
                     paste(all_pred_cols, collapse=" + "))),
    data=train_features
  )

  # --- prepare test ---
  test_features <- merge(embed_wide[ID %in% holdout_ids],
                         gwas_wide[ID %in% holdout_ids],
                         by="ID", all=TRUE)
  test_features <- merge(test_features,
                         test_df[,.(ID,pheno=get(pheno_col))],
                         by="ID")

  for (col in all_pred_cols) {
    if (col %in% names(test_features)) {
      med <- median(test_features[[col]], na.rm=TRUE)
      if (is.na(med)) med <- 0
      test_features[is.na(get(col)), (col):=med]
    }
  }

  test_features$pred <- predict(fit_combined, newdata=test_features)
  test_features <- merge(test_features, id_to_clade, by="ID")

  clade_results <- test_features[, .(
    n=.N,
    spearman = if (.N>=3)
      cor(pred, pheno, method="spearman", use="complete.obs")
      else NA_real_
  ), by=clade]

  clade_results[, rep := r]
  all_clade_results[[r]] <- clade_results

  all_rep_results[[r]] <- list(
    rep=r,
    train_r2=summary(fit_combined)$r.squared,
    mean_within_clade=mean(clade_results$spearman, na.rm=TRUE),
    global_cor=cor(test_features$pred, test_features$pheno,
                   method="spearman", use="complete.obs")
  )
}
clade_all <- rbindlist(all_clade_results)
rep_summary <- rbindlist(all_rep_results)

message(rep_summary)

saveRDS(list(
  rep_summary=rep_summary,
  clade_results=clade_all
), "results/embed_cv_/cv_results.rds")

p_cv <- ggplot(rep_summary, aes(x=mean_within_clade)) +
  geom_histogram(bins=10, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=mean(rep_summary$mean_within_clade),
             linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="Distribution of Mean Within-Clade Spearman (10 splits)",
       x="Mean within-clade Spearman", y="Count")

ggsave("results/embed_cv_/within_clade_distribution.png",
       p_cv, width=6, height=5)

    ```

    ### Output
    ```

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’

The following object is masked from ‘package:arrow’:

    type

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    table, tapply, union, unique, unsplit, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:data.table’:

    first, second

The following object is masked from ‘package:utils’:

    findMatches

The following objects are masked from ‘package:base’:

    expand.grid, I, unname

Loading required package: IRanges

Attaching package: ‘IRanges’

The following object is masked from ‘package:data.table’:

    shift

Loading required package: XVector
Loading required package: GenomeInfoDb

Attaching package: ‘Biostrings’

The following object is masked from ‘package:ape’:

    complement

The following object is masked from ‘package:base’:

    strsplit

Loading previous results...
Orthogonal sites available: 95

   sig_control sig_no_control 
            90              5 

Baseline model performance:
  Mean within-clade Spearman: 0.16

Loading data...
Tree: 6094 tips

Creating GWAS site features...
Using 90 orthogonal sig_control sites
Created 89 GWAS site features
GWAS feature matrix: 6094 x 89

=== CV Replicate 1 ===
Error in compute_gene_cors(gdt, embed_cols) : 
  could not find function "compute_gene_cors"
Execution halted

    ```

    ---

