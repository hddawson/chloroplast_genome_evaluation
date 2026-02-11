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
