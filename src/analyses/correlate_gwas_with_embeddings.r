library(data.table)
library(ape)
library(arrow)
library(Biostrings)
library(ggplot2)
library(pheatmap)

# ---- 1. LOAD GWAS RESULTS AND SITE CLASSIFICATIONS ----
message("Loading GWAS and segregation results...")

seg_results <- readRDS("results/segregation_analysis_results.rds")
site_class <- seg_results$site_class

message("Site classes:")
print(table(site_class$class))

# ---- 2. LOAD EMBEDDING PREDICTORS FROM BEST MODEL ----
message("\nLoading embedding data and model results...")

cv_results <- readRDS("results/phylo_cv_clade_results.rds")
feature_cors <- cv_results$feature_cors

# Get the predictors used in the model
# These are gene__embedding_N format
pred_genes <- unique(feature_cors$gene)
message("Genes in prediction model: ", paste(pred_genes, collapse = ", "))

# Load embeddings
embeds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds)
embed_cols <- grep("embedding", colnames(embeds), value = TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

# ---- 3. LOAD ALIGNMENTS AND CREATE GWAS SITE FEATURES ----
message("\nLoading alignments and creating GWAS site features...")

# Get genes with GWAS sites
gwas_genes <- unique(site_class$Gene)

# Load alignments
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
message("Loaded alignments for ", length(aln_list), " genes")

# Create binary/dummy variables for each GWAS site
# For each site, create indicator for major allele vs others
create_site_features <- function(site_class, aln_list) {
  site_features <- list()
  
  for (i in 1:nrow(site_class)) {
    gene <- site_class$Gene[i]
    pos <- site_class$Position[i]
    site_id <- site_class$site[i]
    site_type <- site_class$class[i]
    
    if (!gene %in% names(aln_list)) next
    aln_mat <- aln_list[[gene]]
    if (pos > ncol(aln_mat)) next
    
    # Get residues at this position
    residues <- aln_mat[, pos]
    residues[residues == "-"] <- NA
    
    # Skip if too few non-gap
    if (sum(!is.na(residues)) < 100) next
    
    # Find major allele
    res_table <- table(residues, useNA = "no")
    if (length(res_table) < 2) next  # monomorphic
    
    major_allele <- names(res_table)[which.max(res_table)]
    
    # Create binary: 1 = major allele, 0 = other
    binary <- as.integer(residues == major_allele)
    binary[is.na(residues)] <- NA
    
    site_features[[site_id]] <- data.table(
      ID = rownames(aln_mat),
      value = binary,
      site = site_id,
      gene = gene,
      position = pos,
      class = site_type,
      major_allele = major_allele,
      maf = 1 - max(res_table) / sum(res_table)
    )
  }
  
  return(site_features)
}

site_features <- create_site_features(site_class, aln_list)
message("Created features for ", length(site_features), " GWAS sites")

# Combine into wide format
site_dt_list <- lapply(names(site_features), function(s) {
  dt <- site_features[[s]][, .(ID, value)]
  setnames(dt, "value", s)
  return(dt)
})

# Merge all site features
site_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), site_dt_list)
message("Site feature matrix: ", nrow(site_wide), " samples x ", ncol(site_wide) - 1, " sites")

# ---- 4. GET EMBEDDING VALUES FOR PREDICTOR GENES ----
message("\nExtracting embedding values for predictor genes...")

# For each gene in the prediction model, get the relevant embedding dimensions
# feature_cors has: gene, embedding, correlation, dim_rank

# Get top embedding dims used
top_features <- feature_cors[dim_rank == 1]  # primary dims per gene
message("Primary embedding dimensions: ", nrow(top_features))

# Extract embedding values for these genes
embed_features <- list()
for (i in 1:nrow(top_features)) {
  gene <- top_features$gene[i]
  embed_dim <- top_features$embedding[i]
  
  gene_embeds <- clean_embeds[Gene == gene, c("ID", embed_dim), with = FALSE]
  if (nrow(gene_embeds) == 0) next
  
  # Rename column to gene__embedding format
  col_name <- paste0(gene, "__", embed_dim)
  setnames(gene_embeds, embed_dim, col_name)
  
  embed_features[[col_name]] <- gene_embeds
}

# Merge embedding features
embed_wide <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), embed_features)
message("Embedding feature matrix: ", nrow(embed_wide), " samples x ", ncol(embed_wide) - 1, " dims")

# ---- 5. MERGE AND COMPUTE CORRELATIONS ----
message("\nComputing correlations between GWAS sites and embeddings...")

# Merge site and embedding features
combined <- merge(site_wide, embed_wide, by = "ID")
message("Combined matrix: ", nrow(combined), " samples")

# Get column names
site_cols <- setdiff(names(site_wide), "ID")
embed_cols_used <- setdiff(names(embed_wide), "ID")

# Compute correlation matrix
site_mat <- as.matrix(combined[, ..site_cols])
embed_mat <- as.matrix(combined[, ..embed_cols_used])

# Pairwise correlations: sites (rows) x embeddings (cols)
cor_mat <- cor(site_mat, embed_mat, use = "pairwise.complete.obs")
message("Correlation matrix: ", nrow(cor_mat), " sites x ", ncol(cor_mat), " embeddings")

# ---- 6. ANALYZE CORRELATIONS ----
message("\n=== CORRELATION ANALYSIS ===")

# Summary stats
cor_vec <- as.vector(cor_mat)
cor_vec <- cor_vec[!is.na(cor_vec)]

message("\nOverall correlation distribution:")
message("  Mean |r|: ", round(mean(abs(cor_vec)), 4))
message("  Median |r|: ", round(median(abs(cor_vec)), 4))
message("  Max |r|: ", round(max(abs(cor_vec)), 4))
message("  % |r| > 0.1: ", round(100 * mean(abs(cor_vec) > 0.1), 1), "%")
message("  % |r| > 0.2: ", round(100 * mean(abs(cor_vec) > 0.2), 1), "%")
message("  % |r| > 0.3: ", round(100 * mean(abs(cor_vec) > 0.3), 1), "%")

# By site class
site_info <- rbindlist(lapply(site_features, function(x) x[1, .(site, gene, class, maf)]))

cor_by_class <- data.table(
  site = rownames(cor_mat),
  max_abs_cor = apply(abs(cor_mat), 1, max, na.rm = TRUE),
  mean_abs_cor = apply(abs(cor_mat), 1, mean, na.rm = TRUE)
)
cor_by_class <- merge(cor_by_class, site_info, by = "site")

message("\nCorrelation by site class:")
print(cor_by_class[, .(
  n_sites = .N,
  mean_max_cor = round(mean(max_abs_cor), 4),
  mean_mean_cor = round(mean(mean_abs_cor), 4),
  pct_max_gt_0.2 = round(100 * mean(max_abs_cor > 0.2), 1)
), by = class])

# ---- 7. FIND ORTHOGONAL GWAS SITES ----
message("\n=== IDENTIFYING ORTHOGONAL GWAS SITES ===")

# Sites with low correlation to all embeddings (potentially useful new predictors)
orthogonal_threshold <- 0.1
cor_by_class[, orthogonal := max_abs_cor < orthogonal_threshold]

message("\nOrthogonal sites (max |r| < ", orthogonal_threshold, "):")
print(cor_by_class[orthogonal == TRUE, .N, by = class])

# Top orthogonal sites by class
message("\nTop orthogonal sig_both sites:")
print(head(cor_by_class[class == "sig_both" & orthogonal == TRUE][order(max_abs_cor)], 10))

message("\nTop orthogonal sig_control sites:")
print(head(cor_by_class[class == "sig_control" & orthogonal == TRUE][order(max_abs_cor)], 10))

# ---- 8. FIND HIGHLY CORRELATED PAIRS ----
message("\n=== HIGHLY CORRELATED SITE-EMBEDDING PAIRS ===")

# Melt correlation matrix to find top pairs
cor_dt <- as.data.table(cor_mat, keep.rownames = "site")
cor_long <- melt(cor_dt, id.vars = "site", variable.name = "embedding", value.name = "correlation")
cor_long <- cor_long[!is.na(correlation)]
cor_long[, abs_cor := abs(correlation)]
cor_long <- merge(cor_long, site_info[, .(site, class, gene)], by = "site")
setnames(cor_long, "gene", "site_gene")
cor_long[, embed_gene := sub("__.*", "", embedding)]

# Top correlated pairs
message("\nTop 20 site-embedding correlations:")
print(cor_long[order(-abs_cor)][1:20])

# Same-gene pairs (site and embedding from same gene)
same_gene <- cor_long[site_gene == embed_gene]
message("\nSame-gene pairs (site gene = embedding gene):")
message("  N pairs: ", nrow(same_gene))
message("  Mean |r|: ", round(mean(same_gene$abs_cor), 4))
message("  Max |r|: ", round(max(same_gene$abs_cor), 4))

diff_gene <- cor_long[site_gene != embed_gene]
message("\nDifferent-gene pairs:")
message("  N pairs: ", nrow(diff_gene))
message("  Mean |r|: ", round(mean(diff_gene$abs_cor), 4))

# ---- 9. PLOTS ----
message("\nGenerating plots...")

# Heatmap of correlations (subset if too large)
n_sites_plot <- min(100, nrow(cor_mat))
n_embeds_plot <- min(20, ncol(cor_mat))

# Select top sites by variance in correlation
site_var <- apply(cor_mat, 1, var, na.rm = TRUE)
top_sites <- names(sort(site_var, decreasing = TRUE))[1:n_sites_plot]

# Add site class annotation
site_annot <- cor_by_class[site %in% top_sites, .(site, class)]
site_annot_df <- data.frame(class = site_annot$class, row.names = site_annot$site)

png("results/gwas_embedding_correlation_heatmap.png", width = 12, height = 10, units = "in", res = 150)
pheatmap(cor_mat[top_sites, , drop = FALSE],
         main = "Correlation: GWAS Sites vs Embedding Dimensions",
         annotation_row = site_annot_df,
         show_rownames = FALSE,
         fontsize = 8)
dev.off()

# Distribution of max correlations by class
p1 <- ggplot(cor_by_class, aes(x = max_abs_cor, fill = class)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~class, ncol = 1) +
  geom_vline(xintercept = orthogonal_threshold, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Max |Correlation| with Any Embedding",
    subtitle = paste0("Red line = orthogonality threshold (", orthogonal_threshold, ")"),
    x = "Max |r| with any embedding dimension",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("results/gwas_embedding_max_cor_dist.png", p1, width = 8, height = 8)

# Scatter: same gene vs different gene correlations
p2 <- ggplot(cor_long, aes(x = site_gene == embed_gene, y = abs_cor)) +
  geom_boxplot(aes(fill = site_gene == embed_gene)) +
  labs(
    title = "Correlation Magnitude: Same vs Different Gene",
    x = "Site and Embedding from Same Gene",
    y = "|Correlation|"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("results/gwas_embedding_same_vs_diff_gene.png", p2, width = 6, height = 5)

# MAF vs correlation
p3 <- ggplot(cor_by_class, aes(x = maf, y = max_abs_cor, color = class)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Minor Allele Frequency vs Embedding Correlation",
    x = "Minor Allele Frequency",
    y = "Max |r| with any embedding"
  ) +
  theme_minimal()

ggsave("results/gwas_embedding_maf_vs_cor.png", p3, width = 8, height = 5)

message("Plots saved to results/")

# ---- 10. SAVE RESULTS ----
results <- list(
  cor_mat = cor_mat,
  cor_by_class = cor_by_class,
  cor_long = cor_long,
  site_info = site_info,
  orthogonal_sites = cor_by_class[orthogonal == TRUE],
  same_gene_cors = same_gene,
  top_pairs = cor_long[order(-abs_cor)][1:100]
)

saveRDS(results, "results/gwas_embedding_correlation_results.rds")
message("\nResults saved to results/gwas_embedding_correlation_results.rds")

# ---- 11. SUMMARY FOR MODELING ----
message("\n=== SUMMARY FOR MODELING ===")

n_orthogonal_both <- sum(cor_by_class$class == "sig_both" & cor_by_class$orthogonal)
n_orthogonal_control <- sum(cor_by_class$class == "sig_control" & cor_by_class$orthogonal)

message("\nOrthogonal GWAS sites that could be added as predictors:")
message("  sig_both: ", n_orthogonal_both, " sites")
message("  sig_control: ", n_orthogonal_control, " sites")
message("  Total: ", n_orthogonal_both + n_orthogonal_control, " potentially independent predictors")

message("\nThese sites have functional evidence (GWAS) but are not captured by embeddings.")
message("Consider adding them to the prediction model as additional features.")
