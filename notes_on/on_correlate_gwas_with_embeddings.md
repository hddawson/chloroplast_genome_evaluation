    ## Run: 2026-02-11T15:17:49

    **Script**: `src/analyses/correlate_gwas_with_embeddings.r`
    **Interpreter**: `Rscript`
    **Hash**: `fcbef3462321`

    ### Summary
    This script analyzes the correlation between GWAS-identified genomic sites and protein embedding features used in a phylogenetic prediction model. It loads GWAS results and gene alignments, creates binary features for major vs. minor alleles at significant sites, and extracts embedding values for genes in the best-performing prediction model. The script computes correlations between these GWAS site features and embedding dimensions, identifying "orthogonal" GWAS sites that show low correlation with existing embeddings and could serve as independent predictors to improve the model. It generates visualizations including heatmaps and distribution plots, and provides recommendations for incorporating functionally-relevant GWAS sites that aren't captured by the current embedding-based features.

    ### Script
    ```
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

Loading GWAS and segregation results...
Site classes:

      sig_both    sig_control sig_no_control 
           213            360           2076 

Loading embedding data and model results...
Genes in prediction model: rbcL, psbB, psaB, psbJ, psbK, psbM, psbT, petA, psbD, atpA, rpoA, psbC, psbF, ycf3, rpl36, ndhC, rpl14, rps18, psbH, atpB, matK, ndhH, psbA, psaC, rpl16, atpF, petG, psaJ, atpH, ndhJ, ndhE, rps11, rpoC1, psaA, rps3, atpE, ndhG, ndhK, cemA, rps14, rps8, rpl23, rps19, psbE, rps4, ndhI, ndhB, rpl20, ycf4, petN, rpl33, psbN, psbZ, ccsA, rpoB, atpI, psbI, ndhA, rps7, ndhD, rpl2

Loading alignments and creating GWAS site features...
Loaded alignments for 61 genes
Created features for 2649 GWAS sites
Site feature matrix: 10857 samples x 2649 sites

Extracting embedding values for predictor genes...
Primary embedding dimensions: 61
Embedding feature matrix: 10857 samples x 61 dims

Computing correlations between GWAS sites and embeddings...
Combined matrix: 10857 samples
Warning message:
In cor(site_mat, embed_mat, use = "pairwise.complete.obs") :
  the standard deviation is zero
Correlation matrix: 2649 sites x 61 embeddings

=== CORRELATION ANALYSIS ===

Overall correlation distribution:
  Mean |r|: 0.1152
  Median |r|: 0.0838
  Max |r|: 0.8224
  % |r| > 0.1: 43.5%
  % |r| > 0.2: 18%
  % |r| > 0.3: 7.2%

Correlation by site class:
            class n_sites mean_max_cor mean_mean_cor pct_max_gt_0.2
           <char>   <int>        <num>         <num>          <num>
1: sig_no_control    2076       0.3851        0.1231           92.8
2:    sig_control     360       0.2153        0.0642           46.7
3:       sig_both     213       0.3811        0.1235           92.0

=== IDENTIFYING ORTHOGONAL GWAS SITES ===

Orthogonal sites (max |r| < 0.1):
            class     N
           <char> <int>
1:    sig_control    90
2: sig_no_control     5

Top orthogonal sig_both sites:
Key: <site>
Empty data.table (0 rows and 7 cols): site,max_abs_cor,mean_abs_cor,gene,class,maf...

Top orthogonal sig_control sites:
        site max_abs_cor mean_abs_cor   gene       class          maf
      <char>       <num>        <num> <char>      <char>        <num>
 1: psbA_132  0.01918970  0.006399972   psbA sig_control 3.715055e-04
 2: psaB_563  0.02256740  0.008145628   psaB sig_control 1.015791e-03
 3: ndhH_199  0.02277642  0.007560248   ndhH sig_control 7.720517e-04
 4: psaB_254  0.02374839  0.007781156   psaB sig_control 6.464124e-04
 5: psbC_397  0.02539581  0.007072917   psbC sig_control 3.696516e-04
 6:  psaC_28  0.02591526  0.007874094   psaC sig_control 9.219139e-05
 7: psbB_530  0.02719691  0.008632667   psbB sig_control 2.767783e-04
 8: psbB_284  0.02801391  0.009546261   psbB sig_control 2.767783e-04
 9: rpoB_893  0.02869112  0.008748383   rpoB sig_control 8.342603e-04
10: psaA_548  0.02991841  0.006062028   psaA sig_control 1.849968e-04
    orthogonal
        <lgcl>
 1:       TRUE
 2:       TRUE
 3:       TRUE
 4:       TRUE
 5:       TRUE
 6:       TRUE
 7:       TRUE
 8:       TRUE
 9:       TRUE
10:       TRUE

=== HIGHLY CORRELATED SITE-EMBEDDING PAIRS ===

Top 20 site-embedding correlations:
         site            embedding correlation   abs_cor          class
       <char>               <fctr>       <num>     <num>         <char>
 1:   psbJ_17   psbJ__embedding_72   0.8223515 0.8223515 sig_no_control
 2:   psbT_32  psbT__embedding_803  -0.7886037 0.7886037       sig_both
 3:  rps7_124 rpl14__embedding_388  -0.7817063 0.7817063 sig_no_control
 4:   psbT_36  psbT__embedding_803  -0.7796128 0.7796128       sig_both
 5:  rps4_131  rpoB__embedding_428   0.7746512 0.7746512       sig_both
 6:  rps7_187  rpoB__embedding_428   0.7533297 0.7533297    sig_control
 7:    atpA_7 rpl14__embedding_388  -0.7530119 0.7530119 sig_no_control
 8: rpoC1_434 rpl14__embedding_388  -0.7514730 0.7514730 sig_no_control
 9:  atpA_456  rpoB__embedding_428   0.7483890 0.7483890 sig_no_control
10: rpoC1_219  rpoB__embedding_428   0.7444133 0.7444133 sig_no_control
11:  rps7_124  rpoB__embedding_428   0.7439109 0.7439109 sig_no_control
12: rpoC1_416  rpoB__embedding_428   0.7403238 0.7403238 sig_no_control
13:   rps8_49  rpoB__embedding_428   0.7386624 0.7386624 sig_no_control
14:   rpoB_74  rpoB__embedding_428   0.7377157 0.7377157 sig_no_control
15:  rps7_187 rpl14__embedding_388  -0.7374190 0.7374190    sig_control
16:  cemA_434  rpoB__embedding_428   0.7357343 0.7357343    sig_control
17: rpoC1_176  rpoB__embedding_428   0.7352302 0.7352302 sig_no_control
18:  rpoB_517  rpoB__embedding_428   0.7332785 0.7332785 sig_no_control
19: rpoB_1045  rpoB__embedding_428   0.7323765 0.7323765 sig_no_control
20:  rpoB_664 rpl14__embedding_388  -0.7318202 0.7318202 sig_no_control
         site            embedding correlation   abs_cor          class
    site_gene embed_gene
       <char>     <char>
 1:      psbJ       psbJ
 2:      psbT       psbT
 3:      rps7      rpl14
 4:      psbT       psbT
 5:      rps4       rpoB
 6:      rps7       rpoB
 7:      atpA      rpl14
 8:     rpoC1      rpl14
 9:      atpA       rpoB
10:     rpoC1       rpoB
11:      rps7       rpoB
12:     rpoC1       rpoB
13:      rps8       rpoB
14:      rpoB       rpoB
15:      rps7      rpl14
16:      cemA       rpoB
17:     rpoC1       rpoB
18:      rpoB       rpoB
19:      rpoB       rpoB
20:      rpoB      rpl14
    site_gene embed_gene

Same-gene pairs (site gene = embedding gene):
  N pairs: 2649
  Mean |r|: 0.1466
  Max |r|: 0.8224

Different-gene pairs:
  N pairs: 158903
  Mean |r|: 0.1146

Generating plots...
null device 
          1 
`geom_smooth()` using formula = 'y ~ x'
Plots saved to results/

Results saved to results/gwas_embedding_correlation_results.rds

=== SUMMARY FOR MODELING ===

Orthogonal GWAS sites that could be added as predictors:
  sig_both: 0 sites
  sig_control: 90 sites
  Total: 90 potentially independent predictors

These sites have functional evidence (GWAS) but are not captured by embeddings.
Consider adding them to the prediction model as additional features.

    ```

    ---

