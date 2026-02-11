library(arrow)
library(dplyr)
library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD DATA ----
data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

# Population structure PCs (for Frisch-Waugh-Lovell residualization)
ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")

n_pop_pcs <- 1000L
stopifnot(ncol(ev_pcs$x) >= n_pop_pcs)

pc_scores <- as.data.table(ev_pcs$x[, seq_len(n_pop_pcs)])
setnames(pc_scores, paste0("popPC", seq_len(n_pop_pcs)))
pc_scores[, ID := aln_index$index]

# Clean IDs (outlier removal)
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds_with_mds)
clean_ids <- unique(embeds_with_mds[ManualOutlier == FALSE, ID])
stopifnot(length(clean_ids) > 0)
rm(embeds_with_mds); gc()

message("Data rows: ", nrow(data))
message("Pop structure PCs: ", n_pop_pcs)
message("Clean IDs: ", length(clean_ids))

# ---- 2. PRE-RESIDUALIZE PHENOTYPE ----
pheno <- data[, .(ID, pheno = get(pheno_col))]
pop_pc_cols <- paste0("popPC", seq_len(n_pop_pcs))

pheno_pop <- merge(pheno, pc_scores, by = "ID")
pheno_pop <- pheno_pop[ID %in% clean_ids & !is.na(pheno)]
stopifnot(nrow(pheno_pop) > 0)

X_pop <- as.matrix(pheno_pop[, ..pop_pc_cols])
fit_pop <- lm(pheno_pop$pheno ~ X_pop)
pheno_pop[, pheno_resid := residuals(fit_pop)]
rm(X_pop, fit_pop, pc_scores); gc()

pheno_lookup <- pheno_pop[, .(ID, pheno, pheno_resid)]
rm(pheno_pop); gc()
message("Pre-residualized phenotype for ", nrow(pheno_lookup), " individuals")

# ---- 3. LOAD NETSURF ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

# ---- 4. MAP NETSURF TO ALIGNMENT POSITIONS ----
aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

netsurf_mapped <- rbindlist(lapply(aln_files, function(aln_file) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  netsurf_gene <- netsurf_df[Gene == gene]
  if (nrow(netsurf_gene) == 0) return(NULL)
  
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*", "", names(aln))
  common_ids <- intersect(unique(netsurf_gene$ID), names(aln))
  if (length(common_ids) == 0) return(NULL)
  
  rbindlist(lapply(common_ids, function(acc) {
    chars <- strsplit(as.character(aln[[acc]]), "")[[1]]
    pos_map <- cumsum(chars != "-")
    pos_map[chars == "-"] <- NA
    pos_lookup <- setNames(seq_along(chars), pos_map)
    
    ns_acc <- netsurf_gene[ID == acc]
    ns_acc[, Position := as.integer(pos_lookup[as.character(n)])]
    ns_acc[!is.na(Position)]
  }))
}), fill = TRUE)

stopifnot(nrow(netsurf_mapped) > 0)
message("NetSurf mapped: ", nrow(netsurf_mapped), " rows")

# ---- 5. GET CONSENSUS Q3 PER POSITION ----
struct_consensus <- netsurf_mapped[, .(
  consensus_q3 = names(which.max(table(q3))),
  consensus_freq = max(table(q3)) / .N,
  n_seqs = .N
), by = .(Gene, Position)]

struct_consensus <- struct_consensus[consensus_freq >= 0.8]
message("High-consensus positions: ", nrow(struct_consensus))

# ---- 6. MAP TO RESIDUE_INDEX ----
pos_maps <- netsurf_mapped[, .(Aligned_Position = Position, Ungapped_Position = n), by = .(Gene, ID)]
pos_consensus <- pos_maps[, .(Residue_Index = as.integer(median(Ungapped_Position))), by = .(Gene, Aligned_Position)]

struct_consensus <- merge(struct_consensus, pos_consensus,
                          by.x = c("Gene", "Position"),
                          by.y = c("Gene", "Aligned_Position"),
                          all.x = TRUE)

struct_consensus <- struct_consensus[!is.na(Residue_Index)]
setorder(struct_consensus, Gene, Residue_Index, -consensus_freq, Position)
struct_consensus <- struct_consensus[, .SD[1], by = .(Gene, Residue_Index)]
stopifnot(nrow(struct_consensus[, .N, by = .(Gene, Residue_Index)][N > 1]) == 0)
message("Positions with Residue_Index: ", nrow(struct_consensus))

rm(netsurf_df, netsurf_mapped, pos_maps, pos_consensus); gc()

# ---- 7. IDENTIFY Q3 RUNS ----
struct_runs <- unique(struct_consensus[, .(Gene, Residue_Index, consensus_q3)])
setkey(struct_runs, Gene, Residue_Index)

struct_runs[, `:=`(
  prev_q3 = shift(consensus_q3, type = "lag"),
  prev_idx = shift(Residue_Index, type = "lag")
), by = Gene]

struct_runs[, new_run := (consensus_q3 != prev_q3) | 
              is.na(prev_q3) | 
              (Residue_Index - prev_idx > 2), 
            by = Gene]

struct_runs[, run_id := cumsum(new_run), by = Gene]
struct_runs[, site_id := paste0(Gene, "_run", run_id)]

message("Total Q3 runs identified: ", length(unique(struct_runs$site_id)))

# ---- 8. GENE-SPECIFIC PCA + ASSOCIATION TESTING ----
# For each gene: load raw embeddings, subset to Q3 run positions,
# average across positions per taxon, PCA within gene, test top PCs.
# This uses ALL taxa in embeddings (not just tree tips).

n_test_pcs <- 10L
embed_dir <- "data/embeddings"

genes_with_runs <- unique(struct_runs$Gene)
embed_files <- list.files(embed_dir, pattern = "_residue_embeddings\\.parquet$", full.names = TRUE)
embed_genes <- sub("_residue_embeddings\\.parquet$", "", basename(embed_files))
gene_file_map <- setNames(embed_files, embed_genes)

genes_available <- intersect(genes_with_runs, embed_genes)
message("Genes with both runs and embeddings: ", length(genes_available))

results_list <- list()
result_idx <- 0L

for (gene in genes_available) {
  
  gene_runs <- struct_runs[Gene == gene]
  gene_sites <- unique(gene_runs$site_id)
  gene_residues <- unique(gene_runs$Residue_Index)
  
  message("\n--- ", gene, ": ", length(gene_sites), " runs, ", 
          length(gene_residues), " positions ---")
  
  # Load embeddings for this gene's relevant positions only
  embed_file <- gene_file_map[[gene]]
  schema <- schema(open_dataset(embed_file))
  all_cols <- schema$names
  embed_cols <- grep("^embedding_", all_cols, value = TRUE)
  meta_cols <- c("ID", "Residue_Index")
  stopifnot(all(meta_cols %in% all_cols))
  
  gene_dt <- as.data.table(
    open_dataset(embed_file) |>
      filter(Residue_Index %in% gene_residues) |>
      select(all_of(c(meta_cols, embed_cols))) |>
      collect()
  )
  
  if (nrow(gene_dt) == 0) {
    message("  No embedding data, skipping")
    next
  }
  
  # Keep only taxa that have phenotype data
  gene_dt <- gene_dt[ID %in% pheno_lookup$ID]
  
  if (nrow(gene_dt) == 0) {
    message("  No taxa overlap with phenotype, skipping")
    next
  }
  
  message("  Loaded ", nrow(gene_dt), " residue embeddings for ", 
          uniqueN(gene_dt$ID), " taxa")
  
  for (site in gene_sites) {
    run_info <- gene_runs[site_id == site]
    run_residues <- run_info$Residue_Index
    q3_class <- run_info$consensus_q3[1]
    rid <- run_info$run_id[1]
    
    run_dt <- gene_dt[Residue_Index %in% run_residues]
    if (nrow(run_dt) == 0) next
    
    # Average embeddings across positions per taxon
    taxon_embeds <- run_dt[, lapply(.SD, mean, na.rm = TRUE),
                           .SDcols = embed_cols, by = ID]
    
    if (nrow(taxon_embeds) < 100) next
    
    # Gene-specific PCA
    X <- as.matrix(taxon_embeds[, ..embed_cols])
    
    col_var <- apply(X, 2, var, na.rm = TRUE)
    X <- X[, col_var > 0, drop = FALSE]
    if (ncol(X) < n_test_pcs) next
    
    X <- scale(X)
    
    n_pcs_compute <- min(n_test_pcs, ncol(X), nrow(X) - 1)
    svd_res <- tryCatch(svd(X, nu = n_pcs_compute, nv = n_pcs_compute), 
                        error = function(e) NULL)
    if (is.null(svd_res)) next
    
    # PC scores = U * D
    pc_scores_site <- svd_res$u[, seq_len(n_pcs_compute), drop = FALSE] *
      rep(svd_res$d[seq_len(n_pcs_compute)], each = nrow(svd_res$u))
    colnames(pc_scores_site) <- paste0("gPC", seq_len(n_pcs_compute))
    
    site_pheno <- merge(
      data.table(ID = taxon_embeds$ID, pc_scores_site),
      pheno_lookup, by = "ID", all.x = FALSE
    )
    site_pheno <- site_pheno[!is.na(pheno_resid)]
    
    if (nrow(site_pheno) < 100) next
    
    gpc_cols <- paste0("gPC", seq_len(n_pcs_compute))
    pc_mat <- as.matrix(site_pheno[, ..gpc_cols])
    
    if (any(apply(pc_mat, 2, var) == 0)) next
    
    y_resid <- site_pheno$pheno_resid
    y_raw <- site_pheno$pheno
    
    fit_resid <- tryCatch(lm(y_resid ~ pc_mat), error = function(e) NULL)
    fit_resid_null <- lm(y_resid ~ 1)
    fit_raw <- tryCatch(lm(y_raw ~ pc_mat), error = function(e) NULL)
    fit_raw_null <- lm(y_raw ~ 1)
    
    if (is.null(fit_resid) || is.null(fit_raw)) next
    
    p_corrected <- tryCatch(anova(fit_resid_null, fit_resid)[2, "Pr(>F)"], 
                            error = function(e) NA_real_)
    p_uncorrected <- tryCatch(anova(fit_raw_null, fit_raw)[2, "Pr(>F)"], 
                              error = function(e) NA_real_)
    
    pvar <- svd_res$d^2 / sum(svd_res$d^2)
    
    result_idx <- result_idx + 1L
    results_list[[result_idx]] <- data.table(
      site_id = site,
      Gene = gene,
      consensus_q3 = q3_class,
      run_id = rid,
      N = nrow(site_pheno),
      n_pcs_used = n_pcs_compute,
      pvar_pc1 = pvar[1],
      pvar_top10 = sum(pvar[seq_len(min(10, length(pvar)))]),
      R2_corrected = summary(fit_resid)$r.squared,
      R2_uncorrected = summary(fit_raw)$r.squared,
      P_corrected = p_corrected,
      P_uncorrected = p_uncorrected
    )
  }
  
  rm(gene_dt); gc()
  n_tested_gene <- sum(sapply(results_list[max(1, result_idx - length(gene_sites) + 1):result_idx], 
                               function(x) !is.null(x)))
  message("  Tested: ", n_tested_gene, " sites")

}

results <- rbindlist(Filter(Negate(is.null), results_list))
stopifnot(nrow(results) > 0)
message("\nTested ", nrow(results), " sites across ", uniqueN(results$Gene), " genes")

# ---- 9. SAVE ----
dir.create("results/q3_gwas", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "results/q3_gwas/q3_run_associations.rds")
fwrite(results, "results/q3_gwas/q3_run_associations.csv")
saveRDS(struct_runs, "results/q3_gwas/struct_runs.rds")

# ---- 10. VISUALIZE ----
results[, neg_log10_P := -log10(P_corrected)]
results[, neg_log10_P_uncorr := -log10(P_uncorrected)]

setorder(results, Gene, run_id)
results[, x_index := .I]
results[, gene_idx := as.integer(factor(Gene, levels = unique(Gene)))]

bonf_threshold <- -log10(0.05 / nrow(results))
gene_labels <- results[, .(x_mid = median(x_index)), by = Gene]

p_manhattan <- ggplot(results, aes(x = x_index, y = neg_log10_P, color = factor(gene_idx %% 2))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("0" = "grey30", "1" = "steelblue"), guide = "none") +
  scale_x_continuous(breaks = gene_labels$x_mid, labels = gene_labels$Gene) +
  labs(title = paste0("Q3 Run Associations (gene-specific PCA, ", n_pop_pcs, " pop PCs corrected)"),
       x = "Gene", y = "-log10(P)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("results/q3_gwas/manhattan_q3_runs_corrected.png", p_manhattan, width = 14, height = 6, dpi = 150)

p_uncorr <- ggplot(results, aes(x = x_index, y = neg_log10_P_uncorr, color = factor(gene_idx %% 2))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("0" = "grey30", "1" = "steelblue"), guide = "none") +
  scale_x_continuous(breaks = gene_labels$x_mid, labels = gene_labels$Gene) +
  labs(title = "Q3 Run Associations (gene-specific PCA, uncorrected)",
       x = "Gene", y = "-log10(P)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("results/q3_gwas/manhattan_q3_runs_uncorrected.png", p_uncorr, width = 14, height = 6, dpi = 150)

p_q3 <- ggplot(results, aes(x = consensus_q3, y = neg_log10_P, fill = consensus_q3)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2, width = 0.2) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Association by Q3 class (corrected)", x = "Q3", y = "-log10(P)") +
  theme_minimal()

ggsave("results/q3_gwas/boxplot_by_q3.png", p_q3, width = 8, height = 6, dpi = 150)

# ---- 11. SUMMARY ----
message("\n=== SUMMARY ===")
message("Total sites tested: ", nrow(results))
message("Bonferroni threshold: ", round(bonf_threshold, 2))
message("Significant (corrected): ", sum(results$P_corrected < 0.05 / nrow(results), na.rm = TRUE))
message("Significant (uncorrected): ", sum(results$P_uncorrected < 0.05 / nrow(results), na.rm = TRUE))
message("Median R2_corrected: ", round(median(results$R2_corrected, na.rm = TRUE), 4))
message("Median pvar_pc1: ", round(median(results$pvar_pc1, na.rm = TRUE), 4))

message("\nTop 10 corrected hits:")
print(head(results[order(P_corrected)], 10))