library(arrow)
library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD DATA ----
pca <- as.data.table(read_parquet("results/embeddings_full_pca.parquet"))
data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

message("PCA rows: ", nrow(pca))
message("Data rows: ", nrow(data))

# ---- 2. LOAD NETSURF ----
netsurf_files <- list.files("data/netsurf_results/", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(netsurf_files) > 0)

netsurf_df <- rbindlist(lapply(netsurf_files, fread), fill = TRUE)
setnames(netsurf_df, function(x) trimws(gsub("\\[|\\]", ".", x)))
netsurf_df[, ID := sub("^>", "", sub("_Gene_.*", "", id))]
netsurf_df[, Gene := str_split_i(str_split_i(id, "_Gene_", 2), "_Taxonomy_", 1)]

# ---- 3. MAP NETSURF TO ALIGNMENT POSITIONS ----
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

# ---- 4. GET CONSENSUS Q3 PER POSITION ----
struct_consensus <- netsurf_mapped[, .(
  consensus_q3 = names(which.max(table(q3))),
  consensus_freq = max(table(q3)) / .N,
  n_seqs = .N
), by = .(Gene, Position)]

# Filter for high consensus
struct_consensus <- struct_consensus[consensus_freq >= 0.8]
message("High-consensus positions: ", nrow(struct_consensus))

# ---- 5. MAP TO RESIDUE_INDEX ----
pos_maps <- netsurf_mapped[, .(Aligned_Position = Position, Ungapped_Position = n), by = .(Gene, ID)]
pos_consensus <- pos_maps[, .(Residue_Index = as.integer(median(Ungapped_Position))), by = .(Gene, Aligned_Position)]

struct_consensus <- merge(struct_consensus, pos_consensus,
                          by.x = c("Gene", "Position"),
                          by.y = c("Gene", "Aligned_Position"),
                          all.x = TRUE)

struct_consensus <- struct_consensus[!is.na(Residue_Index)]
# Multiple aligned positions can map to the same Residue_Index via median;
# keep the one with highest consensus frequency (ties broken by first Position)
setorder(struct_consensus, Gene, Residue_Index, -consensus_freq, Position)
struct_consensus <- struct_consensus[, .SD[1], by = .(Gene, Residue_Index)]
stopifnot(nrow(struct_consensus[, .N, by = .(Gene, Residue_Index)][N > 1]) == 0)

message("Positions with Residue_Index: ", nrow(struct_consensus))

# ---- 6. IDENTIFY Q3 RUNS (FIXED) ----
# Keep only unique Gene/Residue_Index combinations with their consensus Q3
struct_runs <- unique(struct_consensus[, .(Gene, Residue_Index, consensus_q3)])

# Sort and identify runs
setkey(struct_runs, Gene, Residue_Index)

struct_runs[, `:=`(
  prev_q3 = shift(consensus_q3, type = "lag"),
  prev_idx = shift(Residue_Index, type = "lag")
), by = Gene]

# New run if: different Q3 OR gap > 2 residues
struct_runs[, new_run := (consensus_q3 != prev_q3) | 
              is.na(prev_q3) | 
              (Residue_Index - prev_idx > 2), 
            by = Gene]

struct_runs[, run_id := cumsum(new_run), by = Gene]

# Create unique run identifier
struct_runs[, site_id := paste0(Gene, "_run", run_id)]

message("Total Q3 runs identified: ", length(unique(struct_runs$site_id)))

# ---- 7. MERGE PCA WITH STRUCTURE RUNS (FIXED) ----
pca_struct <- merge(pca, 
                    struct_runs[, .(Gene, Residue_Index, site_id, consensus_q3, run_id)],
                    by = c("Gene", "Residue_Index"),
                    all.x = FALSE)

message("PCA residues in structure runs: ", nrow(pca_struct))

# ---- 8. COMPUTE MEAN PC SCORES PER RUN PER TAXON ----
pc_cols <- grep("^PC[0-9]+$", names(pca_struct), value = TRUE)

run_pcs <- pca_struct[, lapply(.SD, mean, na.rm = TRUE), 
                      .SDcols = pc_cols,
                      by = .(ID, Gene, site_id, consensus_q3, run_id)]

message("Run-aggregated PC scores: ", nrow(run_pcs), " rows")
message("Unique sites: ", length(unique(run_pcs$site_id)))
message("Unique taxa: ", length(unique(run_pcs$ID)))

# ---- 9. ADD PHENOTYPE ----
pheno <- data[, .(ID, pheno = get(pheno_col))]
run_pcs <- merge(run_pcs, pheno, by = "ID", all.x = TRUE)

stopifnot(sum(!is.na(run_pcs$pheno)) > 0)
message("Rows with phenotype: ", sum(!is.na(run_pcs$pheno)))

# ---- 10. TEST ASSOCIATION FOR EACH SITE ----
# Use PC1-PC10 as predictors
test_pcs <- paste0("PC", 1:10)

results_list <- list()

for (site in unique(run_pcs$site_id)) {
  site_data <- run_pcs[site_id == site & !is.na(pheno)]
  
  if (nrow(site_data) < 100) next  # minimum sample size
  
  # Prepare predictors
  X <- as.matrix(site_data[, ..test_pcs])
  y <- site_data$pheno
  
  # Fit models
  fit_null <- lm(y ~ 1)
  fit_full <- lm(y ~ X)
  
  # Extract stats
  r2_null <- summary(fit_null)$r.squared
  r2_full <- summary(fit_full)$r.squared
  p_value <- anova(fit_null, fit_full)[2, "Pr(>F)"]
  
  results_list[[site]] <- data.table(
    site_id = site,
    Gene = site_data$Gene[1],
    consensus_q3 = site_data$consensus_q3[1],
    run_id = site_data$run_id[1],
    N = nrow(site_data),
    R2_null = r2_null,
    R2_full = r2_full,
    Delta_R2 = r2_full - r2_null,
    P = p_value
  )
}

results <- rbindlist(results_list)
message("Tested ", nrow(results), " sites")

# ---- 11. SAVE RESULTS ----
dir.create("results/q3_gwas", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "results/q3_gwas/q3_run_associations.rds")
fwrite(results, "results/q3_gwas/q3_run_associations.csv")

# ---- 12. VISUALIZE ----
results[, neg_log10_P := -log10(P)]

# Manhattan plot
p_manhattan <- ggplot(results, aes(x = seq_along(site_id), y = neg_log10_P, color = consensus_q3)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05 / nrow(results)), linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "Set1", name = "Q3 Structure") +
  labs(title = "Q3 Secondary Structure Run Associations with Temperature",
       x = "Q3 Run Index", y = "-log10(P)") +
  theme_minimal()

ggsave("results/q3_gwas/manhattan_q3_runs.png", p_manhattan, width = 12, height = 6, dpi = 150)

# By Q3 type
p_q3 <- ggplot(results, aes(x = consensus_q3, y = neg_log10_P, fill = consensus_q3)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2, width = 0.2) +
  geom_hline(yintercept = -log10(0.05 / nrow(results)), linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Association strength by Q3 structure type",
       x = "Q3 Structure", y = "-log10(P)") +
  theme_minimal()

ggsave("results/q3_gwas/boxplot_by_q3.png", p_q3, width = 8, height = 6, dpi = 150)

message("\n=== SUMMARY ===")
message("Total sites tested: ", nrow(results))
message("Bonferroni threshold: ", -log10(0.05 / nrow(results)))
message("Significant sites (Bonferroni): ", sum(results$P < 0.05 / nrow(results)))
message("Top 5% sites: ", sum(results$P < quantile(results$P, 0.05)))

print(head(results[order(P)], 10))