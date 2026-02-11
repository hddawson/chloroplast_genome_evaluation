library(data.table)
library(Biostrings)
library(ggplot2)

# ---- 1. LOAD RESULTS + STRUCT RUNS ----
results <- fread("results/q3_gwas/q3_run_associations.csv")
struct_runs <- readRDS("results/q3_gwas/struct_runs.rds")

bonf <- 0.05 / nrow(results)
sig_sites <- results[P_corrected < bonf]
message("Significant sites: ", nrow(sig_sites), " (all in: ", 
        paste(unique(sig_sites$Gene), collapse = ", "), ")")

# Focus on rbcL
rbcL_runs <- struct_runs[Gene == "rbcL"]
rbcL_results <- results[Gene == "rbcL"]
rbcL_sig <- sig_sites[Gene == "rbcL"]

# ---- 2. MAP RUNS TO RESIDUE POSITIONS ----
# Each run spans a range of Residue_Index values
run_positions <- rbcL_runs[, .(
  start_pos = min(Residue_Index),
  end_pos = max(Residue_Index),
  n_residues = .N,
  consensus_q3 = consensus_q3[1]
), by = .(site_id, run_id)]

# Add significance info
run_positions <- merge(run_positions, 
                       rbcL_results[, .(site_id, P_corrected, R2_corrected, N)],
                       by = "site_id", all.x = TRUE)
run_positions[, sig := P_corrected < bonf]

setorder(run_positions, start_pos)

message("\n=== All rbcL Q3 runs ===")
message("Total runs: ", nrow(run_positions))
message("Significant runs: ", sum(run_positions$sig, na.rm = TRUE))

message("\n=== Significant runs detail ===")
print(run_positions[sig == TRUE, .(site_id, consensus_q3, start_pos, end_pos, 
                                    n_residues, P_corrected, R2_corrected)])

# ---- 3. LOAD rbcL ALIGNMENT ----
aln_file <- "data/tmp/alignedGenes/rbcL_AA_aligned.fasta"
aln <- readAAStringSet(aln_file)
names(aln) <- sub("\\|.*", "", names(aln))
aln_mat <- as.matrix(aln)
message("\nrbcL alignment: ", nrow(aln_mat), " sequences x ", ncol(aln_mat), " positions")

# ---- 4. CONSENSUS SEQUENCE ----
consensus_seq <- apply(aln_mat, 2, function(col) {
  col <- col[col != "-"]
  if (length(col) == 0) return("-")
  names(which.max(table(col)))
})

# Map alignment columns to ungapped residue index
ungapped_idx <- cumsum(consensus_seq != "-")
ungapped_idx[consensus_seq == "-"] <- NA

message("Consensus length (ungapped): ", max(ungapped_idx, na.rm = TRUE))

# ---- 5. BUILD ANNOTATION TRACK ----
# For each alignment column, assign Q3 class from runs
aln_annot <- data.table(
  aln_col = seq_along(consensus_seq),
  residue = consensus_seq,
  residue_idx = ungapped_idx
)

# Map Residue_Index from runs to annotation
aln_annot <- merge(aln_annot, 
                   rbcL_runs[, .(Residue_Index, consensus_q3, site_id, run_id)],
                   by.x = "residue_idx", by.y = "Residue_Index",
                   all.x = TRUE)

# Add significance
aln_annot <- merge(aln_annot, 
                   rbcL_results[, .(site_id, P_corrected)],
                   by = "site_id", all.x = TRUE)
aln_annot[, sig := !is.na(P_corrected) & P_corrected < bonf]

setorder(aln_annot, aln_col)

# ---- 6. GENE MAP FIGURE ----
# Linear map of rbcL showing Q3 domains + significance

plot_dt <- aln_annot[!is.na(residue_idx)]  # ungapped positions only

p_map <- ggplot(plot_dt, aes(x = residue_idx)) +
  # Q3 track
  geom_tile(aes(y = 1, fill = consensus_q3), height = 0.8, na.rm = TRUE) +
  scale_fill_manual(values = c("C" = "#E41A1C", "E" = "#377EB8", "H" = "#4DAF4A"),
                    name = "Q3 Class", na.value = "grey90") +
  # Significance markers
  geom_point(data = plot_dt[sig == TRUE], 
             aes(y = 2), shape = 25, fill = "black", size = 2) +
  # Significance regions (shaded)
  geom_tile(data = plot_dt[sig == TRUE],
            aes(y = 1), height = 0.8, fill = NA, color = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, max(plot_dt$residue_idx, na.rm = TRUE), by = 50)) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Q3", "Sig"), limits = c(0.5, 2.5)) +
  labs(title = "rbcL: Q3 Secondary Structure Domains and Significant Associations",
       x = "Residue Position", y = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("results/q3_gwas/figures/rbcL_gene_map.png", p_map, width = 14, height = 3, dpi = 300)
ggsave("results/q3_gwas/figures/rbcL_gene_map.pdf", p_map, width = 14, height = 3)

# ---- 7. DETAILED MANHATTAN FOR rbcL ONLY ----
rbcL_results[, neg_log10_P := -log10(P_corrected)]
setorder(rbcL_results, run_id)

# Get start position for each run for x-axis
rbcL_results <- merge(rbcL_results, 
                       run_positions[, .(site_id, start_pos, end_pos)],
                       by = "site_id", all.x = TRUE)

p_rbcL_manh <- ggplot(rbcL_results, aes(x = start_pos, y = neg_log10_P, 
                                          color = consensus_q3)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("C" = "#E41A1C", "E" = "#377EB8", "H" = "#4DAF4A"),
                     name = "Q3 Class") +
  labs(title = "rbcL Q3 Run Associations",
       x = "Residue Position (run start)", y = "-log10(P)") +
  theme_minimal()

ggsave("results/q3_gwas/figures/rbcL_manhattan.png", p_rbcL_manh, width = 12, height = 5, dpi = 300)
ggsave("results/q3_gwas/figures/rbcL_manhattan.pdf", p_rbcL_manh, width = 12, height = 5)

# ---- 8. EXPORT ANNOTATED POSITIONS TABLE ----
sig_detail <- run_positions[sig == TRUE]

# Get the actual consensus residues spanning each significant run
sig_detail[, residues := sapply(seq_len(.N), function(i) {
  idx <- seq(start_pos[i], end_pos[i])
  idx <- idx[idx <= length(consensus_seq) & !is.na(ungapped_idx[idx])]
  # Find alignment columns matching these ungapped positions
  aln_cols <- which(ungapped_idx %in% idx)
  paste(consensus_seq[aln_cols], collapse = "")
})]

fwrite(sig_detail, "results/q3_gwas/rbcL_significant_runs.csv")

message("\n=== Significant rbcL runs with consensus sequences ===")
print(sig_detail[, .(site_id, consensus_q3, start_pos, end_pos, 
                      n_residues, P_corrected, residues)])

message("\nFiles saved to results/q3_gwas/figures/")