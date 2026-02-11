library(data.table)
library(Biostrings)

# ---- 1. LOAD DATA ----
results <- fread("results/q3_gwas/q3_run_associations.csv")
struct_runs <- readRDS("results/q3_gwas/struct_runs.rds")

bonf <- 0.05 / nrow(results)

rbcL_runs <- struct_runs[Gene == "rbcL"]
rbcL_results <- results[Gene == "rbcL"]

# Merge run positions with p-values
rbcL_runs <- merge(rbcL_runs, 
                   rbcL_results[, .(site_id, P_corrected, R2_corrected)],
                   by = "site_id", all.x = TRUE)

# ---- 2. LOAD ALIGNMENT AND BUILD UNGAPPED -> ALIGNED POSITION MAP ----
aln_file <- "data/tmp/alignedGenes/rbcL_AA_aligned.fasta"
aln <- readAAStringSet(aln_file)
names(aln) <- sub("\\|.*", "", names(aln))

# Use first sequence to build a reference ungapped->aligned map
# Better: use consensus (majority non-gap at each column)
aln_mat <- as.matrix(aln)
n_seqs <- nrow(aln_mat)
n_cols <- ncol(aln_mat)

# For Jalview, annotations are per-column of the alignment (1-based)
# We need: Residue_Index (ungapped) -> alignment column
# Use majority: for each alignment column, the ungapped position is
# the most common cumsum(non-gap) value across sequences

# Actually simpler: use consensus sequence to define the mapping
consensus <- apply(aln_mat, 2, function(col) {
  col <- col[col != "-"]
  if (length(col) == 0) return("-")
  names(which.max(table(col)))
})

ungapped_pos <- cumsum(consensus != "-")
ungapped_pos[consensus == "-"] <- NA

# Build lookup: Residue_Index -> alignment column(s)
res_to_aln <- data.table(
  aln_col = seq_len(n_cols),
  residue_idx = ungapped_pos
)
res_to_aln <- res_to_aln[!is.na(residue_idx)]

# ---- 3. MAP RUNS TO ALIGNMENT COLUMNS ----
run_aln <- merge(rbcL_runs[, .(site_id, Residue_Index, consensus_q3, P_corrected, R2_corrected)],
                 res_to_aln, by.x = "Residue_Index", by.y = "residue_idx",
                 all.x = TRUE)
run_aln <- run_aln[!is.na(aln_col)]

# Get run extents in alignment coordinates
run_extents <- run_aln[, .(
  aln_start = min(aln_col),
  aln_end = max(aln_col),
  consensus_q3 = consensus_q3[1],
  P_corrected = P_corrected[1],
  R2_corrected = R2_corrected[1]
), by = site_id]

run_extents[, sig := P_corrected < bonf]
run_extents[, neg_log10_P := -log10(P_corrected)]

message("rbcL runs mapped to alignment: ", nrow(run_extents))
message("Significant: ", sum(run_extents$sig, na.rm = TRUE))

# ---- 4. WRITE JALVIEW FEATURES FILE ----
# Format: https://www.jalview.org/help/html/features/featuresFormat.html
# Group features by type, then list as: description\tSEQID\t-1\tstart\tend\ttype\tscore

# We annotate the FIRST sequence as the reference (or use "ALIGNMENT" for annotation rows)
# Jalview features format:
# Line 1+: feature color definitions
# Then: description<tab>sequenceid<tab>sequenceIndex<tab>start<tab>end<tab>featuretype[<tab>score]

ref_id <- names(aln)[1]  # use first sequence as reference

features_file <- "results/q3_gwas/rbcL_q3_runs.features"

lines <- c(
  # Color definitions
  "Q3_Coil\tcc0000",
  "Q3_Strand\t0000cc",
  "Q3_Helix\t00cc00",
  "Significant_Run\tff00ff",
  # Now the features - Q3 runs
  paste0("STARTGROUP\tQ3_Structure")
)

for (i in seq_len(nrow(run_extents))) {
  r <- run_extents[i]
  q3_type <- switch(r$consensus_q3, "C" = "Q3_Coil", "E" = "Q3_Strand", "H" = "Q3_Helix")
  desc <- paste0(r$site_id, " (", r$consensus_q3, ", p=", signif(r$P_corrected, 3), ")")
  lines <- c(lines, paste(desc, ref_id, -1, r$aln_start, r$aln_end, q3_type, 
                           round(r$neg_log10_P, 2), sep = "\t"))
}

lines <- c(lines, "ENDGROUP\tQ3_Structure")

# Significant runs as separate group
lines <- c(lines, "STARTGROUP\tSignificant")
for (i in which(run_extents$sig)) {
  r <- run_extents[i]
  desc <- paste0(r$site_id, " p=", signif(r$P_corrected, 3))
  lines <- c(lines, paste(desc, ref_id, -1, r$aln_start, r$aln_end, "Significant_Run",
                           round(r$neg_log10_P, 2), sep = "\t"))
}
lines <- c(lines, "ENDGROUP\tSignificant")

writeLines(lines, features_file)
message("Jalview features file: ", features_file)

# ---- 5. WRITE JALVIEW ANNOTATION FILE ----
# This adds annotation rows below the alignment in Jalview
# Format: JALVIEW_ANNOTATION
# Then: NO_GRAPH or BAR_GRAPH or LINE_GRAPH lines

annot_file <- "results/q3_gwas/rbcL_q3_runs.annotations"

annot_lines <- c("JALVIEW_ANNOTATION")

# Row 1: Q3 consensus (character annotation)
q3_per_col <- rep(NA_character_, n_cols)
for (i in seq_len(nrow(run_aln))) {
  q3_per_col[run_aln$aln_col[i]] <- run_aln$consensus_q3[i]
}
# Format: each column separated by |, empty columns are space
q3_vals <- ifelse(is.na(q3_per_col), " ", q3_per_col)
annot_lines <- c(annot_lines, 
                 paste0("NO_GRAPH\tQ3_Structure\tQ3 consensus\t", 
                        paste(q3_vals, collapse = "|")))

# Row 2: -log10(P) as bar graph
p_per_col <- rep(0, n_cols)
for (i in seq_len(nrow(run_aln))) {
  col <- run_aln$aln_col[i]
  p_val <- -log10(run_aln$P_corrected[i])
  if (!is.na(p_val) && p_val > p_per_col[col]) {
    p_per_col[col] <- p_val
  }
}
p_vals <- ifelse(p_per_col == 0, " ", round(p_per_col, 2))
annot_lines <- c(annot_lines, 
                 paste0("BAR_GRAPH\t-log10(P)\tAssociation score\t",
                        paste(p_vals, collapse = "|")))

# Row 3: R2 as line graph
r2_per_col <- rep(0, n_cols)
for (i in seq_len(nrow(run_aln))) {
  col <- run_aln$aln_col[i]
  r2_val <- run_aln$R2_corrected[i]
  if (!is.na(r2_val) && r2_val > r2_per_col[col]) {
    r2_per_col[col] <- r2_val
  }
}
r2_vals <- ifelse(r2_per_col == 0, " ", round(r2_per_col, 4))
annot_lines <- c(annot_lines,
                 paste0("LINE_GRAPH\tR2_corrected\tEffect size\t",
                        paste(r2_vals, collapse = "|")))

# Row 4: Significance binary
sig_per_col <- rep(" ", n_cols)
sig_run_cols <- run_aln[!is.na(P_corrected) & P_corrected < bonf, aln_col]
sig_per_col[sig_run_cols] <- "*"
annot_lines <- c(annot_lines,
                 paste0("NO_GRAPH\tSignificant\tBonferroni sig\t",
                        paste(sig_per_col, collapse = "|")))

writeLines(annot_lines, annot_file)
message("Jalview annotation file: ", annot_file)

# ---- 6. SUMMARY ----
message("\nTo use in Jalview:")
message("  1. Open: ", aln_file)
message("  2. File > Load Features: ", features_file)
message("  3. File > Load Annotations: ", annot_file)
message("\nThe -log10(P) bar graph and Q3 track will appear below the alignment.")