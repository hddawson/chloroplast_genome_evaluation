library(arrow)
library(data.table)
library(Biostrings)

# ---------------------------------------------------------------------
# LOAD COMMON DATA
# ---------------------------------------------------------------------

data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")

# Load embedding QC data
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))

# Create clean ID sets per gene (exclude ManualOutlier == TRUE)
setDT(embeds_with_mds)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

# Prepare PC table
pcs_IDS <- aln_index$index
scores <- as.data.table(ev_pcs$x)
setnames(scores, paste0("PC", seq_len(ncol(scores))))
scores[, ID := pcs_IDS]

n_pcs <- 1000L
pc_names <- paste0("PC", seq_len(n_pcs))
scores <- scores[, c("ID", pc_names), with = FALSE]

# Pre-join phenotype + PCs once
pheno_pcs <- data[, .(ID, pheno = get(pheno_col))][scores, on = "ID", nomatch = 0]
stopifnot(nrow(pheno_pcs) > 0)

# ---------------------------------------------------------------------
# PREPARE FILE LISTS
# ---------------------------------------------------------------------

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes_to_process <- get_gene(aln_files)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

# ---------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------

set.seed(123)

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  
  out_path <- paste0("results/residue_models_clean/", gene, "_effects.rds")
  
  if (file.exists(out_path)) {
    message("Skipping ", gene, " (output already exists)")
    next
  }
  
  # Get clean IDs for this gene
  clean_ids_gene <- clean_ids_by_gene[Gene == gene, ID]
  if (length(clean_ids_gene) == 0) {
    message("Skipping ", gene, ": no clean IDs after filtering ManualOutlier")
    next
  }
  
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # --- Read alignment ---
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # --- Strip IDs ---
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  seq_ids <- names(aln)
  
  # Get common IDs: intersection of (alignment IDs, phenotype IDs, clean IDs for this gene)
  common_ids <- intersect(intersect(seq_ids, pheno_pcs$ID), clean_ids_gene)
  if (length(common_ids) < 8500L) {
    message("Skipping ", gene, ": insufficient overlap after filtering outliers (n=", 
            length(common_ids), ")")
    next
  }
  
  # Subset to common IDs
  aln_mat <- aln_mat[match(common_ids, seq_ids), , drop = FALSE]
  pheno_pcs_sub <- pheno_pcs[ID %in% common_ids]
  setkey(pheno_pcs_sub, ID)
  pheno_pcs_sub <- pheno_pcs_sub[match(common_ids, ID)]
  
  y <- pheno_pcs_sub$pheno
  X_pcs <- scale(as.matrix(pheno_pcs_sub[, ..pc_names]))
  
  # -------------------------------------------------------------------
  # LOOP OVER POSITIONS
  # -------------------------------------------------------------------
  
  positions <- seq_len(ncol(aln_mat))
  results_list <- vector("list", length(positions))
  
  for (pos in positions) {
    
    residues <- aln_mat[, pos]
    residue_table <- table(residues)
    
    # Skip gaps-only or monomorphic positions
    if (length(residue_table) < 2L || all(names(residue_table) == "-")) next
    
    # Skip positions with only gaps
    non_gap <- residues != "-"
    if (sum(non_gap) < 8500L) next
    
    # Check special case
    if (pos == 19 & gene == "psbJ") next
    
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    var0 <- which(apply(X_aa, 2, var) == 0)
    if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # -----------------------------------------------------------------
    # FIT MODELS
    # -----------------------------------------------------------------
    
    # Model 0: intercept only (null model)
    fit_null <- tryCatch(lm(y ~ 1), error = function(e) NULL)
    if (is.null(fit_null)) next
    
    # Model 1: residue variation only (no pop structure control)
    fit_aa <- tryCatch(lm(y ~ X_aa), error = function(e) NULL)
    if (is.null(fit_aa)) next
    
    # Model 2: population structure control only
    fit_pcs <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
    if (is.null(fit_pcs)) next
    
    # Model 3: residue variation with population structure control
    fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_aa <- summary(fit_aa)$r.squared
    r2_pcs <- summary(fit_pcs)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    # P-values for each test
    p_aa_only <- tryCatch(anova(fit_null, fit_aa)[2, "Pr(>F)"], error = function(e) NA_real_)
    p_pcs_only <- tryCatch(anova(fit_null, fit_pcs)[2, "Pr(>F)"], error = function(e) NA_real_)
    p_aa_with_pcs <- tryCatch(anova(fit_pcs, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    # Extract only residue effect sizes and SEs from full model
    coef_summary <- summary(fit_full)$coefficients
    aa_coefs <- grep("^X_aa", rownames(coef_summary), value = TRUE)
    
    if (length(aa_coefs) > 0) {
      effects_dt <- data.table(
        Residue = sub("res_factor", "", aa_coefs),
        Effect = coef_summary[aa_coefs, "Estimate"],
        SE = coef_summary[aa_coefs, "Std. Error"],
        T_value = coef_summary[aa_coefs, "t value"],
        P_value = coef_summary[aa_coefs, "Pr(>|t|)"]
      )
    } else {
      effects_dt <- NULL
    }
    
    results_list[[pos]] <- list(
      Gene = gene,
      Aligned_Position = pos,
      N = sum(non_gap),
      R2_aa = r2_aa,
      R2_pcs = r2_pcs,
      R2_full = r2_full,
      P_aa_only = p_aa_only,
      P_pcs_only = p_pcs_only,
      P_aa_with_pcs = p_aa_with_pcs,
      residue_counts = as.list(residue_table),
      effects = effects_dt,
      IDs = common_ids
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) > 0) {
    dir.create("results/residue_models_clean", showWarnings = FALSE, recursive = TRUE)
    saveRDS(results_list, out_path)
    message("\nCompleted ", gene, ": ", length(results_list), " positions with effects")
  } else {
    message("No positions retained for gene ", gene)
  }
}