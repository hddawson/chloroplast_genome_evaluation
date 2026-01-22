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

embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))

setDT(embeds_with_mds)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

# Prepare PC table - keep ALL PCs
pcs_IDS <- aln_index$index
scores <- as.data.table(ev_pcs$x)
max_pcs <- ncol(scores)
setnames(scores, paste0("PC", seq_len(max_pcs)))
scores[, ID := pcs_IDS]

message("Total PCs available: ", max_pcs)

# Pre-join phenotype + all PCs
all_pc_names <- paste0("PC", seq_len(max_pcs))
pheno_pcs_full <- data[, .(ID, pheno = get(pheno_col))][scores, on = "ID", nomatch = 0]
stopifnot(nrow(pheno_pcs_full) > 0)

# ---------------------------------------------------------------------
# PREPARE FILE LISTS
# ---------------------------------------------------------------------

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)
get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes_to_process <- get_gene(aln_files)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes")

# ---------------------------------------------------------------------
# N_PCS EXPERIMENT SETUP
# ---------------------------------------------------------------------

# Define n_pcs values: 0, 250, 500, ..., up to max_pcs
n_pcs_values <- seq(0, max_pcs, by = 250)
if (max(n_pcs_values) < max_pcs) {
  n_pcs_values <- c(n_pcs_values, max_pcs)
}

message("Testing ", length(n_pcs_values), " n_pcs values: ", 
        paste(range(n_pcs_values), collapse = " to "))

# Create output directory structure
base_out_dir <- "results/npc_experiment"
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# ---------------------------------------------------------------------
# MAIN LOOP: iterate over n_pcs values
# ---------------------------------------------------------------------

for (n_pcs in n_pcs_values) {
  
  message("\n========== n_pcs = ", n_pcs, " ==========")
  
  out_dir <- file.path(base_out_dir, paste0("npc_", n_pcs))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Subset PCs for this run
  if (n_pcs == 0) {
    pc_names_use <- character(0)
  } else {
    pc_names_use <- paste0("PC", seq_len(n_pcs))
  }
  
  for (gene in genes_to_process) {
    
    out_path <- file.path(out_dir, paste0(gene, "_effects.rds"))
    
    if (file.exists(out_path)) {
      message("Skipping ", gene, " (output already exists)")
      next
    }
    
    # Get clean IDs for this gene
    clean_ids_gene <- clean_ids_by_gene[Gene == gene, ID]
    if (length(clean_ids_gene) == 0) {
      message("Skipping ", gene, ": no clean IDs")
      next
    }
    
    aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
    
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (is.null(aln) || length(aln) == 0) {
      message("Skipping ", gene, ": could not read alignment")
      next
    }
    
    names(aln) <- sub("\\|.*", "", names(aln))
    aln_mat <- as.matrix(aln)
    seq_ids <- names(aln)
    
    common_ids <- intersect(intersect(seq_ids, pheno_pcs_full$ID), clean_ids_gene)
    if (length(common_ids) < 8500L) {
      message("Skipping ", gene, ": insufficient overlap (n=", length(common_ids), ")")
      next
    }
    
    aln_mat <- aln_mat[match(common_ids, seq_ids), , drop = FALSE]
    pheno_pcs_sub <- pheno_pcs_full[ID %in% common_ids]
    setkey(pheno_pcs_sub, ID)
    pheno_pcs_sub <- pheno_pcs_sub[match(common_ids, ID)]
    
    y <- pheno_pcs_sub$pheno
    
    # Prepare PC matrix (or NULL if n_pcs == 0)
    if (n_pcs > 0) {
      X_pcs <- scale(as.matrix(pheno_pcs_sub[, ..pc_names_use]))
    } else {
      X_pcs <- NULL
    }
    
    # -------------------------------------------------------------------
    # LOOP OVER POSITIONS
    # -------------------------------------------------------------------
    
    positions <- seq_len(ncol(aln_mat))
    results_list <- vector("list", length(positions))
    
    for (pos in positions) {
      
      residues <- aln_mat[, pos]
      residue_table <- table(residues)
      
      if (length(residue_table) < 2L || all(names(residue_table) == "-")) next
      
      non_gap <- residues != "-"
      if (sum(non_gap) < 8500L) next
      
      res_factor <- factor(residues)
      X_aa <- model.matrix(~ res_factor - 1)
      
      var0 <- which(apply(X_aa, 2, var) == 0)
      if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
      if (ncol(X_aa) == 0) next
      
      # -----------------------------------------------------------------
      # FIT MODELS
      # -----------------------------------------------------------------
      
      fit_null <- tryCatch(lm(y ~ 1), error = function(e) NULL)
      if (is.null(fit_null)) next
      
      fit_aa <- tryCatch(lm(y ~ X_aa), error = function(e) NULL)
      if (is.null(fit_aa)) next
      
      if (n_pcs > 0) {
        fit_pcs <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
        fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
        if (is.null(fit_pcs) || is.null(fit_full)) next
        
        r2_pcs <- summary(fit_pcs)$r.squared
        r2_full <- summary(fit_full)$r.squared
        p_pcs_only <- tryCatch(anova(fit_null, fit_pcs)[2, "Pr(>F)"], error = function(e) NA_real_)
        p_aa_with_pcs <- tryCatch(anova(fit_pcs, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
        
        coef_summary <- summary(fit_full)$coefficients
      } else {
        r2_pcs <- NA_real_
        r2_full <- summary(fit_aa)$r.squared
        p_pcs_only <- NA_real_
        p_aa_with_pcs <- NA_real_
        
        coef_summary <- summary(fit_aa)$coefficients
      }
      
      r2_aa <- summary(fit_aa)$r.squared
      p_aa_only <- tryCatch(anova(fit_null, fit_aa)[2, "Pr(>F)"], error = function(e) NA_real_)
      
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
        n_pcs = n_pcs,
        R2_aa = r2_aa,
        R2_pcs = r2_pcs,
        R2_full = r2_full,
        P_aa_only = p_aa_only,
        P_pcs_only = p_pcs_only,
        P_aa_with_pcs = p_aa_with_pcs,
        residue_counts = as.list(residue_table),
        effects = effects_dt
      )
    }
    
    results_list <- Filter(Negate(is.null), results_list)
    
    if (length(results_list) > 0) {
      saveRDS(results_list, out_path)
      message("  ", gene, ": ", length(results_list), " positions")
    }
  }
}

message("\n========== EXPERIMENT COMPLETE ==========")
message("Results saved to: ", base_out_dir)