# gwas_onehot_baseline.R
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(tibble)

# LOAD COMMON DATA -----------------------------------------------------------
data <- read_parquet("data/processed_data.parquet")
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")
pcs_IDS <- aln_index$index
scores <- as.data.frame(ev_pcs$x)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
scores <- cbind(ID = pcs_IDS, scores)
n_pcs <- 1000
pc_names <- paste0("PC", seq_len(n_pcs))
scores <- scores %>% select(ID, all_of(pc_names))

# Pre-join phenotype and PCs once
pheno_pcs <- data %>% 
  select(ID, pheno = !!pheno_col) %>%
  inner_join(scores, by = "ID")
stopifnot(nrow(pheno_pcs) > 0)

# find files
emb_files <- list.files("data/embeddings/", pattern = "_residue_embeddings\\.parquet$", full.names = TRUE)
aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)
get_gene <- function(path) sub("_residue_embeddings\\.parquet|_AA_aligned\\.fasta", "", basename(path))
emb_genes <- get_gene(emb_files)
aln_genes <- get_gene(aln_files)
genes_to_process <- intersect(emb_genes, aln_genes)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

all_results_list <- list()

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  emb_file <- file.path("data/embeddings/", paste0(gene, "_residue_embeddings.parquet"))
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # Only read ID column from embeddings to check which IDs exist
  df_emb <- tryCatch(
    read_parquet(emb_file, col_select = "ID"), 
    error = function(e) NULL
  )
  if (is.null(df_emb)) {
    message("Skipping ", gene, ": could not read embeddings file")
    next
  }
  
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # Strip IDs
  aln_names <- names(aln)
  names(aln) <- sub("\\|.*", "", aln_names)
  emb_ids_base <- sub("\\|.*", "", df_emb$ID)
  
  # Find common IDs once
  common_ids <- intersect(emb_ids_base, pheno_pcs$ID)
  if (length(common_ids) == 0) {
    message("Skipping ", gene, ": no common IDs")
    next
  }
  
  # Subset alignment to common IDs
  aln_sub <- aln[names(aln) %in% common_ids]
  aln_mat <- as.matrix(aln_sub)
  
  # Pre-scale PCs once for this gene
  pheno_pcs_sub <- pheno_pcs %>% filter(ID %in% common_ids)
  X_pcs_scaled <- scale(as.matrix(pheno_pcs_sub[, pc_names]))
  
  results_list <- vector("list", ncol(aln_mat))
  
  for (pos in seq_len(ncol(aln_mat))) {
    residues <- aln_mat[, pos]
    
    # Skip gaps and check sample size
    non_gap <- residues != "-"
    if (sum(non_gap) < 8500) next
    
    # Check variation
    residue_table <- table(residues[non_gap])
    if (length(residue_table) < 2) next
    
    # Subset data
    y <- pheno_pcs_sub$pheno[non_gap]
    X_pcs <- X_pcs_scaled[non_gap, , drop = FALSE]
    res_factor <- factor(residues[non_gap])
    
    # One-hot encode
    X_aa <- model.matrix(~ res_factor - 1)
    zero_var_cols <- which(apply(X_aa, 2, function(x) var(as.numeric(x))) == 0)
    if (length(zero_var_cols) > 0) X_aa <- X_aa[, -zero_var_cols, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # Fit models
    df_reduced <- as.data.frame(cbind(y = y, X_pcs))
    fit_reduced <- tryCatch(lm(y ~ ., data = df_reduced), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    df_full <- as.data.frame(cbind(y = y, X_aa, X_pcs))
    fit_full <- tryCatch(lm(y ~ ., data = df_full), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_full <- summary(fit_full)$r.squared
    p_value <- tryCatch(anova(fit_reduced, fit_full)$`Pr(>F)`[2], error = function(e) NA_real_)
    loglik_ratio <- as.numeric(logLik(fit_full)[1] - logLik(fit_reduced)[1])
    
    results_list[[pos]] <- tibble(
      Gene = gene,
      Aligned_Position = pos,
      N = sum(non_gap),
      R2_full = r2_full,
      R2_reduced = r2_reduced,
      P_value = p_value,
      LogLik_Ratio = loglik_ratio
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) > 0) {
    all_results_list[[gene]] <- results_list
    saveRDS(results_list, paste0("results/tmp_results_onehot_", gene, ".rds"))
    message("Completed ", gene, ": ", length(results_list), " positions analyzed")
  } else {
    message("No positions retained for gene ", gene)
  }
}

saveRDS(all_results_list, "results/results_list_onehot_all_genes.rds")
message("\nSaved one-hot baseline results for ", length(all_results_list), " genes")
gene <- "atpA"
res_list <- readRDS(paste0("results/tmp_results_onehot_", gene, ".rds"))

# Combine list of per-site summaries into one data frame
res_df <- bind_rows(res_list) %>%
  mutate(log10p = -log10(P_value))

a <- res_df %>% 
  select(Aligned_Position, P_value, LogLik_Ratio, log10p) %>%
  distinct()

summary(a)
plot(a$log10p)
bf <- 0.001 / max(a$Aligned_Position)
abline(h=-log10(bf), col="coral")
#257 is the max