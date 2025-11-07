library(arrow)
library(dplyr)
library(Biostrings)

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

# Find files
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
  
  # Read embeddings
  df_emb <- tryCatch(read_parquet(emb_file), error = function(e) NULL)
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
  names(aln) <- sub("\\|.*", "", names(aln))
  df_emb <- df_emb %>% mutate(ID = sub("\\|.*", "", ID))
  
  # Build residue index lookup
  aln_mat <- as.matrix(aln)
  residue_index <- lapply(seq_len(nrow(aln_mat)), function(i) {
    seq_id <- names(aln)[i]
    residues <- aln_mat[i, ]
    non_gap <- which(residues != "-")
    data.frame(
      ID = seq_id,
      Residue_Index = seq_along(non_gap),
      Aligned_Position = non_gap,
      stringsAsFactors = FALSE
    )
  })
  residue_lookup <- do.call(rbind, residue_index)
  
  # Join embeddings with alignment positions
  df_joined <- df_emb %>%
    inner_join(residue_lookup, by = c("ID", "Residue_Index")) %>%
    inner_join(pheno_pcs, by = "ID")
  
  if (nrow(df_joined) == 0) {
    message("Skipping ", gene, ": no joined data")
    next
  }
  
  # Get embedding columns
  emb_cols <- grep("^embedding_", colnames(df_joined), value = TRUE)
  stopifnot(length(emb_cols) > 0)
  
  # Pre-scale PCs for this gene
  common_ids <- unique(df_joined$ID)
  pheno_pcs_sub <- pheno_pcs %>% filter(ID %in% common_ids)
  X_pcs_scaled <- scale(as.matrix(pheno_pcs_sub[, pc_names]))
  rownames(X_pcs_scaled) <- pheno_pcs_sub$ID
  
  # Get unique aligned positions
  positions <- sort(unique(df_joined$Aligned_Position))
  results_list <- vector("list", length(positions))
  
  for (j in seq_along(positions)) {
    pos <- positions[j]
    sub <- df_joined %>% filter(Aligned_Position == pos)
    
    # Check sample size
    if (nrow(sub) < 8500) next
    
    # Check residue variation
    residues <- aln_mat[match(sub$ID, names(aln)), pos]
    residue_table <- table(residues)
    if (length(residue_table) < 2) next
    
    # Prepare data
    y <- sub$pheno
    X_pcs <- X_pcs_scaled[sub$ID, , drop = FALSE]
    X_emb <- scale(as.matrix(sub[, emb_cols]))
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    zero_var_cols <- which(apply(X_aa, 2, function(x) var(as.numeric(x))) == 0)
    if (length(zero_var_cols) > 0) X_aa <- X_aa[, -zero_var_cols, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # Fit three nested models:
    # Reduced: PCs only
    # Partial: PCs + Residues
    # Full: PCs + Residues + Embeddings
    
    df_reduced <- as.data.frame(cbind(y = y, X_pcs))
    fit_reduced <- tryCatch(lm(y ~ ., data = df_reduced), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    df_partial <- as.data.frame(cbind(y = y, X_aa, X_pcs))
    fit_partial <- tryCatch(lm(y ~ ., data = df_partial), error = function(e) NULL)
    if (is.null(fit_partial)) next
    
    df_full <- as.data.frame(cbind(y = y, X_aa, X_emb, X_pcs))
    fit_full <- tryCatch(lm(y ~ ., data = df_full), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    # Extract metrics
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_partial <- summary(fit_partial)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    # Test reduced vs partial (do residues help?)
    p_res <- tryCatch(anova(fit_reduced, fit_partial)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    # Test partial vs full (do embeddings help beyond residues?)
    p_emb <- tryCatch(anova(fit_partial, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    results_list[[j]] <- data.frame(
      Gene = gene,
      Aligned_Position = pos,
      N = nrow(sub),
      R2_reduced = r2_reduced,
      R2_partial = r2_partial,
      R2_full = r2_full,
      P_res = p_res,
      P_emb = p_emb,
      stringsAsFactors = FALSE
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) > 0) {
    all_results_list[[gene]] <- results_list
    saveRDS(results_list, paste0("results/tmp_results_nested_", gene, ".rds"))
    message("Completed ", gene, ": ", length(results_list), " positions analyzed")
  } else {
    message("No positions retained for gene ", gene)
  }
}

saveRDS(all_results_list, "results/results_list_nested_all_genes.rds")
message("\nSaved nested model results for ", length(all_results_list), " genes")