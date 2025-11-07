library(arrow)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(Biostrings)
library(tibble)

# Load shared data once
data <- read_parquet("data/processed_data.parquet")
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln <- read_parquet("data/tmp/majMinor_aln.pq")
pcs_IDS <- aln$index
scores <- as.data.frame(ev_pcs$x)
scores <- cbind(ID = pcs_IDS, scores)
n_pcs <- 1000
pc_names <- paste0("PC", seq_len(n_pcs))
colnames(scores)[-1] <- paste0("PC", seq_len(ncol(scores)-1))
scores <- scores %>% select(ID, all_of(pc_names))

# Find matching embedding and alignment files
emb_files <- list.files("data/embeddings/", pattern = "_residue_embeddings\\.parquet$", full.names = TRUE)
aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

# Extract gene names
get_gene <- function(path) sub("_residue_embeddings\\.parquet|_AA_aligned\\.fasta", "", basename(path))
emb_genes <- get_gene(emb_files)
aln_genes <- get_gene(aln_files)

# Only process genes with both files
genes_to_process <- intersect(emb_genes, aln_genes)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

all_results_list <- list()

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  
  # Load gene-specific files
  emb_file <- file.path("data/embeddings/", paste0(gene, "_residue_embeddings.parquet"))
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  df <- tryCatch(read_parquet(emb_file), error = function(e) NULL)
  if (is.null(df)) {
    message("Skipping ", gene, ": could not read embeddings")
    next
  }
  
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # Create alignment lookup
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  colnames(aln_mat) <- seq_len(ncol(aln_mat))
  
  aln_lookup <- as.data.frame(aln_mat) %>%
    tibble::rownames_to_column("ID") %>%
    tidyr::pivot_longer(-ID, names_to = "Aligned_Position", values_to = "Residue") %>%
    mutate(Aligned_Position = as.integer(gsub("V", "", Aligned_Position))) %>%
    group_by(ID) %>%
    mutate(Residue_Index = ifelse(Residue == "-", NA_integer_, cumsum(Residue != "-"))) %>%
    ungroup() %>%
    filter(!is.na(Residue_Index)) %>%
    select(ID, Residue_Index, Aligned_Position)
  
  # Join data
  common_ids <- intersect(df$ID, pcs_IDS)
  if (length(common_ids) == 0) {
    message("Skipping ", gene, ": no common IDs")
    next
  }
  
  df_joined <- df %>%
    mutate(ID = sub("\\|.*", "", ID)) %>%
    inner_join(aln_lookup, by = c("ID", "Residue_Index")) %>%
    inner_join(data %>% select(ID, pheno = !!pheno_col), by = "ID") %>%
    inner_join(scores, by = "ID") %>%
    filter(ID %in% common_ids)
  
  if (nrow(df_joined) == 0) {
    message("Skipping ", gene, ": no joined data")
    next
  }
  
  df_joined <- df_joined %>% 
    mutate(GroupID = paste(Gene, Aligned_Position, sep = "_"))
  
  groups <- unique(df_joined$GroupID)
  results_list <- vector("list", length(groups))
  emb_cols <- grep("^embedding_", colnames(df_joined), value = TRUE)
  
  for (i in seq_along(groups)) {
    gid <- groups[i]
    sub <- df_joined %>% filter(GroupID == gid)
    
    if (nrow(sub) < 8500) next
    
    y <- sub$pheno
    X_emb <- scale(as.matrix(sub[, emb_cols]))
    X_pcs <- scale(as.matrix(sub[, pc_names]))
    
    df_reduced <- as.data.frame(cbind(y = y, X_pcs))
    fit_reduced <- tryCatch(lm(y ~ ., data = df_reduced), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    df_full <- as.data.frame(cbind(y = y, X_emb, X_pcs))
    fit_full <- tryCatch(lm(y ~ ., data = df_full), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_full <- summary(fit_full)$r.squared
    p_value <- anova(fit_reduced, fit_full)$`Pr(>F)`[2]
    loglik_ratio <- logLik(fit_full)[1] - logLik(fit_reduced)[1]
    
    coefs <- as.data.frame(coef(fit_full)) %>%
      tibble::rownames_to_column("Predictor")
    colnames(coefs)[2] <- "Estimate"
    
    coefs <- coefs %>%
      dplyr::filter(Predictor != "(Intercept)") %>%
      dplyr::mutate(
        Gene             = sub$Gene[1],
        Aligned_Position = sub$Aligned_Position[1],
        Residue          = sub$Residue[1],
        N                = nrow(sub),
        R2_full          = r2_full,
        r2_reduced       = r2_reduced,
        P_value          = p_value,
        LogLik_Ratio     = loglik_ratio
      )
    
    results_list[[i]] <- coefs
  }
  saveRDS(results_list, paste0("results/tmp_results_", gene, "_cve.rds"))
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) > 0) {
    all_results_list[[gene]] <- results_list
    message("Completed ", gene, ": ", length(results_list), " positions analyzed")
  }
}

saveRDS(all_results_list, "results/results_list_cve_all_genes.rds")
message("\nSaved results for ", length(all_results_list), " genes")

# Combine into single dataframe
results_df <- do.call(rbind, unlist(all_results_list, recursive = FALSE))
saveRDS(results_df, "results/results_df_cve_all_genes.rds")


###atpA 
# Load results for atpA
gene <- "atpA"
atpA_gwas <- readRDS(paste0("results/tmp_results_", gene, "_cve.rds"))
atpA_gwas <- Filter(Negate(is.null), atpA_gwas)
atpA_gwas_res <- do.call(rbind, atpA_gwas)


res_df <- atpA_gwas_res %>% 
  select(LogLik_Ratio, Aligned_Position, Residue, R2_full,
         r2_reduced, P_value) %>% 
  distinct() 

res_df$R2_diff <- res_df$R2_full - res_df$r2_reduced
plot(res_df$Aligned_Position, res_df$R2_diff)
lines(res_df$Aligned_Position, res_df$R2_diff)

#
gene <- "atpB"
gwas <- readRDS(paste0("results/tmp_results_", gene, "_cve.rds"))
gwas <- Filter(Negate(is.null), gwas)
gwas_res <- do.call(rbind, gwas)


res_df <- gwas_res %>% 
  select(LogLik_Ratio, Aligned_Position, Residue, R2_full,
         r2_reduced, P_value) %>% 
  distinct() 

res_df$R2_diff <- res_df$R2_full - res_df$r2_reduced
hist(res_df$R2_diff, main=gene)
plot(res_df$Aligned_Position, res_df$R2_diff,main=gene)
lines(res_df$Aligned_Position, res_df$R2_diff)

lines(res_df$Aligned_Position, res_df$LogLik_Ratio)

psbA <- do.call(rbind,readRDS("results/results_list_cve.rds"))
hist(psbA$LogLik_Ratio)


psbA_res_df <- psbA %>% 
  select(LogLik_Ratio, Aligned_Position, Residue, R2_full, r2_reduced, P_value) %>% 
  distinct() 
psbA_res_df$R2_diff <- psbA_res_df$R2_full - psbA_res_df$r2_reduced
plot(psbA_res_df$Aligned_Position,psbA_res_df$R2_diff)
lines(psbA_res_df$Aligned_Position,psbA_res_df$R2_diff)

hist(-log10(psbA_res_df$P_value))
hist(-log10(res_df$P_value))

hist(psbA_res_df$R2_diff)
hist(res_df$R2_diff)


library(dplyr)
library(ggplot2)

# Find all intermediate result files
result_files <- list.files("results/", pattern = "^tmp_results_.*_cve\\.rds$", full.names = TRUE)
stopifnot(length(result_files) > 0)

# Extract gene names
gene_names <- sub("^tmp_results_", "", sub("_cve\\.rds$", "", basename(result_files)))

message("Found ", length(result_files), " genes: ", paste(gene_names, collapse = ", "))

# Process each gene
all_gene_summaries <- list()

for (i in seq_along(result_files)) {
  gene <- gene_names[i]
  message("Processing ", gene)
  
  gwas <- readRDS(result_files[i])
  gwas <- Filter(Negate(is.null), gwas)
  
  if (length(gwas) == 0) {
    message("Skipping ", gene, ": no results")
    next
  }
  
  gwas_res <- do.call(rbind, gwas)
  
  res_df <- gwas_res %>% 
    select(LogLik_Ratio, Aligned_Position, Residue, R2_full, r2_reduced, P_value) %>% 
    distinct() %>%
    mutate(
      R2_diff = R2_full - r2_reduced,
      neg_log10_P = -log10(P_value),
      Gene = gene
    )
  
  stopifnot(nrow(res_df) > 0)
  
  all_gene_summaries[[gene]] <- res_df
  
  # Create plots for this gene
  pdf(paste0("results/EDA_", gene, ".pdf"), width = 12, height = 8)
  par(mfrow = c(2, 3))
  
  # R2_diff by position
  plot(res_df$Aligned_Position, res_df$R2_diff, 
       main = paste(gene, "- R² improvement by position"),
       xlab = "Aligned Position", ylab = "ΔR²")
  lines(res_df$Aligned_Position, res_df$R2_diff)
  
  # LogLik by position
  plot(res_df$Aligned_Position, res_df$LogLik_Ratio,
       main = paste(gene, "- Log-likelihood ratio by position"),
       xlab = "Aligned Position", ylab = "Log-Lik Ratio")
  lines(res_df$Aligned_Position, res_df$LogLik_Ratio)
  
  # -log10(P) by position
  plot(res_df$Aligned_Position, res_df$neg_log10_P,
       main = paste(gene, "- Significance by position"),
       xlab = "Aligned Position", ylab = "-log10(P)")
  lines(res_df$Aligned_Position, res_df$neg_log10_P)
  abline(h = -log10(0.05), col = "red", lty = 2)
  
  # Histograms
  hist(res_df$R2_diff, main = paste(gene, "- ΔR² distribution"), xlab = "ΔR²")
  hist(res_df$LogLik_Ratio, main = paste(gene, "- Log-Lik ratio distribution"), xlab = "Log-Lik Ratio")
  hist(res_df$neg_log10_P, main = paste(gene, "- -log10(P) distribution"), xlab = "-log10(P)")
  
  dev.off()
  message("Saved plots for ", gene)
}

# Combine all genes
all_res_df <- do.call(rbind, all_gene_summaries)
saveRDS(all_res_df, "results/all_genes_summary_stats.rds")

# Cross-gene comparison plots
pdf("results/EDA_all_genes_comparison.pdf", width = 14, height = 10)

# R2_diff comparison
boxplot(R2_diff ~ Gene, data = all_res_df, las = 2,
        main = "ΔR² distribution by gene", ylab = "ΔR²")

# LogLik comparison
boxplot(LogLik_Ratio ~ Gene, data = all_res_df, las = 2,
        main = "Log-likelihood ratio by gene", ylab = "Log-Lik Ratio")

# P-value comparison
boxplot(neg_log10_P ~ Gene, data = all_res_df, las = 2,
        main = "-log10(P) by gene", ylab = "-log10(P)")
abline(h = -log10(0.05), col = "red", lty = 2)

# Summary table
summary_table <- all_res_df %>%
  group_by(Gene) %>%
  summarise(
    N_positions = n(),
    Mean_R2_diff = mean(R2_diff),
    Max_R2_diff = max(R2_diff),
    N_sig = sum(P_value < 0.05),
    Prop_sig = mean(P_value < 0.05)
  ) %>%
  arrange(desc(Mean_R2_diff))

print(summary_table)

dev.off()

# Save summary table
saveRDS(summary_table, "results/gene_summary_table.rds")
write.csv(summary_table, "results/gene_summary_table.csv", row.names = FALSE)

message("\nCompleted EDA for ", length(all_gene_summaries), " genes")
message("Summary statistics saved to results/all_genes_summary_stats.rds")