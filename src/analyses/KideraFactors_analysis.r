# ==== KIDERA FACTORS BY GENE: COMPREHENSIVE ANALYSIS ====
library(ape)
library(data.table)
library(Peptides)
library(ggplot2)

# ==== CONFIGURATION ====
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
partition_map_file <- "raxml_input/partitionMap.rds"
gwas_dir <- "results/residue_models_triple/"

# ==== LOAD DATA ====
cat("Loading data...\n")
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)

data <- as.data.table(arrow::read_parquet(data_file))
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID

partition_map <- readRDS(partition_map_file)
setDT(partition_map)

# ==== LOAD AND CLASSIFY GWAS EFFECTS ====
model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_effects <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    data.table(Gene = m$Gene, GenePos = m$Aligned_Position,
               P_aa_with_pcs = m$P_aa_with_pcs, P_aa_only = m$P_aa_only)
  }), fill = TRUE)
}), fill = TRUE)

gwas_effects <- unique(gwas_effects)
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

thresh_control <- quantile(gwas_effects$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(gwas_effects$P_aa_only, 0.20, na.rm = TRUE)
gwas_effects[, sig_class := fcase(
  P_aa_with_pcs < thresh_control & P_aa_only < thresh_nocontrol, "sig_both",
  P_aa_with_pcs < thresh_control & P_aa_only >= thresh_nocontrol, "sig_control",
  P_aa_with_pcs >= thresh_control & P_aa_only < thresh_nocontrol, "sig_nocontrol",
  default = "not_sig"
)]

# ==== FUNCTION: CALCULATE SINGLE KF BY GENE ====
calc_kf_by_gene <- function(aln_mat, positions, partition_map, kf_index) {
  stopifnot(all(positions <= ncol(aln_mat)))
  pos_genes <- partition_map[GlobalPos %in% positions, .(Gene, GlobalPos)]
  
  gene_kf <- rbindlist(lapply(unique(pos_genes$Gene), function(g) {
    gene_pos <- pos_genes[Gene == g, GlobalPos]
    if (length(gene_pos) < 1) return(NULL)
    
    sub_mat <- aln_mat[, gene_pos, drop = FALSE]
    seqs <- apply(sub_mat, 1, function(row) {
      row[!row %in% c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- ""
      paste(row, collapse = "")
    })
    
    kf_vals <- sapply(seqs, function(s) {
      if (nchar(s) < 1) return(NA)
      kideraFactors(s)[[1]][kf_index]
    })
    
    data.table(taxon = names(kf_vals), Gene = g, KF_value = kf_vals, n_sites = length(gene_pos))
  }))
  gene_kf
}

# ==== KF LABELS ====
kf_info <- data.table(
  kf_index = 1:10,
  kf_name = paste0("KF", 1:10),
  kf_desc = c("Helix/bend preference", "Side-chain size", "Extended structure preference",
              "Hydrophobicity", "Double-bend preference", "Partial specific volume",
              "Flat extended preference", "Occurrence in alpha region", "pK-C",
              "Surrounding hydrophobicity")
)

# ==== COMPUTE ALL KF BY GENE FOR sig_both AND not_sig ====
cat("Computing Kidera factors by gene...\n")

all_kf_data <- rbindlist(lapply(1:10, function(kf_idx) {
  cat("  KF", kf_idx, "\n")
  
  kf_by_gene <- rbindlist(lapply(c("not_sig", "sig_both"), function(sc) {
    pos <- gwas_effects[sig_class == sc, GlobalPos]
    pos <- pos[pos <= ncol(aln_mat)]
    if (length(pos) < 10) return(NULL)
    
    gene_kf <- calc_kf_by_gene(aln_mat, pos, partition_map, kf_idx)
    gene_kf[, sig_class := sc]
    gene_kf
  }))
  
  kf_by_gene$taxon <- gsub(paste0("\\.KF", kf_idx, "$"), "", kf_by_gene$taxon)
  kf_by_gene[, KF := paste0("KF", kf_idx)]
  kf_by_gene
}))

all_kf_data[, phenotype := pheno[taxon]]
all_kf_data <- all_kf_data[!is.na(KF_value) & !is.na(phenotype)]

cat("Total observations:", nrow(all_kf_data), "\n")

# ==== STATISTICAL TESTS ====
cat("\nRunning correlation tests...\n")

kf_cors <- all_kf_data[, .(
  r = cor(KF_value, phenotype, use = "complete.obs"),
  p = tryCatch(cor.test(KF_value, phenotype)$p.value, error = function(e) NA),
  n = .N
), by = .(Gene, sig_class, KF)]

# ANOVA: KF differs by sig_class within genes?
anova_results <- all_kf_data[, {
  if (length(unique(sig_class)) > 1 && .N > 30) {
    fit <- aov(KF_value ~ sig_class, data = .SD)
    s <- summary(fit)[[1]]
    .(F_val = s$`F value`[1], anova_p = s$`Pr(>F)`[1])
  } else {
    .(F_val = NA_real_, anova_p = NA_real_)
  }
}, by = .(Gene, KF)]

# ==== SUMMARY STATISTICS BY KF ====
cat("\n=== Summary by Kidera Factor ===\n")

kf_summary <- kf_cors[, .(
  n_sig_cors = sum(p < 0.05, na.rm = TRUE),
  n_genes = .N,
  pct_sig = 100 * sum(p < 0.05, na.rm = TRUE) / .N,
  mean_abs_r = mean(abs(r), na.rm = TRUE),
  max_abs_r = max(abs(r), na.rm = TRUE),
  median_p = median(p, na.rm = TRUE)
), by = KF]

anova_summary <- anova_results[, .(
  n_sig_anova = sum(anova_p < 0.05, na.rm = TRUE),
  mean_F = mean(F_val, na.rm = TRUE)
), by = KF]

kf_summary <- merge(kf_summary, anova_summary, by = "KF")
kf_summary <- merge(kf_summary, kf_info[, .(kf_name, kf_desc)], by.x = "KF", by.y = "kf_name")
kf_summary <- kf_summary[order(-pct_sig)]

print(kf_summary)

# ==== VISUALIZATIONS ====

# 1. Faceted scatter plots (like _x.pdf)
cat("\nGenerating plots...\n")

for (kf_idx in 1:10) {
  kf_name <- paste0("KF", kf_idx)
  kf_desc <- kf_info[kf_index == kf_idx, kf_desc]
  
  plot_data <- all_kf_data[KF == kf_name]
  
  p <- ggplot(plot_data, aes(x = KF_value, y = phenotype, color = sig_class)) +
    geom_point(alpha = 0.3, size = 0.8) +
    facet_wrap(~ Gene, scales = "free_x") +
    labs(x = paste0(kf_name, ": ", kf_desc), y = "Phenotype (bio_8_p50)",
         title = paste0("Gene-specific ", kf_name, " vs Phenotype"),
         color = "Sig Class") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(sprintf("results/kf%d_by_gene.pdf", kf_idx), p, width = 12, height = 10)
}

# 2. Summary heatmap: correlation strength by KF and Gene
cor_mat <- dcast(kf_cors[sig_class == "sig_both"], Gene ~ KF, value.var = "r", fill = 0)
genes <- cor_mat$Gene
cor_mat[, Gene := NULL]
cor_mat <- as.matrix(cor_mat)
rownames(cor_mat) <- genes

pdf("results/kf_gene_correlation_heatmap.pdf", width = 10, height = 14)
pheatmap::pheatmap(cor_mat,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   breaks = seq(-0.3, 0.3, length.out = 101),
                   main = "KF-Phenotype Correlations by Gene (sig_both sites)")
dev.off()

# 3. Summary bar plot: biological signal strength
p_summary <- ggplot(kf_summary, aes(x = reorder(KF, -pct_sig), y = pct_sig)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", pct_sig)), vjust = -0.3, size = 3) +
  labs(x = "Kidera Factor", y = "% Genes with Significant Correlation (p<0.05)",
       title = "Biological Signal Strength by Kidera Factor") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/kf_signal_strength.pdf", p_summary, width = 8, height = 6)

# ==== DETAILED REPORT ====
cat("\n========================================\n")
cat("KIDERA FACTOR BIOLOGICAL SIGNAL REPORT\n")
cat("========================================\n\n")

for (i in 1:nrow(kf_summary)) {
  row <- kf_summary[i]
  cat(sprintf("%s: %s\n", row$KF, row$kf_desc))
  cat(sprintf("  Significant correlations: %d/%d genes (%.1f%%)\n", 
              row$n_sig_cors, row$n_genes, row$pct_sig))
  cat(sprintf("  Mean |r|: %.3f, Max |r|: %.3f\n", row$mean_abs_r, row$max_abs_r))
  cat(sprintf("  Genes with sig. ANOVA (sig_class effect): %d\n", row$n_sig_anova))
  
  # Top genes for this KF
  top <- kf_cors[KF == row$KF & p < 0.05][order(p)][1:5]
  if (nrow(top) > 0) {
    cat("  Top genes:\n")
    for (j in 1:nrow(top)) {
      cat(sprintf("    %s (%s): r=%.3f, p=%.2e\n", 
                  top[j, Gene], top[j, sig_class], top[j, r], top[j, p]))
    }
  }
  cat("\n")
}

# ==== SAVE RESULTS ====
saveRDS(list(
  kf_cors = kf_cors,
  anova_results = anova_results,
  kf_summary = kf_summary,
  all_kf_data = all_kf_data
), "data/tmp/all_kf_gene_analysis.rds")

cat("\nSaved: data/tmp/all_kf_gene_analysis.rds\n")
cat("Saved: results/kf[1-10]_by_gene.pdf\n")
cat("Saved: results/kf_gene_correlation_heatmap.pdf\n")
cat("Saved: results/kf_signal_strength.pdf\n")