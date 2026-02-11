library(ape)
library(arrow)
library(data.table)
library(phylolm)

# ---- 1. LOAD DATA ----
tree <- read.tree("results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree")
data <- as.data.table(read_parquet("data/processed_data.parquet"))

# Drop gymnosperms
pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree$tip.label)
#tree <- drop.tip(tree, pinales_in_tree)
message("Tree: ", Ntip(tree), " tips after dropping gymnosperms")

# ---- 2. PREPARE PHENOTYPE ----
pheno <- setNames(data$pheno_wc2.1_2.5m_bio_8_p50, data$ID)  # Using MAT as phenotype - adjust if needed
pheno <- pheno[tree$tip.label]
pheno <- pheno[!is.na(pheno)]

# Keep only tips with phenotype data
tree <- keep.tip(tree, names(pheno))
message("After filtering: ", Ntip(tree), " tips with phenotype data")
is.rooted(tree)
stopifnot(length(pheno) == Ntip(tree))
stopifnot(all(names(pheno) %in% tree$tip.label))

# ---- 3. EXTRACT GENES FROM data$Genes ----
all_genes <- unique(unlist(lapply(data$Genes, function(x) {
  # Remove brackets and quotes, split by comma
  genes <- gsub("\\[|\\]|'", "", x)
  genes <- strsplit(genes, ", ")[[1]]
  genes[genes != ""]
})))

message("Found ", length(all_genes), " unique genes")

# Create presence/absence matrix
gene_wide <- data.table(ID = data$ID)
for (gene in all_genes) {
  gene_wide[[gene]] <- sapply(data$Genes, function(x) {
    as.integer(grepl(paste0("'", gene, "'"), x, fixed = TRUE))
  })
}

# Match to tree tips
gene_wide <- gene_wide[ID %in% tree$tip.label]
setkey(gene_wide, ID)
gene_wide <- gene_wide[match(tree$tip.label, ID)]  # Order by tree

stopifnot(all(gene_wide$ID == tree$tip.label))
stopifnot(all(gene_wide$ID == names(pheno)))

# ---- 4. RUN PHYLOLM FOR EACH GENE ----
gene_cols <- all_genes
message("\nTesting ", length(gene_cols), " genes with phylolm...")

results <- rbindlist(lapply(gene_cols, function(gene) {
  test_data <- data.frame(
    pheno = pheno,
    predictor = gene_wide[[gene]]
  )
  
  # Skip if no variation
  if (length(unique(test_data$predictor)) < 2) {
    return(data.table(gene = gene, coefficient = NA, se = NA, tvalue = NA, 
                      pvalue = NA, lambda = NA, note = "no_variation"))
  }
  
  tryCatch({
    fit <- phylolm(pheno ~ predictor, data = test_data, phy = tree, model = "lambda")
    print(summary(fit))
    coef_summary <- summary(fit)$coefficients
    
    data.table(
      gene = gene,
      coefficient = coef_summary["predictor", "Estimate"],
      se = coef_summary["predictor", "StdErr"],
      tvalue = coef_summary["predictor", "t.value"],
      pvalue = coef_summary["predictor", "p.value"],
      lambda = fit$optpar,
      note = NA_character_
    )
  }, error = function(e) {
    data.table(gene = gene, coefficient = NA, se = NA, tvalue = NA, 
               pvalue = NA, lambda = NA, note = as.character(e$message))
  })
}))

# ---- 5. MULTIPLE TEST CORRECTION ----
results[, padj := p.adjust(pvalue, method = "fdr")]
setorder(results, pvalue)

print(results[1:min(20, nrow(results))])

# Build full model with all genes
gene_matrix <- as.matrix(gene_wide[, ..gene_cols])
rownames(gene_matrix) <- gene_wide$ID

test_data <- data.frame(
  pheno = pheno,
  gene_matrix
)

# Remove genes with no variation
gene_vars <- sapply(gene_wide[, ..gene_cols], var)
variable_genes <- gene_cols[gene_vars > 0]
message("Kept ", length(variable_genes), " / ", length(gene_cols), " genes with variation")

# Remove highly correlated genes
gene_matrix <- as.matrix(gene_wide[, ..variable_genes])
cor_matrix <- cor(gene_matrix)
library(caret)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.95)  # requires caret package

if (length(high_cor) > 0) {
  variable_genes <- variable_genes[-high_cor]
  message("Removed ", length(high_cor), " genes with |r| > 0.95")
}

# Now fit with reduced gene set
test_data <- data.frame(
  pheno = pheno,
  gene_wide[, ..variable_genes]
)

fit_full <- phylolm(pheno ~ ., data = test_data, phy = tree, model = "lambda")
summary(fit_full)

# Or use model selection
library(MuMIn)
fit_full <- phylolm(pheno ~ ., data = test_data, phy = tree, model = "lambda", na.action = "na.fail")
best_models <- dredge(fit_full, rank = "AIC")

# ---- 6. SAVE RESULTS ----
dir.create("results", showWarnings = FALSE)
saveRDS(results, "results/phylolm_gene_associations.rds")
fwrite(results, "results/phylolm_gene_associations.csv")

message("\n=== TOP ASSOCIATIONS ===")
print(results[1:min(20, nrow(results))])

message("\nSignificant hits (padj < 0.05): ", sum(results$padj < 0.05, na.rm = TRUE))
message("Genes with errors: ", sum(!is.na(results$note)))
message("\nResults saved to results/phylolm_gene_associations.csv")

# ---- Generate summary plot ----
library(ggplot2)

# Extract coefficients and p-values
coef_df <- data.frame(
  gene = rownames(summary(fit_full)$coefficients)[-1],  # exclude intercept
  coefficient = summary(fit_full)$coefficients[-1, "Estimate"],
  pvalue = summary(fit_full)$coefficients[-1, "p.value"]
)
coef_df$neglog10p <- -log10(coef_df$pvalue)

# Volcano plot
p <- ggplot(coef_df, aes(x = coefficient, y = neglog10p)) +
  geom_point(aes(color = pvalue < 0.05), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(data = subset(coef_df, pvalue < 0.05), 
            aes(label = gene), size = 3, hjust = -0.1, vjust = 0.5) +
  scale_color_manual(values = c("grey60", "steelblue"), guide = "none") +
  labs(title = "Chloroplast Gene Presence/Absence vs Phenotype",
       subtitle = paste0("Phylogenetic GLM (lambda=", round(fit_full$optpar, 2), 
                         ") | R² = ", round(fit_full$r.squared, 3)),
       x = "Coefficient Estimate",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
p
ggsave("results/gene_phenotype_volcano.png", p, width = 8, height = 6, dpi = 150)

# ---- Generate table for Slack ----
sig_genes <- coef_df[coef_df$pvalue < 0.05, ]
sig_genes <- sig_genes[order(sig_genes$pvalue), ]

cat("\n=== COPY THIS FOR SLACK ===\n\n")
cat("```\n")
cat(sprintf("%-15s %10s %10s\n", "Gene", "Coef", "P-value"))
cat(sprintf("%-15s %10s %10s\n", "---------------", "----------", "----------"))
for (i in 1:nrow(sig_genes)) {
  cat(sprintf("%-15s %10.3f %10.4f\n", 
              sig_genes$gene[i], 
              sig_genes$coefficient[i], 
              sig_genes$pvalue[i]))
}
cat("```\n")