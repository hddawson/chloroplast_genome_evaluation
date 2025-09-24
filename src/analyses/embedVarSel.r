#!/usr/bin/env Rscript

# Performance analysis of embedding dimensions per gene
# Usage: Rscript embedding_performance_analysis.R

library(data.table)
library(ggplot2)
library(arrow)

# Set working directory 
setwd("/local/workdir/hdd29/chloroplast_genome_evaluation")

# Load data
cat("Loading data...\n")
embeds <- readRDS("data/tmp/embeds_with_mds.rds")
df <- read_parquet("data/processed_data.parquet")

# Data prep
setDT(embeds)
setDT(df)

clean_embeds <- embeds[ManualOutlier == FALSE]
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
genes <- unique(clean_embeds$Gene)
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

# Validation functions
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm=TRUE))

# Basic assertions
stopifnot(nrow(clean_embeds) > 0)
stopifnot(length(embed_cols) > 0)
stopifnot(pheno_col %in% colnames(df))

# Results storage
results <- data.table(n_top = integer(),
                      split = character(),
                      spearman = numeric(),
                      rmse = numeric())

set.seed(123)

cat("Running analysis for n_top 1 to 100...\n")

for (n_top in 1:100) {
  if (n_top %% 10 == 0) cat("n_top =", n_top, "\n")
  cat(n_top)  
  # Gene-level feature selection
  sel_list <- vector("list", length(genes))
  names(sel_list) <- genes
  
  for (gene in genes) {
    gdt <- clean_embeds[Gene == gene, c("ID", embed_cols), with = FALSE]
    gdt <- merge(gdt, df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = FALSE, all.y = FALSE)
    
    if (nrow(gdt) < 5) next
    
    # Correlation-based feature selection
    cors <- sapply(embed_cols, function(col) cor(gdt[[col]], gdt$pheno, use = "complete.obs"))
    available <- names(cors)[!is.na(cors)]
    k <- min(n_top, length(available))
    top_dims <- names(sort(abs(cors[available]), decreasing = TRUE))[1:k]
    
    # Select and rename features
    sel <- gdt[, c("ID", top_dims), with = FALSE]
    setnames(sel, old = top_dims, new = paste0(gene, "__", top_dims))
    sel_list[[gene]] <- sel
  }
  
  # Remove empty genes
  sel_list <- sel_list[!sapply(sel_list, is.null)]
  if (length(sel_list) == 0) next
  
  # Combine all gene features
  combined <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)
  combined <- merge(combined, df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE, all.y = FALSE)
  
  # Median imputation for missing values
  pred_cols <- setdiff(names(combined), c("ID", "pheno"))
  for (col in pred_cols) {
    med <- median(combined[[col]], na.rm = TRUE)
    combined[is.na(get(col)), (col) := med]
  }
  
  # Random 10% holdout
  idx_test <- sample(seq_len(nrow(combined)), size = ceiling(0.1 * nrow(combined)))
  train1 <- combined[-idx_test]
  test1  <- combined[idx_test]
  
  fit1 <- lm(pheno ~ ., data = train1[, c("pheno", pred_cols), with = FALSE])
  pred1 <- predict(fit1, newdata = test1)
  
  results <- rbind(results, data.table(
    n_top = n_top,
    split = "random10",
    spearman = cor(pred1, test1$pheno, method = "spearman", use="complete.obs"),
    rmse = rmse(pred1, test1$pheno)
  ))
  
  # Poaceae holdout
  poaceae_ids <- df[grepl("Poaceae", Taxonomy), ID]
  train2 <- combined[!ID %in% poaceae_ids]
  test2  <- combined[ID %in% poaceae_ids]
  
  if (nrow(test2) > 0 && nrow(train2) > 10) {
    fit2 <- lm(pheno ~ ., data = train2[, c("pheno", pred_cols), with = FALSE])
    pred2 <- predict(fit2, newdata = test2)
    
    results <- rbind(results, data.table(
      n_top = n_top,
      split = "poaceae",
      spearman = cor(pred2, test2$pheno, method = "spearman", use="complete.obs"),
      rmse = rmse(pred2, test2$pheno)
    ))
  }
}

# Save results
write.csv(results, "results/embedding_performance.csv", row.names=FALSE)

# Generate plot
p <- ggplot(results, aes(x = n_top, y = spearman, color = split)) +
  geom_line() +
  geom_point() +
  labs(title = "Effect of number of embedding dims per gene",
       x = "Number of top embedding dims per gene",
       y = "Spearman correlation (Predicted vs Observed)") +
  theme_minimal()

ggsave("plots/embedding_performance.png", p, width=8, height=6)

# Summary stats
cat("\nSummary of results:\n")
print(results[, .(max_spearman = max(spearman, na.rm=TRUE)), by=split])

best_random <- results[split=="random10"][which.max(spearman)]
best_poaceae <- results[split=="poaceae"][which.max(spearman)]

cat("\nBest random10 performance: n_top =", best_random$n_top, 
    ", spearman =", round(best_random$spearman, 3), "\n")
cat("Best poaceae performance: n_top =", best_poaceae$n_top, 
    ", spearman =", round(best_poaceae$spearman, 3), "\n")

cat("\nResults saved to results/embedding_performance.csv\n")
cat("Plot saved to plots/embedding_performance.png\n")
