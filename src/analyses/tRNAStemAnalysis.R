library(data.table)
library(arrow)

data <- read_parquet("data/processed_data.parquet")

trna_data <- fread("results/tRNA_stem_gc.tsv")

table(trna_data$gene)
length(unique(trna_data$gene))

summary(trna_data)

trna_data <- merge(trna_data, 
                   data[, c("ID", "pheno_wc2.1_2.5m_bio_8_p50")], 
                   by.x = "accession", 
                   by.y = "ID", 
                   all.x = TRUE)
# Quick check
sum(is.na(trna_data$pheno_wc2.1_2.5m_bio_8_p50))  # how many missing temp?

numeric_cols <- sapply(trna_data, is.numeric)

hist(trna_data$total_stem_gc)
hist(trna_data$total_loop_gc)


# Or correlation matrix heatmap
cor_mat <- cor(trna_data[, ..numeric_cols], use = "complete.obs")
heatmap(cor_mat, scale="row")

cor_mat

# Or just print correlations with temp
cor(trna_data[, ..numeric_cols], trna_data$pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs")
data_numeric_cols <- sapply(data, is.numeric)
data_cors <- cor(data[, data_numeric_cols], data$pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs")
hist(data_cors)


fit <- lm(pheno_wc2.1_2.5m_bio_8_p50~., data=trna_data[,..numeric_cols])
summary(fit)

setDT(trna_data)

# Mean GC across all genes per species
species_means <- trna_data[, .(
  total_stem_gc = mean(total_stem_gc, na.rm = TRUE),
  total_loop_gc = mean(total_loop_gc, na.rm = TRUE),
  acceptor_stem_gc = mean(acceptor_stem_gc, na.rm = TRUE),
  d_stem_gc = mean(d_stem_gc, na.rm = TRUE),
  anticodon_stem_gc = mean(anticodon_stem_gc, na.rm = TRUE),
  t_stem_gc = mean(t_stem_gc, na.rm = TRUE),
  anticodon_loop_gc = mean(anticodon_loop_gc, na.rm = TRUE),
  t_loop_gc = mean(t_loop_gc, na.rm = TRUE),
  temp = mean(pheno_wc2.1_2.5m_bio_8_p50, na.rm = TRUE),
  n_genes = .N
), by = accession]

# Check
nrow(species_means)
summary(species_means$n_genes)

mean_numeric_cols <- sapply(species_means, is.numeric)


# Correlations at species level
cor(species_means[, ..mean_numeric_cols], use = "complete.obs")

# Quick plot
plot(species_means$temp, species_means$total_stem_gc, 
     xlab = "Temperature", ylab = "Mean Stem GC",
     pch = 16, cex = 0.5)
abline(lm(total_stem_gc ~ temp, data = species_means), col = "red")

# Fit
summary(lm(temp ~ total_stem_gc + total_loop_gc, data = species_means))


# Per-gene correlations
gene_cors <- trna_data[, .(
  cor_stem = cor(total_stem_gc, pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs"),
  cor_loop = cor(total_loop_gc, pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs"),
  cor_acceptor = cor(acceptor_stem_gc, pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs"),
  cor_t_stem = cor(t_stem_gc, pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs"),
  cor_t_loop = cor(t_loop_gc, pheno_wc2.1_2.5m_bio_8_p50, use = "complete.obs"),
  n = .N
), by = gene]

# Sort by strongest stem correlation
gene_cors[order(-abs(cor_stem))]
summary(gene_cors)
# Or visualize
barplot(gene_cors$cor_stem, names.arg = gene_cors$gene, las = 2, 
        main = "Stem GC ~ Temp correlation by gene",
        ylab = "Pearson r")
abline(h = 0)

cor(data$pheno_wc2.1_2.5m_bio_8_p50, data$geno_genomeLength)
cor(data$pheno_wc2.1_2.5m_bio_8_p50, data$geno_genomicGC)

library(lme4)
fit_mixed <- lmer(pheno_wc2.1_2.5m_bio_8_p50 ~ total_stem_gc + total_loop_gc + (1|gene), 
                  data = trna_data)
summary(fit_mixed)

# Focus on trnL-UAG (strongest)
trnL_UAG <- trna_data[gene == "trnL-UAG"]

# Scatterplot
par(mfrow = c(1, 2))
plot(trnL_UAG$pheno_wc2.1_2.5m_bio_8_p50, trnL_UAG$total_stem_gc,
     xlab = "Temp", ylab = "Stem GC", main = "trnL-UAG", pch = 16, cex = 0.3)
abline(lm(total_stem_gc ~ pheno_wc2.1_2.5m_bio_8_p50, data = trnL_UAG), col = "red")

# Compare to trnC-GCA (opposite direction)
trnC <- trna_data[gene == "trnC-GCA"]
plot(trnC$pheno_wc2.1_2.5m_bio_8_p50, trnC$total_stem_gc,
     xlab = "Temp", ylab = "Stem GC", main = "trnC-GCA", pch = 16, cex = 0.3)
abline(lm(total_stem_gc ~ pheno_wc2.1_2.5m_bio_8_p50, data = trnC), col = "red")

# Formal test for top gene
summary(lm(total_stem_gc ~ pheno_wc2.1_2.5m_bio_8_p50, data = trnL_UAG))

# PCA on species means (just the GC columns)
gc_cols <- c("total_stem_gc", "total_loop_gc", "acceptor_stem_gc", 
             "d_stem_gc", "anticodon_stem_gc", "t_stem_gc", 
             "anticodon_loop_gc", "t_loop_gc")

# Remove NAs for PCA
pca_data <- species_means[complete.cases(species_means[, ..gc_cols])]

# Run PCA
pca <- prcomp(pca_data[, ..gc_cols], scale = TRUE)
summary(pca)

pca$sdev^2 / sum(pca$sdev^2)

# Merge Order from data
pca_data <- merge(pca_data, data[, c("ID", "Order")], by.x = "accession", by.y = "ID", all.x = TRUE)

# Plot by Order
par(mfrow = c(1, 2))

# Color by Order (top orders only, others grey)
top_orders <- names(sort(table(pca_data$Order), decreasing = TRUE))
order_col <- ifelse(pca_data$Order %in% top_orders, as.numeric(factor(pca_data$Order)), NA)
order_col[is.na(order_col)] <- "grey80"

plot(pca$x[,1], pca$x[,2], 
     col = order_col, pch = 16, cex = 0.5,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "Colored by Order")
legend("topright", legend = top_orders, col = 1:8, pch = 16, cex = 0.6)

# Color by temperature
temp_col <- colorRampPalette(c("blue", "white", "red"))(100)
temp_idx <- cut(pca_data$temp, breaks = 100, labels = FALSE)
plot(pca$x[,1], pca$x[,5],
     col = temp_col[temp_idx], pch = 16, cex = 0.5,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"),
     main = "Colored by Temperature")
par(mfrow=c(1,1))

dinuc_cols <- grep("dinucleotide", colnames(data))
summary(lm(pheno_wc2.1_2.5m_bio_8_p50 ~ .,
           data[,c("pheno_wc2.1_2.5m_bio_8_p50", dinuc_cols)]))


