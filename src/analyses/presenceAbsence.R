library(data.table)
library(stringr)
data <- fread("data/tmp/rbcL_aln/merged_aa_counts.csv")

binary_cols <- names(data)[sapply(data, function(col) length(unique(col[!is.na(col)])) == 2)]
binary_cols <- binary_cols[!str_detect(binary_cols, "_outlier")]
binary_cols <- binary_cols[1:119]

genes <- unique(unlist(data$Genes))
str(data$Genes)

pheno_cols <- names(data)[grep("pheno",names(data))]

pca <- prcomp(data[,..binary_cols],scale.=TRUE)
summary(pca)
orders <- unique(unlist(str_split(data$Taxonomy,";")))
orders <- orders[grep("ales", orders)]
par(mfrow=c(1,2))
plot(pca$x[,1],pca$x[,2], main= "Presence Absence PCs")
points(pca$x[grep("Brassicales", data$Taxonomy),1],pca$x[grep("Brassicales", data$Taxonomy),2],
       col="seagreen")
points(pca$x[grep("Liliales", data$Taxonomy),1],pca$x[grep("Liliales", data$Taxonomy),2],
       col="yellow")
barplot(pca$sdev^2 * 100 / sum(pca$sdev^2),main="% Variance explained")


loadings <- pca$rotation

top_pc1 <- sort(loadings[, "PC1"], decreasing = TRUE)[1:10]
top_pc2 <- sort(loadings[, "PC2"], decreasing = TRUE)[1:10]

aa_pcs <- readRDS("data/tmp/aa_supermatrix_expanded_3pcs_pca.rds")
identical(rownames(aa_pcs$x), rownames(pca$x))

cor(aa_pcs$x[,2],pca$x[,2])

cor_mat <- cor(aa_pcs$x[,1:5], pca$x[,1:5])
round(cor_mat, 3)

pheatmap(
  cor_mat,
  display_numbers = TRUE,         # show numbers in cells
  number_format = "%.2f",         # two decimals
  cluster_rows = FALSE,           # keep PC1–PC5 order
  cluster_cols = FALSE,
  color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(100),
  main = "Correlation between AA PCs and Presence/Absence PCs"
)

results <- list()

for (pheno in pheno_cols) {
  
  # escape column names safely
  safe_terms <- paste0("`", binary_cols, "`")
  
  fmla <- as.formula(
    paste("`", pheno, "` ~ ", paste(safe_terms, collapse = " + "), sep = "")
  )
  
  # fit linear model
  fit <- lm(fmla, data = data)
  
  # extract tidy stats
  tidy_stats <- tidy(fit) %>%
    mutate(Phenotype = pheno)
  
  # model fit (R² etc.)
  gl <- glance(fit) %>%
    mutate(Phenotype = pheno)
  
  # store
  results[[pheno]] <- list(coeffs = tidy_stats, fit = gl)
}

# combine into big dataframes
all_coeffs <- bind_rows(lapply(results, `[[`, "coeffs"))
all_fits   <- bind_rows(lapply(results, `[[`, "fit"))

# Example: top signals by phenotype
all_coeffs %>%
  filter(term != "(Intercept)") %>%
  arrange(p.value) %>%
  head(20)

hist(all_fits$r.squared, main="Fit of WC Bio ~ Gain/loss over angiosperms")

results_list <- list()

for (ord in orders) {
  # subset df to samples in this Order
  subdf <- data %>% filter(grepl(ord, Taxonomy))
  
  # skip if too small
  if (nrow(subdf) < 20) next   # adjust cutoff as needed
  
  message("Fitting models for ", ord, " with ", nrow(subdf), " samples")
  
  for (pheno in pheno_cols) {
    # build safe formula
    safe_terms <- paste0("`", binary_cols, "`")
    
    fmla <- as.formula(
      paste("`", pheno, "` ~ ", paste(safe_terms, collapse = " + "), sep = "")
    )
    
    fit <- lm(fmla, data = subdf)
    
    # extract tidy results
    tidy_res <- broom::tidy(fit) %>%
      filter(term != "(Intercept)") %>%
      mutate(Phenotype = pheno,
             Order = ord,
             Rsq = summary(fit)$r.squared,
             Adj_Rsq = summary(fit)$adj.r.squared,
             N = nrow(subdf))
    
    results_list[[length(results_list) + 1]] <- tidy_res
  }
}

all_results <- bind_rows(results_list)

par(mfrow=c(1,2))
hist(all_fits$r.squared, main="Fits of WC_Bio ~ Gain/loss across Orders")
hist(all_results$Rsq[all_results$p.value < 0.05], main="Fits of WC_Bio ~ Gain/loss within Orders")

#let's fit order as a covariate
orders <- orders[-58] #drop Halesia
pattern <- paste(orders, collapse = "|")  # regex OR pattern
data$Order <- str_extract(data$Taxonomy, pattern)
table(data$Order)

library(broom.mixed)
library(lme4)
data$Order <- as.factor(data$Order)

results <- list()

for (pheno in pheno_cols) {
  
  # escape column names safely
  safe_terms <- paste0("`", binary_cols, "`")
  
  # fit linear model
  fit <- lmer(as.formula(
    paste0("`", pheno, "` ~ ", paste(safe_terms, collapse = " + "), " + (1|Order)")
  ), data = data)
  
  # extract tidy stats (fixed effects only)
  tidy_stats <- tidy(fit, effects = "fixed") %>%
    mutate(Phenotype = pheno)
  
  # model fit (R² etc.)
  gl <- glance(fit) %>%
    mutate(Phenotype = pheno)
  
  # store results
  results[[pheno]] <- list(coeffs = tidy_stats, fit = gl)
}

all_coeffs <- bind_rows(lapply(results, `[[`, "coeffs"))
all_fits   <- bind_rows(lapply(results, `[[`, "fit"))

gene_summary <- all_coeffs %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  summarise(
    mean_est = mean(estimate, na.rm = TRUE),
    median_est = median(estimate, na.rm = TRUE),
    prop_sig = mean(abs(statistic) > 1.96, na.rm = TRUE),  # ~ p<0.05
    n_sig = sum(abs(statistic) > 1.96, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_sig))

head(gene_summary, 15)

ggplot(all_coeffs %>% filter(term != "(Intercept)"),
       aes(x = Phenotype, y = term, fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Gene Effects Across Phenotypes",
       x = "Phenotype", y = "Gene")

fit_example <- results[["pheno_Topt_site_p10"]]$fit
VarCorr(fit_example)

top10_pheno <- all_results %>%
  group_by(Phenotype) %>%
  summarise(
    max_Rsq = max(Rsq, na.rm = TRUE),
    avg_Rsq = mean(Rsq, na.rm = TRUE),
    n = length(Rsq)
  ) %>%
  arrange(desc(max_Rsq)) %>%
  slice_head(n = 10)

top10_pheno

tech_means <- tapply(data$geneCount, data$SequencingTech, mean, na.rm = TRUE)

par(mfrow=c(1,1))
library(ggplot2)
tech_counts <- data[, .N, by = SequencingTech]
long_dt <- melt(
  data,
  id.vars = "SequencingTech",
  measure.vars = binary_cols,
  variable.name = "Gene",
  value.name = "Presence"
)

gc_dt <- melt(
  data,
  id.vars = "SequencingTech",
  measure.vars = "genomicGC")
ggplot(gc_dt, aes(x = SequencingTech, y = value, fill = SequencingTech)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # legend repeats across facets, remove
  ) +
  labs(
    title = "genomic GC across sequencing technologies",
    x = "Sequencing Technology",
    y = "GC content"
  ) +
  geom_text(data = tech_counts,
            aes(x = SequencingTech,
                y = max(gc_dt$value, na.rm = TRUE) * 1.05,
                label = paste0("n=", N)),
            inherit.aes = FALSE, vjust = 0)



ggplot(data, aes(x = SequencingTech, y = ycf4, fill = SequencingTech)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3) +
  geom_text(data = tech_counts,
            aes(x = SequencingTech,
                y = max(data$GeneCount, na.rm = TRUE) * 1.05,
                label = paste0("n=", N)),
            inherit.aes = FALSE, vjust = 0) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Average Gene Count by Sequencing Technology",
       x = "Sequencing Technology", y = "Mean Gene Count")

p <- ggplot(long_dt, aes(x = SequencingTech, y = Presence, fill = SequencingTech)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # legend repeats across facets, remove
  ) +
  labs(
    title = "Gene presence across sequencing technologies",
    x = "Sequencing Technology",
    y = "Proportion present"
  ) +
  facet_wrap(~Gene, scales = "free_y")  # one panel per gene

# Save to PDF
ggsave("plots/gene_presence_by_tech.pdf", p, width = 20, height = 25, units = "in")

#is there an association with gain/loss and sequencing tech/year
table(data$Year)
table(data$SequencingTech)
head(binary_cols)
colMeans(data[,..binary_cols])


# Calculate mean presence for each gene
gene_means <- data.table(
  Gene = binary_cols,
  MeanPresence = colMeans(data[, ..binary_cols], na.rm = TRUE)
)

# Extract gene prefix
gene_means[, Prefix := substr(Gene, 1, 3)]

# Order genes by mean presence
gene_means <- gene_means[order(MeanPresence)]
gene_means[, Gene := factor(Gene, levels = Gene)]  # keep sort order

n_prefixes <- length(unique(gene_means$Prefix))
prefix_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_prefixes)

p <- ggplot(gene_means, aes(x = Gene, y = MeanPresence, color = Prefix)) +
  geom_point(size = 3) +
  scale_color_manual(values = setNames(prefix_colors, unique(gene_means$Prefix))) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
    legend.title = element_blank()
  ) +
  labs(
    title = "Mean Gene Presence Across Samples",
    x = "Gene",
    y = "Mean Presence (Proportion)"
  )

print(p)





assoc_results_tech <- lapply(binary_cols, function(gene) {
  tab <- table(data[[gene]], data$SequencingTech)
  if (all(dim(tab) > 1)) {
    test <- chisq.test(tab)
    data.frame(Gene = gene, 
               P_value = test$p.value, 
               stringsAsFactors = FALSE)
  } else {
    data.frame(Gene = gene, P_value = NA)  # no variation
  }
}) %>% bind_rows() %>%
  mutate(P_adj = p.adjust(P_value, method = "fdr")) %>%
  arrange(P_adj)
head(assoc_results_tech, 10)

