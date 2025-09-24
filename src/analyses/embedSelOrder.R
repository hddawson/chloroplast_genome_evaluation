library(data.table)
library(arrow)


rmse     <- function(a,b) sqrt(mean((a - b)^2))


embeds <- readRDS("data/tmp/embeds_with_mds.rds")
table(embeds$ManualOutlier)
dim(embeds)
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)

df <- read_parquet("data/processed_data.parquet")
setDT(df)

# variable selection: pick embedding dims by mean correlation across Orders
setDT(clean_embeds)

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
orders <- unique(df$Order)
n_top <- 1
genes <- unique(clean_embeds$Gene)
sel_list <- vector("list", length(genes))
names(sel_list) <- genes

# Pre-merge data once outside the loop
merged_data <- merge(clean_embeds[, c("ID", "Gene", embed_cols), with = FALSE], 
                     df[, .(ID, Order, pheno = get(pheno_col))], 
                     by = "ID")

# Pre-filter genes with sufficient samples
gene_counts <- merged_data[, .N, by = Gene]
valid_genes <- gene_counts[N >= 5]$Gene
stopifnot(length(valid_genes) > 0)

# Vectorized correlation computation
compute_gene_cors <- function(gene_data, embed_cols) {
  # Split by order for correlation computation
  order_groups <- split(gene_data, gene_data$Order)
  
  # Compute correlations for all embedding columns at once per order
  cors_by_order <- lapply(order_groups, function(group) {
    # Vectorized correlation: cor() can handle matrices
    embed_matrix <- as.matrix(group[, embed_cols, with = FALSE])
    cors <- cor(embed_matrix, group$pheno, use = "complete.obs")[, 1]
    return(cors)
  })
  
  # Average correlations across orders
  cors_matrix <- do.call(rbind, cors_by_order)
  mean_cors <- colMeans(cors_matrix, na.rm = TRUE)
  names(mean_cors) <- embed_cols
  
  return(mean_cors)
}

# Main loop - now much faster
sel_list <- vector("list", length(valid_genes))
names(sel_list) <- valid_genes

for (i in seq_along(valid_genes)) {
  cat(gene)
  gene <- valid_genes[i]

  gdt <- merged_data[Gene == gene]

  # Fast correlation computation
  cors_by_order <- compute_gene_cors(gdt, embed_cols)
  
  # Optional: still create histogram
  #hist(cors_by_order, main = gene)
  
  # Select top dimensions
  available <- names(cors_by_order)[!is.na(cors_by_order)]
  k <- min(n_top, length(available))
  top_dims <- names(sort(abs(cors_by_order[available]), decreasing = TRUE))[1:k]
  
  sel <- gdt[, c("ID", top_dims), with = FALSE]
  sel <- unique(sel)  # Remove potential duplicates
  setnames(sel, old = top_dims, new = paste0(gene, "__", top_dims))
  sel_list[[gene]] <- sel
}

# Remove any NULL entries
sel_list <- sel_list[!sapply(sel_list, is.null)]
# drop skipped genes
sel_list <- sel_list[!sapply(sel_list, is.null)]

# merge across all IDs (full outer join)
combined <- Reduce(function(a, b) merge(a, b, by = "ID", all = TRUE), sel_list)

# bring phenotype back
combined <- merge(combined, df[, .(ID, pheno = get(pheno_col))], by = "ID", all.x = TRUE)

# median impute missing predictors
pred_cols <- setdiff(names(combined), c("ID", "pheno"))
for (col in pred_cols) {
  med <- median(combined[[col]], na.rm = TRUE)
  combined[is.na(get(col)), (col) := med]
}

cat("Samples in combined matrix:", nrow(combined), "\n")
cat("Number of predictors:", length(pred_cols), "\n")

pred_cols <- grep("embedding", colnames(combined), value=T)
formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))

fit <- lm(as.formula(formula_str), data=combined)
summary(fit)
plot(fit$fitted.values, combined$pheno)

poaceae_ids <- df[grepl("Poales", Taxonomy), ID]

train <- combined[!ID %in% poaceae_ids]
test  <- combined[ID %in% poaceae_ids]

fit2 <- lm(train$pheno ~ ., data = train[, ..pred_cols])
summary(fit2)
pred2 <- predict(fit2, newdata = test)

plot(pred2, test$pheno, main= "bio8 ~ top 1 emb/gene/order, holdout Poales")

cor(pred2, test$pheno, method="spearman")
abline(a=0,b=1, col="coral")
text(30, 5,
     paste0("spearman=",round(cor(pred2,test$pheno,method="spearman"), 3)),col="red")
text(30,2,
     paste0("RMSE=",round(rmse(pred2,test$pheno), 3)),col="red")

LSD::heatscatter(pred2, test$pheno)

library(ape)
tree <- read.tree("tree/T2.raxml.rba.raxml.lastTree.TMP")

head(tree$tip.label)

