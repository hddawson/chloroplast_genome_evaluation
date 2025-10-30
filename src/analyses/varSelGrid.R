library(data.table)
library(arrow)
library(ape)

rmse <- function(a,b) sqrt(mean((a - b)^2))

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

df <- read_parquet("data/processed_data.parquet")
setDT(df)

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
poaceae_ids <- df[grepl("Poales", Taxonomy), ID]
train_df <- df[!ID %in% poaceae_ids]
test_df  <- df[ ID %in% poaceae_ids]
train_embeds <- clean_embeds[ID %in% train_df$ID]
test_embeds  <- clean_embeds[ID %in% test_df$ID]

compute_gene_cors_mean <- function(gd, embed_cols) {
  em <- as.matrix(gd[, embed_cols, with = FALSE])
  cors <- cor(em, gd$pheno, use="complete.obs")[,1]
  names(cors) <- embed_cols
  cors
}
compute_gene_cors_by_order <- function(gd, embed_cols) {
  og <- split(gd, gd$Order)
  cors_by_order <- lapply(og, function(g){
    em <- as.matrix(g[, embed_cols, with = FALSE])
    cor(em, g$pheno, use="complete.obs")[,1]
  })
  cors_mat <- do.call(rbind, cors_by_order)
  cm <- colMeans(cors_mat, na.rm=TRUE)
  names(cm) <- embed_cols
  cm
}

compute_gene_cors_by_order_weighted <- function(gd, embed_cols) {
  og <- split(gd, gd$Order)
  cors_by_order <- lapply(og, function(g){
    em <- as.matrix(g[, embed_cols, with = FALSE])
    cor(em, g$pheno, use="complete.obs")[,1]
  })
  n_per_order <- sapply(og, nrow)
  cors_mat <- do.call(rbind, cors_by_order)
  cm <- colSums(cors_mat * n_per_order, na.rm=TRUE) / colSums(!is.na(cors_mat) * n_per_order)
  names(cm) <- embed_cols
  cm
}


merged_train_base <- merge(train_embeds[, c("ID","Gene",embed_cols), with=FALSE],
                           train_df[, .(ID, Order, pheno = get(pheno_col))], by="ID")
gene_counts <- merged_train_base[, .N, by=Gene]
valid_genes <- gene_counts[N >= 5]$Gene

dir.create("plots", showWarnings = FALSE)
grid_n_top <- c(1,2,3,4,5,6,7,8,9,10)
results <- list()
res_dt <- data.table()

methods <- c("mean","by_order","by_order_weighted")

for(n_top in grid_n_top){
  for(m in methods){
    sel_list <- vector("list", length(valid_genes)); names(sel_list) <- valid_genes
    for(g in valid_genes){
      gdt <- merged_train_base[Gene == g]
      if(nrow(gdt) < 5) next
      corss <- switch(m,
                      mean = compute_gene_cors_mean(gdt, embed_cols),
                      by_order = compute_gene_cors_by_order(gdt, embed_cols),
                      by_order_weighted = compute_gene_cors_by_order_weighted(gdt, embed_cols))
      available <- names(corss)[!is.na(corss)]
      if(length(available)==0) next
      k <- min(n_top, length(available))
      top_dims <- names(sort(abs(corss[available]), decreasing = TRUE))[1:k]
      sel <- unique(gdt[, c("ID", top_dims), with=FALSE])
      setnames(sel, old = top_dims, new = paste0(g, "__", top_dims))
      sel_list[[g]] <- sel
    }
    sel_list <- sel_list[!sapply(sel_list, is.null)]
    
    if(length(sel_list)==0) next
    combined_train <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), sel_list)
    combined_train <- merge(combined_train, train_df[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
    
    pred_cols <- setdiff(names(combined_train), c("ID","pheno"))
    train_medians <- lapply(pred_cols, function(col) median(combined_train[[col]], na.rm=TRUE))
    names(train_medians) <- pred_cols
    
    for(col in pred_cols) combined_train[is.na(get(col)), (col) := train_medians[[col]]]
    
    formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))
    fit <- lm(as.formula(formula_str), data=combined_train)
    
    train_preds <- fit$fitted.values
    train_r2 <- summary(fit)$r.squared
    
    png(sprintf("plots/train_%s_top%d.png", m, n_top), width=800, height=600)
    plot(train_preds, combined_train$pheno, main=sprintf("Train: %s top%d", m, n_top), xlab="Predicted", ylab="Observed")
    abline(0,1)
    text(mean(train_preds, na.rm=TRUE), min(combined_train$pheno, na.rm=TRUE), paste0("R2=", round(train_r2,3)))
    dev.off()
    
    test_wide_list <- list()
    tc <- merge(test_embeds[, c("ID","Gene",embed_cols), with=FALSE], test_df[, .(ID, pheno = get(pheno_col))], by="ID")
    for(g in unique(tc$Gene)){
      gd <- tc[Gene == g, c("ID", embed_cols), with=FALSE]
      if(nrow(gd)==0) next
      setnames(gd, old = embed_cols, new = paste0(g, "__", embed_cols))
      test_wide_list[[g]] <- unique(gd)
    }
    if(length(test_wide_list)==0) next
    
    combined_test <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), test_wide_list)
    combined_test <- merge(combined_test, test_df[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
    for(col in pred_cols){
      if(!col %in% names(combined_test)) combined_test[, (col) := train_medians[[col]]]
      combined_test[is.na(get(col)), (col) := train_medians[[col]]]
    }
    
    pred_test <- predict(fit, newdata = combined_test)
    pear <- cor(pred_test, combined_test$pheno, method="pearson")
    spea <- cor(pred_test, combined_test$pheno, method="spearman")
    test_rmse <- rmse(pred_test, combined_test$pheno)
    png(sprintf("plots/test_%s_top%d.png", m, n_top), width=800, height=600)
    plot(pred_test, combined_test$pheno, main=sprintf("Test: %s top%d", m, n_top), xlab="Predicted", ylab="Observed")
    abline(0,1, col="coral")
    legend("topleft", legend=c(paste0("Pearson=", round(pear,3)), paste0("Spearman=", round(spea,3)), paste0("RMSE=", round(test_rmse,3))), bty="n")
    dev.off()
    res_dt <- rbind(res_dt, data.table(method=m, n_top=n_top, train_R2=train_r2, test_Pearson=pear, test_Spearman=spea, test_RMSE=test_rmse, train_plot=paste0("plots/train_",m,"_top",n_top,".png"), test_plot=paste0("plots/test_",m,"_top",n_top,".png")))
  }
}

fwrite(res_dt, "results/varsel_grid_results.csv")
res_dt

library(ggplot2)
res <- read.csv("results/varsel_grid_results.csv")
boxplot(test_Pearson ~ method + n_top, data=res)
boxplot(test_Pearson ~ method + n_top, data=res)

ggplot(res, aes(x=factor(n_top), y=test_Pearson, color=method, group=method)) +
  geom_point(size=3) +
  geom_line(aes(linetype=method)) +
  labs(x="Top N embeddings per gene", y="Test Pearson", title="Effect of variable selection and n_embeddings on test performance") +
  theme_minimal()

ggplot(res, aes(x=factor(n_top), y=test_RMSE, color=method, group=method)) +
  geom_point(size=3) +
  geom_line(aes(linetype=method)) +
  labs(x="Top N embeddings per gene", y="Test Pearson", title="Effect of variable selection and n_embeddings on test performance") +
  theme_minimal()

plot(res$test_Spearman, res$test_Pearson)

abline(line(res$test_Spearman, res$test_Pearson))

