library(data.table)
library(arrow)
library(ape)

rmse <- function(a,b) sqrt(mean((a - b)^2))

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

df <- read_parquet("data/processed_data.parquet")
plot(scale(df$pheno_wc2.1_2.5m_bio_8_p50),
     scale(df$pheno_wc2.1_2.5m_bio_13_p50),main="Mean temp vs Mean precipitation, wettest quarter",
     xlab="scaled bio8", ylab="scaled bio13")
cor(df$pheno_wc2.1_2.5m_bio_8_p50,df$pheno_wc2.1_2.5m_bio_13_p50)
abline(a=0,b=1)
setDT(df)
LSD::heatscatter(df$pheno_wc2.1_2.5m_bio_8_p50,
     df$pheno_wc2.1_2.5m_bio_13_p50,main="Mean temp vs Mean precipitation, wettest quarter",
     xlab=" bio8", ylab=" bio13")

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

unique_orders <- unique(df$Order)
dir.create("fits", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

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

grid_n_top <- c(1,2,5)
methods <- c("mean","by_order")
res_dt <- data.table()

for(holdout_order in unique_orders){
  cat(holdout_order, "\n")
  train_df_ho <- df[Order != holdout_order]
  val_df_ho   <- df[Order == holdout_order]
  train_embeds_ho <- clean_embeds[ID %in% train_df_ho$ID]
  val_embeds_ho   <- clean_embeds[ID %in% val_df_ho$ID]
  
  merged_train_base <- merge(train_embeds_ho[, c("ID","Gene",embed_cols), with=FALSE],
                             train_df_ho[, .(ID, Order, pheno = get(pheno_col))], by="ID")
  gene_counts <- merged_train_base[, .N, by=Gene]
  valid_genes <- gene_counts[N >= 5]$Gene
  if(length(valid_genes)==0) next
  
  for(n_top in grid_n_top){
    for(m in methods){
      sel_list <- vector("list", length(valid_genes)); names(sel_list) <- valid_genes
      for(g in valid_genes){
        gdt <- merged_train_base[Gene == g]
        if(nrow(gdt) < 5) next
        corss <- switch(m,
                        mean = compute_gene_cors_mean(gdt, embed_cols),
                        by_order = compute_gene_cors_by_order(gdt, embed_cols))
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
      combined_train <- merge(combined_train, train_df_ho[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
      pred_cols <- setdiff(names(combined_train), c("ID","pheno"))
      if(length(pred_cols)==0) next
      train_medians <- lapply(pred_cols, function(col) median(combined_train[[col]], na.rm=TRUE))
      names(train_medians) <- pred_cols
      for(col in pred_cols) combined_train[is.na(get(col)), (col) := train_medians[[col]]]
      
      formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))
      fit <- lm(as.formula(formula_str), data=combined_train)
      
      fit_file <- sprintf("fits/fit_%s_top%d_order_%s.rds", m, n_top, make.names(holdout_order))
      saveRDS(fit, fit_file)
      
      val_wide_list <- list()
      vc <- merge(val_embeds_ho[, c("ID","Gene",embed_cols), with=FALSE], val_df_ho[, .(ID, pheno = get(pheno_col))], by="ID")
      for(g in unique(vc$Gene)){
        gd <- vc[Gene == g, c("ID", embed_cols), with=FALSE]
        if(nrow(gd)==0) next
        setnames(gd, old = embed_cols, new = paste0(g, "__", embed_cols))
        val_wide_list[[g]] <- unique(gd)
      }
      if(length(val_wide_list)==0) next
      combined_val <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), val_wide_list)
      combined_val <- merge(combined_val, val_df_ho[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
      
      for(col in pred_cols){
        if(!(col %in% names(combined_val))) combined_val[, (col) := train_medians[[col]]]
        combined_val[is.na(get(col)), (col) := train_medians[[col]]]
      }
      
      pred_val <- predict(fit, newdata = combined_val)
      pear <- if(all(is.na(pred_val)) || all(is.na(combined_val$pheno))) NA else cor(pred_val, combined_val$pheno, method="pearson", use="complete.obs")
      spea <- if(all(is.na(pred_val)) || all(is.na(combined_val$pheno))) NA else cor(pred_val, combined_val$pheno, method="spearman", use="complete.obs")
      test_rmse <- if(all(is.na(pred_val)) || all(is.na(combined_val$pheno))) NA else rmse(pred_val, combined_val$pheno)
      
      png(sprintf("plots/otest_%s_%s_top%d.png", make.names(holdout_order), m, n_top), width=800, height=600)
      plot(pred_val, combined_val$pheno, main=sprintf("Test: %s %s top%d", make.names(holdout_order), m, n_top), xlab="Predicted", ylab="Observed")
      abline(0,1, col="coral")
      legend("topleft", legend=c(paste0("Pearson=", round(pear,3)), paste0("Spearman=", round(spea,3)), paste0("RMSE=", round(test_rmse,3))), bty="n")
      dev.off()
      
      res_dt <- rbind(res_dt, data.table(order=holdout_order, method=m, n_top=n_top,
                                         train_R2=summary(fit)$r.squared,
                                         val_Pearson=pear, val_Spearman=spea,
                                         val_RMSE=test_rmse, fit_file=fit_file))
    }
  }
}

fwrite(res_dt, "results/varsel_grid_results_orders.csv")
res_dt <- read.csv("results/varsel_grid_results_orders.csv")

head(res_dt)

hist(res_dt$n_top)
hist(res_dt$train_R2)
hist(res_dt$val_Pearson)
hist(res_dt$val_Spearman)
hist(res_dt$val_RMSE)

ggplot(res_dt, aes(x=factor(n_top), y=val_Pearson, color=method, group=method)) +
  geom_point(size=3) +
  labs(x="Top N embeddings per gene",
       y="Test Pearson",
       title="Effect of variable selection and n_embeddings on test Pearson") +
  facet_wrap(~method + order) +
  theme_minimal()

barplot(train_R2 ~ order + n_top, data=res_dt, res_dt$method=="mean")
barplot(train_R2 ~ order + n_top, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE)

barplot(val_RMSE ~ order + n_top, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE)

barplot(val_Pearson ~ order + n_top, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE)

barplot(val_Spearman ~ order + n_top, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE)

par(mfrow=c(1,2))
barplot(train_R2 ~ order + n_top, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE,
        main="mean/gene")
barplot(train_R2 ~ order + n_top, data=res_dt,
        subset=res_dt$method=="by_order", beside=TRUE,
        main="mean/gene/order")

par(mfrow=c(1,2))
barplot(val_Spearman ~ order + method, data=res_dt,
        subset=res_dt$method=="mean", beside=TRUE,
        main="mean/gene")
barplot(val_Spearman ~ order + n_top, data=res_dt,
        subset=res_dt$method=="by_order", beside=TRUE,
        main="mean/gene/order")
dev.off()
barplot(val_Spearman ~ order + method, data=res_dt,
        subset=res_dt$n_top==1, beside=TRUE,
        main="mean/gene")

barplot(val_Spearman ~ order + method, data=res_dt,
        subset=res_dt$n_top==2, beside=TRUE,
        main="mean/gene")
barplot(val_Spearman ~ order + method, data=res_dt,
        subset=res_dt$n_top==5, beside=TRUE,
        main="mean/gene")

pca <- prcomp(res_dt[,c("val_Spearman","val_Pearson", "val_RMSE")])
res_fit <- lm(val_Spearman ~ method + n_top, data=res_dt)
plot(res_fit)
summary(res_fit)
mean(res_dt$val_Spearman)

res_fit <- lm(val_Spearman ~ method + n_top + order, data=res_dt)
plot(res_fit)
summary(res_fit)

plot(res_dt[which(res_dt$method=="by_order"),"val_Spearman"],
     res_dt[which(res_dt$method=="mean"),"val_Spearman"])
line(res_dt[which(res_dt$method=="by_order"),"val_Spearman"],
     res_dt[which(res_dt$method=="mean"),"val_Spearman"])
plot(res_dt[which(res_dt$method=="by_order"),"val_Pearson"],
     res_dt[which(res_dt$method=="mean"),"val_Pearson"])

df <- subset(res_dt, method == "mean")
# make a table of val_Spearman by n_top and order
tab <- with(df, tapply(val_Spearman, list(n_top, order), mean))
# barplot expects matrix input
barplot(tab, beside = TRUE, legend.text = F,
        las = 2)
abline(h = 0.3, col = "red", lty = 2) 

plot(res_dt$val_Pearson, res_dt$val_Spearman)

hist(res_dt$val_Pearson)

barplot(val_Spearman ~ n_top + order, data=res_dt, subset=res_dt$method=="by_order",
        beside=T,las=2)
barplot(val_Spearman ~ n_top + order, data=res_dt, subset=res_dt$method=="by_mean",
        beside=T,las=2)

library(ggplot2)

ggplot(res_dt, aes(x=val_Spearman, y=val_RMSE, color=method, shape=method)) +
  geom_point() +
  facet_grid(order ~ n_top) +
  theme_bw()

plot(res_dt[res_dt$method=="mean",]$val_Spearman,
     res_dt[res_dt$method=="by_order",]$val_Spearman)
abline(a=0,b=1)
plot(res_dt[res_dt$method=="mean",]$val_Pearson,
     res_dt[res_dt$method=="by_order",]$val_Pearson,col=colors()[as.integer(as.factor(res_dt$order))],
     pch=19)

aov_fit <- aov(val_Spearman ~ order * n_top, data=res_dt[res_dt$method=="mean",])
summary(aov_fit)

aggregate(val_Spearman ~ n_top, data=res_dt[res_dt$method=="mean",], mean, na.rm=TRUE)

