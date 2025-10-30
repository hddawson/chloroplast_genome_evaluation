library(data.table)
library(arrow)
library(ape)

rmse <- function(a,b) sqrt(mean((a - b)^2))

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

df <- read_parquet("data/processed_data.parquet")
setDT(df)

plot(df$pheno_Topt_site_p50,df$pheno_wc2.1_2.5m_bio_8_p50)

unique_orders <- unique(df$Order)
pheno_cols <- grep("pheno", colnames(df), value=TRUE)
pheno_cols <- setdiff(pheno_cols, c("pheno_n_occurrences"))

dir.create("fits", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

compute_gene_cors_mean <- function(gd, embed_cols) {
  em <- as.matrix(gd[, embed_cols, with = FALSE])
  cors <- cor(em, gd$pheno, use="complete.obs")[,1]
  names(cors) <- embed_cols
  cors
}

res_dt <- data.table()

for(pheno_col in pheno_cols){
  cat("Running for", pheno_col, "\n")
  for(holdout_order in unique_orders){
    train_df_ho <- df[Order != holdout_order]
    val_df_ho   <- df[Order == holdout_order]
    train_embeds_ho <- clean_embeds[ID %in% train_df_ho$ID]
    val_embeds_ho   <- clean_embeds[ID %in% val_df_ho$ID]
    
    merged_train_base <- merge(train_embeds_ho[, c("ID","Gene",embed_cols), with=FALSE],
                               train_df_ho[, .(ID, Order, pheno = get(pheno_col))], by="ID")
    gene_counts <- merged_train_base[, .N, by=Gene]
    valid_genes <- gene_counts[N >= 5]$Gene
    if(length(valid_genes)==0) next
    
    n_top <- 2
    m <- "mean"
    sel_list <- vector("list", length(valid_genes)); names(sel_list) <- valid_genes
    for(g in valid_genes){
      gdt <- merged_train_base[Gene == g]
      if(nrow(gdt) < 5) next
      corss <- compute_gene_cors_mean(gdt, embed_cols)
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
    
    fit_file <- sprintf("fits/phenos/fit_%s_top%d_order_%s_%s.rds", m, n_top, make.names(holdout_order), pheno_col)
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
    
    png(sprintf("plots/phenos/otest_%s_%s_top%d_%s.png", make.names(holdout_order), m, n_top, pheno_col), width=800, height=600)
    plot(pred_val, combined_val$pheno, main=sprintf("Test: %s %s top%d %s", make.names(holdout_order), m, n_top, pheno_col),
         xlab="Predicted", ylab="Observed")
    abline(0,1, col="coral")
    legend("topleft", legend=c(paste0("Pearson=", round(pear,3)), paste0("Spearman=", round(spea,3)), paste0("RMSE=", round(test_rmse,3))), bty="n")
    dev.off()
    
    res_dt <- rbind(res_dt, data.table(order=holdout_order, method=m, n_top=n_top,
                                       phenotype=pheno_col,
                                       train_R2=summary(fit)$r.squared,
                                       val_Pearson=pear, val_Spearman=spea,
                                       val_RMSE=test_rmse, fit_file=fit_file))
  }
}

fwrite(res_dt, "results/varsel_grid_results_orders_allphenos.csv")


#yay!

file_dir <- "results/varSelPhenoScan/predictions/"
files <-  list.files(file_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
cors <-c()
for (file in files) {cors <- append(cors, cor(read.csv(file)$y_true,read.csv(file)$y_pred))}

max(cors)

library(data.table)

file_dir <- "results/varSelPhenoScan/predictions/"
files <- list.files(file_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

perf_dt <- data.table()

for (file in files) {
  dat <- fread(file)
  
  # compute correlations
  pear <- cor(dat$y_true, dat$y_pred, method = "pearson",  use = "complete.obs")
  spea <- cor(dat$y_true, dat$y_pred, method = "spearman", use = "complete.obs")
  rmse_val <- sqrt(mean((dat$y_true - dat$y_pred)^2, na.rm = TRUE))
  
  # parse order, phenotype, percentile from filename
  fname <- basename(file)
  # example: preds_mean_top2_order_Apiales_pheno_Topt_site_p10_seed1.csv
  m <- regexec("preds_mean_top2_order_(.+?)_pheno_(.+?)_(p\\d+)_seed\\d+\\.csv", fname)
  parts <- regmatches(fname, m)[[1]]
  holdout_order <- parts[2]
  pheno_col <- parts[3]
  percentile <- parts[4]
  
  if (is.na(holdout_order)) {print(fname)}
  
  perf_dt <- rbind(
    perf_dt,
    data.table(order = holdout_order,
               phenotype = pheno_col,
               percentile = percentile,
               val_Pearson = pear,
               val_Spearman = spea,
               val_RMSE = rmse_val,
               file = fname)
  )
}

# Summarize performance by order × phenotype × percentile
perf_summary <- perf_dt[, .(
  mean_Pearson  = mean(val_Pearson, na.rm = TRUE),
  mean_Spearman = mean(val_Spearman, na.rm = TRUE),
  mean_RMSE     = mean(val_RMSE, na.rm = TRUE),
  n_models      = .N
), by = .(order, phenotype, percentile)]

fwrite(perf_dt,     "results/perf_by_file.csv")
fwrite(perf_summary,"results/perf_summary_by_order_pheno_percentile.csv")

print(perf_summary)

table(perf_summary$order)
table(perf_summary$phenotype)

barplot(mean_Pearson ~ phenotype + order, data=perf_summary,
        subset=perf_summary$percentile=="p10",
        beside=T,las=2)


