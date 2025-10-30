#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
  library(arrow)
})

option_list <- list(
  make_option(c("-d","--data-dir"),   type="character", default="data"),
  make_option(c("-o","--output-dir"), type="character", default="results"),
  make_option(c("-c","--cores"),      type="integer",  default=10)
)

opt        <- parse_args(OptionParser(option_list=option_list))
data_dir   <- opt$`data-dir`
output_dir <- opt$`output-dir`
n_cores    <- opt$`cores`
dirs       <- c("predictions","models","plots")
lapply(file.path(output_dir, dirs), dir.create, recursive=TRUE, showWarnings=FALSE)

rmse <- function(a,b) sqrt(mean((a - b)^2))
mae <- function(a,b) mean(abs(a - b))
r2 <- function(a,b) if(all(is.na(a))||all(is.na(b))) NA else cor(a,b)^2
spearman <- function(a,b) if(all(is.na(a))||all(is.na(b))) NA else cor(a,b, method="spearman")
median_impute <- function(x) { x[is.na(x)] <- median(x, na.rm=TRUE); x }

embeds <- as.data.table(read_parquet(file.path(data_dir, "clean_embeds.parquet")))
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds
data_all <- as.data.table(read_parquet(file.path(data_dir, "processed_data.parquet")))
setDT(data_all)
pheno_cols <- grep("pheno", colnames(data_all), value=TRUE)
pheno_cols <- setdiff(pheno_cols, c("pheno_n_occurrences"))
unique_orders <- unique(data_all$Order)

sum(is.na(data_all[,..pheno_cols])) #0

compute_gene_cors_mean <- function(gd, embed_cols) {
  em <- as.matrix(gd[, embed_cols, with = FALSE])
  cors <- cor(em, gd$pheno, use="complete.obs")[,1]
  names(cors) <- embed_cols
  cors
}

n_top <- 2
method_name <- "mean"

holdout_strategies <- list(
  list(name="order_holdout", get_idx=function(order_name) which(data_all$Order==order_name))
)

seeds <- c(1)
grid <- CJ(seed=seeds, pheno=pheno_cols, holdout_idx=seq_along(unique_orders), stringsAsFactors=FALSE)
grid[, task := .I]

run_task <- function(row_index) {
  tryCatch({
    row <- grid[row_index]
    seed <- row$seed
    pheno_col <- row$pheno
    holdout_order <- unique_orders[row$holdout_idx]
    cat("\n⏳ raunning loh",holdout_order,pheno_col,"\n")
    set.seed(seed)
    train_df_ho <- data_all[Order != holdout_order]
    val_df_ho   <- data_all[Order == holdout_order]
    train_embeds_ho <- clean_embeds[ID %in% train_df_ho$ID]
    val_embeds_ho   <- clean_embeds[ID %in% val_df_ho$ID]
    merged_train_base <- merge(train_embeds_ho[, c("ID","Gene",embed_cols), with=FALSE],
                               train_df_ho[, .(ID, Order, pheno = get(pheno_col))], by="ID")
    gene_counts <- merged_train_base[, .N, by=Gene]
    valid_genes <- gene_counts[N >= 5, Gene]
    if(length(valid_genes)==0) return(NULL)
    sel_list <- vector("list", length(valid_genes)); names(sel_list) <- valid_genes
    for(g in valid_genes){
      gdt <- merged_train_base[Gene==g]
      if(nrow(gdt)<5) next
      corss <- compute_gene_cors_mean(gdt, embed_cols)
      available <- names(corss)[!is.na(corss)]
      if(length(available)==0) next
      k <- min(n_top, length(available))
      top_dims <- names(sort(abs(corss[available]), decreasing=TRUE))[1:k]
      sel <- unique(gdt[, c("ID", top_dims), with=FALSE])
      setnames(sel, old=top_dims, new=paste0(g,"__",top_dims))
      sel_list[[g]] <- sel
    }
    sel_list <- sel_list[!sapply(sel_list, is.null)]
    if(length(sel_list)==0) return(NULL)
    combined_train <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), sel_list)
    combined_train <- merge(combined_train, train_df_ho[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
    pred_cols <- setdiff(names(combined_train), c("ID","pheno"))
    if(length(pred_cols)==0) return(NULL)
    train_medians <- lapply(pred_cols, function(col) median(combined_train[[col]], na.rm=TRUE))
    names(train_medians) <- pred_cols
    for(col in pred_cols) combined_train[is.na(get(col)), (col) := train_medians[[col]]]
    formula_str <- paste("pheno ~", paste(pred_cols, collapse=" + "))
    fit <- tryCatch(lm(as.formula(formula_str), data=combined_train), error=function(e) NULL)
    if(is.null(fit)) return(NULL)
    fit_file <- file.path(output_dir,"models",paste0("fit_",method_name,"_top",n_top,"_order_",make.names(holdout_order),"_",pheno_col,"_seed",seed,".rds"))
    coefs <- summary(fit)$coefficients
    saveRDS(list(
      coef = coefs[,1],
      se   = coefs[,2],
      tval = coefs[,3],
      pval = coefs[,4]
    ), fit_file)
    vc <- merge(val_embeds_ho[, c("ID","Gene",embed_cols), with=FALSE], val_df_ho[, .(ID, pheno = get(pheno_col))], by="ID")
    val_wide_list <- list()
    for(g in unique(vc$Gene)){
      gd <- vc[Gene==g, c("ID", embed_cols), with=FALSE]
      if(nrow(gd)==0) next
      setnames(gd, old=embed_cols, new=paste0(g,"__",embed_cols))
      val_wide_list[[g]] <- unique(gd)
    }
    val_wide_list <- val_wide_list[!sapply(val_wide_list, is.null)]
    if(length(val_wide_list)==0) return(NULL)
    combined_val <- Reduce(function(a,b) merge(a,b,by="ID",all=TRUE), val_wide_list)
    combined_val <- merge(combined_val, val_df_ho[, .(ID, pheno = get(pheno_col))], by="ID", all.x=TRUE)
    test_medians <- sapply(pred_cols, function(col) median(combined_val[[col]], na.rm=TRUE), simplify=TRUE)
    for(col in pred_cols){
      if(!(col %in% names(combined_val))) combined_val[, (col) := test_medians[col]]
      combined_val[is.na(get(col)), (col) := test_medians[col]]
    }
    pred_val <- predict(fit, newdata = combined_val)
    truth <- combined_val$pheno
    ids <- combined_val$ID
    pred_df <- data.frame(ID=ids, y_true=truth, y_pred=pred_val)
    write.csv(pred_df, file.path(output_dir,"predictions",paste0("preds_",method_name,"_top",n_top,"_order_",make.names(holdout_order),"_",pheno_col,"_seed",seed,".csv")), row.names=FALSE)
    metrics <- data.frame(combo=paste0("holdout=",holdout_order,"_pheno=",pheno_col,"_seed=",seed),
                          seed=seed, holdout=holdout_order, phenotype=pheno_col,
                          method=method_name, n_top=n_top,
                          n_test=sum(!is.na(truth)),
                          RMSE=if(all(is.na(pred_val))||all(is.na(truth))) NA else rmse(truth,pred_val),
                          MAE=if(all(is.na(pred_val))||all(is.na(truth))) NA else mae(truth,pred_val),
                          R2=if(all(is.na(pred_val))||all(is.na(truth))) NA else r2(truth,pred_val),
                          Spearman=if(all(is.na(pred_val))||all(is.na(truth))) NA else spearman(truth,pred_val),
                          fit_file=fit_file, stringsAsFactors=FALSE)
    png(file.path(output_dir,"plots",paste0("plot_",method_name,"_top",n_top,"_order_",make.names(holdout_order),"_",pheno_col,"_seed",seed,".png")),600,600)
    plot(pred_val, truth, pch=19, cex=1, xlab="Predicted", ylab="Observed", main=paste(pheno_col, holdout_order))
    abline(0,1,col="blue",lwd=2); dev.off()
    metrics}, error = function(e) {
    cat("Task", row_index, "failed:", conditionMessage(e), "\n")
    return(NULL)
  })
}
results <- mclapply(grid$task, run_task, mc.cores=n_cores)
results <- results[!sapply(results, is.null)]
summary_df <- rbindlist(results, fill=TRUE)
write.csv(summary_df, file.path(output_dir,"embeds_varsel_grid_summary.csv"), row.names=FALSE)
cat("Completed. Wrote summary to", file.path(output_dir,"embeds_varsel_grid_summary.csv"), "\n")


res <- read.csv("results/varSelPhenoScan/embeds_varsel_grid_summary.csv")
table(res$holdout)

table(res$pheno)

clean_embeds$ma
