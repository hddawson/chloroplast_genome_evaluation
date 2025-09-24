#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(BGLR)
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

rmse     <- function(a,b) sqrt(mean((a - b)^2))
mae      <- function(a,b) mean(abs(a - b))
r2       <- function(a,b) cor(a,b)^2
spearman <- function(a,b) cor(a,b, method="spearman")
median_impute <- function(x) { x[is.na(x)] <- median(x, na.rm=TRUE); x }

embeds <- as.data.table(read_parquet(file.path(data_dir, "clean_embeds.parquet")))
pheno_data <- as.data.table(read_parquet(file.path(data_dir, "pheno_topt_clean.parquet")))
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)

# fast reshape: ID × Gene → wide
#embeds_long <- melt(embeds, id.vars=c("ID","Gene"), measure.vars=embed_cols,
#                    variable.name="dim", value.name="val")
#embeds_long[, feature := paste(Gene, dim, sep="_")]
#species_embeddings <- dcast(embeds_long, ID ~ feature, value.var="val")
#dim(species_embeddings)
#sum(is.na(species_embeddings))
# impute
#embed_matrix_cols <- setdiff(names(species_embeddings), "ID")
#species_embeddings[, (embed_matrix_cols) := lapply(.SD, median_impute), .SDcols=embed_matrix_cols]
#saveRDS(species_embeddings, "data/tmp/species_embeddings.rds")
species_embeddings <- readRDS("data/tmp/species_embeddings.rds")
embed_matrix_cols <- setdiff(names(species_embeddings), "ID")

data <- read_parquet("data/processed_data.parquet")
setDT(data)
final_data <- merge(pheno_data, species_embeddings, by="ID", all.x=TRUE)
final_data <- final_data[!is.na(pheno_Topt_site_p50)]
taxa_df <- data[, .(ID, Taxonomy, Order, pheno_wc2.1_2.5m_bio_8_p50)]
final_data <- merge(final_data, taxa_df, by = "ID", all.x = TRUE)

y <- final_data$pheno_wc2.1_2.5m_bio_8_p50
Z <- as.matrix(final_data[, ..embed_matrix_cols])

holdout_strategies <- list(
  list(name="order", get_idx=function() grep("Poaceae", final_data$Taxonomy)),
  list(name="random", get_idx=function(seed){ set.seed(seed); sample(seq_len(length(y)), size=round(0.1*length(y))) })
)
model_configs <- list(list(name="BayesC", type="marker", model="BayesC", probIn=0.05))
                      #list(name="BL", type="marker", model="BL", probIn=0.05))
data_configs  <- list(list(name="embeds", X=Z, fixed_effects=NULL))
seeds <- c(1,2,3,4)

grid <- expand.grid(seed=seeds,
                    model_idx=seq_along(model_configs),
                    holdout_idx=seq_along(holdout_strategies),
                    data_idx=seq_along(data_configs),
                    stringsAsFactors=FALSE)
grid$task <- seq_len(nrow(grid))

run_task <- function(row_index) {
  row <- grid[row_index,]
  seed <- row$seed
  model_config <- model_configs[[row$model_idx]]
  holdout_cfg <- holdout_strategies[[row$holdout_idx]]
  data_config <- data_configs[[row$data_idx]]
  set.seed(seed)
  
  test_idx <- if (holdout_cfg$name=="random") holdout_cfg$get_idx(seed) else holdout_cfg$get_idx()
  yNA <- y; yNA[test_idx] <- NA
  
  combo_name <- paste0("Holdout-",holdout_cfg$name,"_Data-",data_config$name,
                       "_Model-",model_config$name,"_Seed-",seed)
  cat("\n⏳ Running:", combo_name, "\n")
  
  ETA <- list(list(X=data_config$X, model=model_config$model, probIn=model_config$probIn))
  fit <- BGLR(y=yNA, ETA=ETA, nIter=200000, burnIn=40000, verbose=FALSE,
              saveAt=file.path(output_dir,"models",paste0("fit_",combo_name,"_")))
  saveRDS(fit, file.path(output_dir,"models",paste0("fit_",combo_name,".rds")))
  
  preds <- fit$yHat[test_idx]; truth <- y[test_idx]; ids <- final_data$ID[test_idx]
  write.csv(data.frame(ID=ids, y_true=truth, y_pred=preds),
            file.path(output_dir,"predictions",paste0(combo_name,"_preds.csv")), row.names=FALSE)
  
  metrics <- data.frame(combo=combo_name, seed=seed, holdout=holdout_cfg$name,
                        data_type=data_config$name, model=model_config$name,
                        n_test=length(test_idx),
                        RMSE=rmse(truth,preds), MAE=mae(truth,preds),
                        R2=r2(truth,preds), Spearman=spearman(truth,preds))
  png(file.path(output_dir,"plots",paste0(combo_name,"_plot.png")),600,600)
  plot(preds, truth, pch=19, cex=1.2,
       xlab=paste("Predicted (",holdout_cfg$name," holdout)"),
       ylab="Observed", main=paste(data_config$name,model_config$name))
  abline(0,1,col="blue",lwd=2); dev.off()
  metrics
}

cat("\n▶ Launching grid search with", nrow(grid), "models using", n_cores, "cores\n")
results <- mclapply(grid$task, run_task, mc.cores=n_cores)
summary_df <- rbindlist(results)
write.csv(summary_df, file.path(output_dir,"embeds_grid_results_summary.csv"), row.names=FALSE)
cat("\n✅ Grid search complete!\n")

library(coda)
m1 <- "results/embedsAll10kIts/models/fit_Holdout-order_Data-embeds_Model-BayesC_Seed-1_mu.dat"
mu1 <- as.mcmc(scan(m1))
m2 <- "results/embedsAll10kIts/models/fit_Holdout-order_Data-embeds_Model-BayesC_Seed-2_mu.dat"
mu2 <- as.mcmc(scan(m2))
m3 <- "results/embedsAll10kIts/models/fit_Holdout-order_Data-embeds_Model-BayesC_Seed-3_mu.dat"
mu3 <- as.mcmc(scan(m3))
m4 <- "results/embedsAll10kIts/models/fit_Holdout-order_Data-embeds_Model-BayesC_Seed-4_mu.dat"
mu4 <- as.mcmc(scan(m4))

chains <- mcmc.list(mu1, mu2, mu3,mu4) 
plot(chains)
gelman.diag(chains)
gelman.plot(chains)
