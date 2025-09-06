#!/usr/bin/env Rscript
# bayes_alphabet_gridsearch_parallel.R

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(BGLR)
  library(parallel)
  library(MASS)
  library(arrow)
})

# ─── parse args ──────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d","--data-dir"),   type="character", default="data",
              help="where processed_data.csv lives"),
  make_option(c("-o","--output-dir"), type="character", default="results",
              help="where to write all outputs"),
  make_option(c("-c","--cores"),      type="integer",  default=10,
              help="number of CPU cores to use")
)
opt        <- parse_args(OptionParser(option_list=option_list))
data_dir   <- opt$`data-dir`
output_dir <- opt$`output-dir`
n_cores    <- opt$`cores`
dirs       <- c("predictions","models","plots")
lapply(file.path(output_dir, dirs), dir.create, recursive=TRUE, showWarnings=FALSE)

# ─── helpers ─────────────────────────────────────────────────────────────────
rmse     <- function(a,b) sqrt(mean((a - b)^2))
mae      <- function(a,b) mean(abs(a - b))
r2       <- function(a,b) cor(a,b)^2
spearman <- function(a,b) cor(a,b, method="spearman")

# ─── load data ───────────────────────────────────────────────────────────────
df        <- fread(file.path(data_dir, "alignment_input_data.csv"))
y         <- df$pheno_Topt_site_p50
hold      <- df$Order %in% c("Fabales", "Poales")
yNA       <- y
yNA[hold] <- NA

data <- as.data.table(read_parquet("data/tmp/onehot_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X <- as.matrix(data[,-1])
dim(X)

data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X_alt <- as.matrix(data[,-1])
dim(X_alt) #X_alt <- majMinor 

# K <- ... # Your kinship matrix for RKHS

# ─── define data types ───────────────────────────────────────────────────────
data_types <- list(
  list(name="binary", X=X_alt),
  list(name="onehot", X=X)  
  # list(name="onehot", X=X_alt, K=K)  # Uncomment and modify when you add K
)

# ─── define grid ─────────────────────────────────────────────────────────────
marker_models <- list(
  list(name="BayesA", model="BayesA"),
  list(name="BayesB", model="BayesB", probIn=0.05),
  list(name="BayesC", model="BayesC", probIn=0.05),
  list(name="BL", model="BL")
  #list(name="RKHS", model="RKHS")  # Added RKHS
)
seeds <- c(1,2,3,4,5)

# ─── define tasks ────────────────────────────────────────────────────────────
grid <- expand.grid(
  marker = seq_along(marker_models),
  data_type = seq_along(data_types),
  seed = seeds,
  stringsAsFactors = FALSE
)

grid$task <- seq_len(nrow(grid))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  mk  <- marker_models[[row$marker]]
  dt  <- data_types[[row$data_type]]
  seed <- row$seed
  set.seed(seed)
  
  combo_name <- paste0("Data-", dt$name, "_Mk-", mk$name, "_Seed-", seed)
  cat("\n⏳ Running:", combo_name)
  
  # Build ETA list based on model type
  if (mk$model == "RKHS") {
    # For RKHS, use kinship matrix if available, otherwise compute from X
    eta_list <- list(list(K = dt$K, model = "RKHS"))
  } else {
    # For marker-based models
    eta_list <- list(
      {
        args <- list(X = dt$X, model = mk$model)
        if (!is.null(mk$probIn)) args$probIn <- mk$probIn
        args
      }
    )
  }

  fit <- BGLR(
    y       = yNA,
    ETA     = eta_list,
    nIter   = 50000,
    burnIn  = 10000,
    verbose = TRUE,
    saveAt  = file.path(output_dir, "models", paste0("fit_", combo_name, "_"))
  )
  
  saveRDS(fit, file.path(output_dir, "models", paste0("fit_", combo_name, ".rds")))
  
  preds <- fit$yHat[hold]
  truth <- y[hold]
  
  write.csv(
    data.frame(ID = df$ID[hold], y_true = truth, y_pred = preds),
    file.path(output_dir, "predictions", paste0(combo_name, "_preds.csv")),
    row.names = FALSE
  )
  
  metrics <- data.frame(
    combo     = combo_name,
    seed      = seed,
    data_type = dt$name,
    marker    = mk$name,
    RMSE      = rmse(truth, preds),
    MAE       = mae(truth, preds),
    R2        = r2(truth, preds),
    Spearman  = spearman(truth, preds),
    stringsAsFactors = FALSE
  )
  
  png(file.path(output_dir, "plots", paste0(combo_name, "_plot.png")),
      width = 600, height = 600)
  plot(preds, truth,
       pch = 19, cex = 1.2,
       xlab = "Predicted q50 (fabales+poales)",
       ylab = "Observed q50 (fabales+poales)",
       main = combo_name)
  abline(a = 0, b = 1, col = "blue", lwd = 2)
  dev.off()
  
  return(metrics)
}

# ─── run in parallel ─────────────────────────────────────────────────────────
cat("\n▶ Launching grid search with", nrow(grid), "models using", n_cores, "cores\n")
results <- mclapply(grid$task, run_task, mc.cores = n_cores)
summary_df <- rbindlist(results)

write.csv(summary_df, file.path(output_dir, "poales_fabales_grid_results_summary.csv"), row.names = FALSE)
cat("\n✅ Grid search complete! Summary at poales_fabales_grid_results_summary.csv\n")


res <- read.csv("results/predictions/Data-onehot_Mk-BayesB_Seed-1_preds.csv")
plot(res$y_pred,res$y_true)
cor(res$y_pred,res$y_true)
pred_files <- list.files("results/predictions/", pattern = "\\.csv$", full.names = TRUE)

library(coda)
par(mfrow=c(2,2))
mu1 <-  "results/models/fit_Data-binary_Mk-BayesA_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-binary_Mk-BayesA_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-binary_Mk-BayesA_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-binary_Mk-BayesA_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-binary_Mk-BayesA_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesA Topt majMinor")
gelman.diag(chains)
gelman.plot(chains, main="BayesA Topt majMinor")

mu1 <-  "results/models/fit_Data-onehot_Mk-BayesA_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-onehot_Mk-BayesA_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-onehot_Mk-BayesA_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-onehot_Mk-BayesA_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-onehot_Mk-BayesA_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesA Topt onehot")
gelman.diag(chains)
gelman.plot(chains, main="BayesA Topt onehot")

mu1 <-  "results/models/fit_Data-binary_Mk-BayesB_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-binary_Mk-BayesB_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-binary_Mk-BayesB_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-binary_Mk-BayesB_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-binary_Mk-BayesB_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesB Topt majrMinor")
gelman.diag(chains)
gelman.plot(chains, main="BayesB Topt majMinor")

mu1 <-  "results/models/fit_Data-onehot_Mk-BayesC_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-onehot_Mk-BayesC_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-onehot_Mk-BayesC_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-onehot_Mk-BayesC_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-onehot_Mk-BayesC_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt one hot")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt one hot")

mu1 <-  "results/models/fit_Data-binary_Mk-BayesC_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-binary_Mk-BayesC_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-binary_Mk-BayesC_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-binary_Mk-BayesC_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-binary_Mk-BayesC_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt majminor")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt majMinor")

mu1 <-  "results/models/fit_Data-onehot_Mk-BL_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-onehot_Mk-BL_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-onehot_Mk-BL_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-onehot_Mk-BL_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-onehot_Mk-BL_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BL Topt one hot")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt one hot")

mu1 <-  "results/models/fit_Data-binary_Mk-BL_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/models/fit_Data-binary_Mk-BL_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/models/fit_Data-binary_Mk-BL_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/models/fit_Data-binary_Mk-BL_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/models/fit_Data-binary_Mk-BL_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BL Topt majminor")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt majMinor")

for (file in pred_files) {
  print(file)
  df <- fread(file)
  row <- strsplit(basename(file), "_")
  encoding <-strsplit(row[[1]][1],"-")[[1]][2]
  mk <-strsplit(row[[1]][2],"-")[[1]][2]
  seed <-strsplit(row[[1]][3],"-")[[1]][2]
  corr <- cor(df$y_pred, df$y_true, method="pearson")
  rmse_val <- rmse(df$y_pred, df$y_true)
  
}

res <- read.csv("results/poales_fabales_grid_results_summary.csv")
res$data_type[res$data_type=="onehot"] <- "majMinor"
res$data_type[res$data_type=="binary"] <- "onehot"
par(mfrow=c(1,1))
boxplot(Spearman ~ data_type + marker, data=res,
        main="MajMinor vs Onehot SNP encodings (T_opt, eval=Poales+Fabales)")


site_ids <- sub("^(.*_)([0-9]+)_[ATCG]$", "\\2", colnames(X))

# Split the column indices by site
site_cols <- split(seq_along(site_ids), site_ids)

# Function to count non-zero columns (alleles) per site
allele_counts <- sapply(site_cols, function(cols) {
  site_data <- X[, cols, drop = FALSE]
  allele_per_row <- rowSums(site_data) # only one 1 per row in one-hot
  sum(colSums(site_data) > 0)  # number of alleles present at this site
})

# Tally how many sites have 2, 3, or 4 alleles
table(allele_counts[allele_counts %in% 2:4])
dim(data)




pca <- read






