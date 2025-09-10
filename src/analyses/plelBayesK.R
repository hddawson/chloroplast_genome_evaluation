#!/usr/bin/env Rscript
# bayes_alphabet_gridsearch_parallel.R
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(BGLR)
  library(parallel)
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
df        <- as.data.table(read_parquet(file.path(data_dir, "processed_data.parquet")))
y         <- df$pheno_Topt_site_p50
data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X <- as.matrix(data[,-1])

# ─── testing subsetting
#y <- y[1:500]
#df <- df[1:500,]
#X <- X[1:500,]

# ─── calculate kinship matrix ────────────────────────────────────────────────
cat("Computing kinship matrix K...\n")
X_centered <- scale(X, center=TRUE, scale=FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

# ─── define data configurations ──────────────────────────────────────────────
data_configs <- list(
  list(name="markers_only", X=X),
  list(name="markers_plus_K", X=X, K=K),
  list(name="K_only", K=K)
)

# ─── define model configurations ─────────────────────────────────────────────
model_configs <- list(
  list(name="BayesC", type="marker", model="BayesC", probIn=0.05),
  list(name="BL", type="marker", model="BL"),
  list(name="BayesC_RKHS", type="combined", marker_model="BayesC", kinship_model="RKHS"),
  list(name="BL_RKHS", type="combined", marker_model="BL", kinship_model="RKHS")
)
seeds <- c(1,2,3,4)

# ─── create grid ─────────────────────────────────────────────────────────────
grid <- expand.grid(
  data_config = sapply(data_configs, function(x) x$name),
  model_config = sapply(model_configs, function(x) x$name),
  seed = seeds,
  stringsAsFactors = FALSE
)

# Filter valid combinations
valid_combinations <- data.frame(
  data_config = c("markers_only", "markers_only", "markers_plus_K", "markers_plus_K", "K_only"),
  model_config = c("BayesC", "BL", "BayesC_RKHS", "BL_RKHS", "BL_RKHS"),
  stringsAsFactors = FALSE
)

grid <- merge(grid, valid_combinations, by = c("data_config", "model_config"))
grid$task <- seq_len(nrow(grid))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  data_name <- row$data_config
  model_name <- row$model_config
  seed <- row$seed
  
  set.seed(seed)
  
  combo_name <- paste0(data_name, "_", model_name, "_seed", seed)
  cat("\n⏳ Running:", combo_name, "\n")
  
  # Get data config
  data_config <- data_configs[[which(sapply(data_configs, function(x) x$name) == data_name)]]
  model_config <- model_configs[[which(sapply(model_configs, function(x) x$name) == model_name)]]
  
  # Build ETA list
  ETA <- list()
  # Handle K_only case
  if (data_name == "K_only") {
    stopifnot(!is.null(data_config$K))
    ETA <- list(kinship = list(K = data_config$K, model = "RKHS"))
  } else if (model_config$type == "marker") {
    stopifnot(!is.null(data_config$X))
    ETA$markers <- list(X = data_config$X, model = model_config$model)
    if (!is.null(model_config$probIn)) {
      ETA$markers$probIn <- model_config$probIn
    }
  } else if (model_config$type == "combined") {
    stopifnot(!is.null(data_config$X) && !is.null(data_config$K))
    ETA$markers <- list(X = data_config$X, model = model_config$marker_model)
    ETA$kinship <- list(K = data_config$K, model = model_config$kinship_model)
  }
  
  
  
  fit <- BGLR(
    y       = y,
    ETA     = ETA,
    nIter   = 40000,
    burnIn  = 8000,
    verbose = FALSE,
    saveAt  = file.path(output_dir, "models", paste0("fit_", combo_name, "_"))
  )
  
  saveRDS(fit, file.path(output_dir, "models", paste0("fit_", combo_name, ".rds")))
  
  preds <- fit$yHat
  truth <- y
  
  write.csv(
    data.frame(ID = seq_along(y), y_true = truth, y_pred = preds),
    file.path(output_dir, "predictions", paste0(combo_name, "_preds.csv")),
    row.names = FALSE
  )
  
  metrics <- data.frame(
    combo        = combo_name,
    data_config  = data_name,
    model_config = model_name,
    seed         = seed,
    RMSE         = rmse(truth, preds),
    MAE          = mae(truth, preds),
    R2           = r2(truth, preds),
    Spearman     = spearman(truth, preds),
    stringsAsFactors = FALSE
  )
  
  png(file.path(output_dir, "plots", paste0(combo_name, "_plot.png")),
      width = 600, height = 600)
  plot(preds, truth,
       pch = 19, cex = 1.2,
       xlab = "Predicted Topt",
       ylab = "Observed Topt",
       main = combo_name)
  abline(a = 0, b = 1, col = "blue", lwd = 2)
  dev.off()
  
  return(metrics)
}

# ─── run in parallel ─────────────────────────────────────────────────────────
cat("\n▶ Launching grid search with", nrow(grid), "models using", n_cores, "cores\n")
results <- mclapply(grid$task, run_task, mc.cores = n_cores)
summary_df <- rbindlist(results)

write.csv(summary_df, file.path(output_dir, "bayes_alphabet_grid_results.csv"), row.names = FALSE)

# ─── summary statistics ──────────────────────────────────────────────────────
cat("\n=== Summary Statistics ===\n")
summary_stats <- summary_df[, .(
  mean_RMSE = mean(RMSE),
  sd_RMSE = sd(RMSE),
  mean_R2 = mean(R2),
  sd_R2 = sd(R2)
), by = .(data_config, model_config)]

print(summary_stats)

cat("\n Grid search complete! Summary at bayes_alphabet_grid_results.csv\n")

res <- read.csv("results/K/bayes_alphabet_grid_results.csv")

ggplot(res, aes(x = interaction(data_config, model_config), y = R2)) +
  geom_boxplot() +
  xlab("data_config : model_config") +
  ylab("R²") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))


K_res <- read_csv("results/K/predictions/K_only_BL_RKHS_seed1_preds.csv")
K_res <- K_res * 0.01
plot(K_res$y_pred, K_res$y_true)
K_res$residuals <- K_res$y_true - K_res$y_pred
plot(K_res$residuals)
plot(K_res$residuals,K_res$y_true)
cor(K_res$residuals,K_res$y_true)
abline(a=0,b=1, col="red")

mu1 <-  "results/K/models/fit_K_only_BL_RKHS_seed1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <-"results/K/models/fit_K_only_BL_RKHS_seed2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/K/models/fit_K_only_BL_RKHS_seed3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/K/models/fit_K_only_BL_RKHS_seed4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
chains <- mcmc.list(m1, m2, m3,m4) 
plot(chains, main="RKHS Topt Kinship")
gelman.diag(chains)
gelman.plot(chains, main="RKHS Topt Markers + Kinship")


mu1 <-  "results/K/models/fit_markers_only_BayesC_seed1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <-"results/K/models/fit_markers_only_BayesC_seed2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/K/models/fit_markers_only_BayesC_seed3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/K/models/fit_markers_only_BayesC_seed4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
chains <- mcmc.list(m1, m2, m3,m4) 
plot(chains, main="RKHS Topt Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt Markers + Kinship")


mu1 <-  "results/K/models/fit_markers_only_BL_seed1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <-"results/K/models/fit_markers_only_BL_seed2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/K/models/fit_markers_only_BL_seed3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/K/models/fit_markers_only_BL_seed4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
chains <- mcmc.list(m1, m2, m3,m4) 
plot(chains, main="markers_only_BL")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt Markers + Kinship")

mu1 <-  "results/K/models/fit_markers_plus_K_BL_RKHS_seed1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <-"results/K/models/fit_markers_plus_K_BL_RKHS_seed2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/K/models/fit_markers_plus_K_BL_RKHS_seed3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/K/models/fit_markers_plus_K_BL_RKHS_seed4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
chains <- mcmc.list(m1, m2, m3,m4) 
plot(chains, main="markers_plus_K_BL_RKHS")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt Markers + Kinship")










