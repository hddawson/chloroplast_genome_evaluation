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
df        <- fread(file.path(data_dir, "alignment_input_data.csv"))
y         <- df$pheno_Topt_site_p50
#order_mat <- model.matrix(~ Order, data=df)[, -1, drop=FALSE]

aln_cols  <- grep("supermatrix", names(df), value=TRUE)
df <- read_parquet("data/tmp/onehot_aln.pq")
X <- as.matrix(setDT(df[,-1]))
hold      <- df$Order %in% c("Fabales", "Poales")
yNA       <- y
yNA[hold] <- NA

df <- read_parquet("data/tmp/majMinor_aln.pq")
X_alt <- as.matrix(setDT(df[,-1]))
# K <- ... # Your kinship matrix for RKHS

# ─── define data types ───────────────────────────────────────────────────────
data_types <- list(
  list(name="binary", X=X),
  list(name="onehot", X=X_alt)  
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
  
  combo_name <- paste0("_Data-", dt$name, "_Mk-", mk$name, "_Seed-", seed)
  cat("\n⏳ Running:", combo_name)
  
  # Build ETA list based on model type
  if (mk$model == "RKHS") {
    # For RKHS, use kinship matrix if available, otherwise compute from X
    if (!is.null(dt$K)) {
      eta_list <- list(list(K = dt$K, model = "RKHS"))
    } else {
      # Compute kinship matrix from X
      X_centered <- scale(dt$X, center=TRUE, scale=FALSE)
      X_centered[is.na(X_centered)] <- 0
      K_computed <- tcrossprod(X_centered) / ncol(X_centered)
      eta_list <- list(list(K = K_computed, model = "RKHS"))
    }
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
    nIter   = 100000,
    burnIn  = 25000,
    verbose = FALSE,
    saveAt  = file.path(output_dir, "models", paste0("fit_", combo_name, "_"))
  )
  
  saveRDS(fit, paste0("data/tmp/",combo_name))
  
  preds <- fit$yHat[hold]
  truth <- y[hold]
  
  write.csv(
    data.frame(LUI = df$LUI[hold], y_true = truth, y_pred = preds),
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
  
  saveRDS(fit, file.path(output_dir, "models", paste0("fit_", combo_name, ".rds")))
  
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