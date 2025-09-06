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

data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X <- as.matrix(data[,-1])

pca <- readRDS("data/tmp/onehot_aln_pca.rds")
pca_scores <- pca$x

# --- testing config
#y <- y[1:500]
#X <- X[1:500,]
#pca_scores <- pca_scores[1:500,]

# ─── calculate kinship matrix ────────────────────────────────────────────────
cat("Computing kinship matrix K...\n")
X_centered <- scale(X, center=TRUE, scale=FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

# ─── define holdout strategies ──────────────────────────────────────────────
holdout_strategies <- list(
  list(name="order", test_idx=which(df$Order %in% c("Fabales", "Poales"))),
  list(name="random", test_idx=NULL)  # Will be set per seed
)

# ─── define data configurations ──────────────────────────────────────────────
data_configs <- list(
  list(name="markers_only", X=X, fixed_effects=NULL),
  list(name="markers_pcfixed", X=X, fixed_effects=pca_scores[,1:20])  # PCs as fixed effects
)

# ─── define model configurations ─────────────────────────────────────────────
model_configs <- list(
  list(name="BayesB", type="marker", model="BayesB", probIn=0.05),
  list(name="BayesC", type="marker", model="BayesC", probIn=0.05),
  list(name="BL", type="marker", model="BL"),
  list(name="RKHS", type="kinship", model="RKHS"),
  list(name="BayesA_RKHS", type="combined", marker_model="BayesA", kinship_model="RKHS"),
  list(name="BL_RKHS", type="combined", marker_model="BL", kinship_model="RKHS")
)

seeds <- c(1,2,3)

# ─── create grid ─────────────────────────────────────────────────────────────
grid_list <- list()
for(h in seq_along(holdout_strategies)) {
  holdout <- holdout_strategies[[h]]
  
  for(i in seq_along(model_configs)) {
    model <- model_configs[[i]]
    
    if(model$type == "kinship") {
      # RKHS only - test relevant data configs
      relevant_configs <- c(1, 4, 5)  # markers_only, markers_pcfixed, pcs_pcfixed
      for(j in relevant_configs) {
        for(seed in seeds) {
          grid_list <- append(grid_list, list(data.frame(
            holdout_idx = h,
            model_idx = i,
            data_idx = j,
            seed = seed,
            stringsAsFactors = FALSE
          )))
        }
      }
    } else {
      # Marker or combined models - test all data configs
      for(j in seq_along(data_configs)) {
        for(seed in seeds) {
          grid_list <- append(grid_list, list(data.frame(
            holdout_idx = h,
            model_idx = i,
            data_idx = j,
            seed = seed,
            stringsAsFactors = FALSE
          )))
        }
      }
    }
  }
}

grid <- do.call(rbind, grid_list)
grid$task <- seq_len(nrow(grid))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  holdout_config <- holdout_strategies[[row$holdout_idx]]
  model_config <- model_configs[[row$model_idx]]
  data_config <- data_configs[[row$data_idx]]
  seed <- row$seed
  set.seed(seed)
  
  # Set holdout indices
  if(holdout_config$name == "random") {
    n_test <- round(0.1 * length(y))
    test_idx <- sample(length(y), n_test)
  } else {
    test_idx <- holdout_config$test_idx
  }
  
  # Create masked y
  yNA <- y
  yNA[test_idx] <- NA
  
  combo_name <- paste0("Holdout-", holdout_config$name, 
                       "_Data-", data_config$name, 
                       "_Model-", model_config$name, 
                       "_Seed-", seed)
  cat("\n⏳ Running:", combo_name)
  
  # Build ETA list based on model type
  eta_list <- list()
  
  # Add fixed effects if specified
  if(!is.null(data_config$fixed_effects)) {
    eta_list <- append(eta_list, list(list(X=data_config$fixed_effects, model="FIXED")))
  }
  
  if (model_config$type == "marker") {
    # Pure marker model
    args <- list(X = data_config$X, model = model_config$model)
    if (!is.null(model_config$probIn)) args$probIn <- model_config$probIn
    eta_list <- append(eta_list, list(args))
    
  } else if (model_config$type == "kinship") {
    # Pure kinship model
    eta_list <- append(eta_list, list(list(K = K, model = "RKHS")))
    
  } else if (model_config$type == "combined") {
    # Combined marker + kinship
    marker_args <- list(X = data_config$X, model = model_config$marker_model)
    kinship_args <- list(K = K, model = model_config$kinship_model)
    eta_list <- append(eta_list, list(marker_args, kinship_args))
  }
  
  # Ensure we have at least one component
  #stopifnot(length(eta_list) > 0, "ETA list cannot be empty")
  
  fit <- BGLR(
    y       = yNA,
    ETA     = eta_list,
    nIter   = 50000,
    burnIn  = 10000,
    verbose = FALSE,
    saveAt  = file.path(output_dir, "models", paste0("fit_", combo_name, "_"))
  )
  
  saveRDS(fit, file.path(output_dir, "models", paste0("fit_", combo_name, ".rds")))
  
  preds <- fit$yHat[test_idx]
  truth <- y[test_idx]
  
  write.csv(
    data.frame(ID = df$ID[test_idx], y_true = truth, y_pred = preds),
    file.path(output_dir, "predictions", paste0(combo_name, "_preds.csv")),
    row.names = FALSE
  )
  
  metrics <- data.frame(
    combo     = combo_name,
    seed      = seed,
    holdout   = holdout_config$name,
    data_type = data_config$name,
    model     = model_config$name,
    n_test    = length(test_idx),
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
       xlab = paste("Predicted (", holdout_config$name, " holdout)"),
       ylab = "Observed",
       main = paste(data_config$name, model_config$name, sep=" + "))
  abline(a = 0, b = 1, col = "blue", lwd = 2)
  dev.off()
  
  return(metrics)
}

# ─── run in parallel ─────────────────────────────────────────────────────────
cat("\n▶ Launching grid search with", nrow(grid), "models using", n_cores, "cores\n")
results <- mclapply(grid$task, run_task, mc.cores = n_cores)
summary_df <- rbindlist(results)

write.csv(summary_df, file.path(output_dir, "popstruct_grid_results_summary.csv"), row.names = FALSE)
cat("\n✅ Grid search complete! Summary at popstruct_grid_results_summary.csv\n")


res <- read.csv("results/popstruct_grid_results_summary.csv")
