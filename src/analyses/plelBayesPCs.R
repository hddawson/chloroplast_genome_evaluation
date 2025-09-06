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
hold      <- df$Order %in% c("Fabales", "Poales")
yNA       <- y
yNA[hold] <- NA

data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, 0)
invariant_cols <- names(data)[vapply(data, function(x) all(x == x[1L]), logical(1))]
data[, (invariant_cols) := NULL]
X <- as.matrix(data[,-1])

pca <- readRDS("data/tmp/onehot_aln_pca.rds")
pca_scores <- pca$x

# --- testing config

# y <- y[1:500]
# yNA <- yNA[1:500]
# X <- X[1:500,]
# pca_scores <- pca_scores[1:500,]

# ─── calculate kinship matrix ────────────────────────────────────────────────
cat("Computing kinship matrix K...\n")
X_centered <- scale(X, center=TRUE, scale=FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

# ─── define data configurations ──────────────────────────────────────────────
data_configs <- list(
  list(name="markers_only", X=X),
  list(name="pcs_only", X=pca_scores[,1:20]),  # First 20 PCs
  list(name="markers_pcs", X=cbind(X, pca_scores[,1:20]))
)

# ─── define model configurations ─────────────────────────────────────────────
model_configs <- list(
  # Marker-based models
  list(name="BayesA", type="marker", model="BayesA"),
  list(name="BayesB", type="marker", model="BayesB", probIn=0.05),
  list(name="BayesC", type="marker", model="BayesC", probIn=0.05),
  list(name="BL", type="marker", model="BL"),
  
  # RKHS model
  list(name="RKHS", type="kinship", model="RKHS"),
  
  # Combined models (marker + RKHS)
  list(name="BayesA_RKHS", type="combined", marker_model="BayesA", kinship_model="RKHS"),
  list(name="BL_RKHS", type="combined", marker_model="BL", kinship_model="RKHS")
)

seeds <- c(1,2,3,4)

# ─── create grid ─────────────────────────────────────────────────────────────
# For marker and combined models, test all data configs
# For kinship-only models, only use first data config (markers not used anyway)
grid_list <- list()
for(i in seq_along(model_configs)) {
  model <- model_configs[[i]]
  
  if(model$type == "kinship") {
    # RKHS only - data config doesn't matter much, just use first
    for(seed in seeds) {
      grid_list <- append(grid_list, list(data.frame(
        model_idx = i,
        data_idx = 1,
        seed = seed,
        stringsAsFactors = FALSE
      )))
    }
  } else {
    # Marker or combined models - test all data configs
    for(j in seq_along(data_configs)) {
      for(seed in seeds) {
        grid_list <- append(grid_list, list(data.frame(
          model_idx = i,
          data_idx = j,
          seed = seed,
          stringsAsFactors = FALSE
        )))
      }
    }
  }
}

grid <- do.call(rbind, grid_list)
grid$task <- seq_len(nrow(grid))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  model_config <- model_configs[[row$model_idx]]
  data_config <- data_configs[[row$data_idx]]
  seed <- row$seed
  set.seed(seed)
  
  combo_name <- paste0("Data-", data_config$name, "_Model-", model_config$name, "_Seed-", seed)
  cat("\n⏳ Running:", combo_name)
  
  # Build ETA list based on model type
  eta_list <- list()
  
  if (model_config$type == "marker") {
    # Pure marker model
    args <- list(X = data_config$X, model = model_config$model)
    if (!is.null(model_config$probIn)) args$probIn <- model_config$probIn
    eta_list <- list(args)
    
  } else if (model_config$type == "kinship") {
    # Pure kinship model
    eta_list <- list(list(K = K, model = "RKHS"))
    
  } else if (model_config$type == "combined") {
    # Combined marker + kinship
    marker_args <- list(X = data_config$X, model = model_config$marker_model)
    kinship_args <- list(K = K, model = model_config$kinship_model)
    eta_list <- list(marker_args, kinship_args)
  }
  
  fit <- BGLR(
    y       = yNA,
    ETA     = eta_list,
    nIter   = 50000,
    burnIn  = 10000,
    verbose = FALSE,
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
    data_type = data_config$name,
    model     = model_config$name,
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


res <- read.csv("results/pcs/poales_fabales_grid_results_summary.csv")
boxplot(Spearman ~ model,data=res)
boxplot(Spearman ~ model + data_type,data=res,las=2)
plot(res$R2,res$MAE)
plot(res$R2,res$Spearman)

library(ggplot2)

res[res$model=="RKHS","data_type"] <- "K"

ggplot(res,
       aes(x=data_type, y=Spearman, fill=model)) +
  geom_boxplot(outlier.size=0.5) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Data type",
    y = "Spearman correlation",
    fill = "Model"
  )


mu1 <-  "results/pcs/models/fit_Data-binary_Mk-BL_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE))
mu2 <- "results/pcs/models/fit_Data-binary_Mk-BL_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE))
mu3 <-  "results/pcs/models/fit_Data-binary_Mk-BL_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE))
mu4 <-  "results/pcs/models/fit_Data-binary_Mk-BL_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE))
mu5 <-  "results/pcs/models/fit_Data-binary_Mk-BL_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE))

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BL Topt majminor")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt majMinor")

library(coda)



# Iterate over all unique model/data_type combos
for (combo_base in unique(sub("_Seed-[0-9]+", "", res$combo))) {
  
  # Get seeds for this combo
  seeds <- res$seed[grep(combo_base, res$combo)]
  
  # Read all chains
  chains_list <- list()
  for (s in seeds) {
    mu_file <- paste0("results/pcs/models/fit_", combo_base, "_Seed-", s,"_mu.dat")
    if (file.exists(mu_file)) {
      cat("Reading:", mu_file, "\n")
      chains_list[[length(chains_list)+1]] <- as.mcmc(scan(mu_file, quiet=TRUE))
    } else {
      warning("File not found: ", mu_file)
    }
  }
  
  if (length(chains_list) > 1) {
    chains <- mcmc.list(chains_list)
    
    # Extract nice labels for plotting
    parts <- strsplit(combo_base, "_")[[1]]
    dt_label <- gsub("Data-", "", parts[1])
    model_label <- gsub("Model-", "", parts[2])
    title_str <- paste(model_label, "|", dt_label)
    
    # Trace + density
    plot(chains, main=combo_base)
    
    # Gelman diagnostics
    print(gelman.diag(chains))
    gelman.plot(chains, main=combo_base)
  }
}

res %>%
  group_by(model, data_type) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_Spearman = mean(Spearman, na.rm = TRUE)
  )







