#!/usr/bin/env Rscript
# bayesian_embeddings_gridsearch.R

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(BGLR)
  library(parallel)
  library(arrow)
  library(Matrix)
  library(stringr)
})

# ─── parse args ──────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d","--data-dir"),   type="character", default="data",
              help="where processed_data.parquet lives"),
  make_option(c("-o","--output-dir"), type="character", default="results",
              help="where to write all outputs"),
  make_option(c("-c","--cores"),      type="integer",  default=5,
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

# ─── load and prepare data ───────────────────────────────────────────────────
cat("Loading processed data...\n")
processed_data <- as.data.table(read_parquet(file.path(data_dir, "processed_data.parquet")))

cat("Loading embeddings...\n")
embed_all <- fread(file.path(data_dir, "combined_embeddings.csv"))
embed_all$LUI <- embed_all$V1
embed_all <- embed_all[,-1]
embed_all[, Species := sub(".*[0-9]{7}(.*)$", "\\1", LUI)]
embed_all$Species <- str_split_i(embed_all$Species,"_",1)

cat("Merging data...\n")
processed_data_subset <- merge(
  processed_data,
  embed_all,
  by.x = "mergeSpecies",
  by.y = "Species",
  all.x = FALSE,
  sort = FALSE
)

embedding_cols <- setdiff(names(processed_data_subset), names(processed_data))  
embedding_cols <- embedding_cols[sapply(processed_data_subset[, ..embedding_cols], is.numeric)]
Z <- as.matrix(processed_data_subset[, ..embedding_cols])

cat("Loading alignment data...\n")
aln <- as.data.table(read_parquet(file.path(data_dir, "tmp/majMinor_aln.pq")))
aln_subset <- aln[match(processed_data_subset$ID, aln[[1]]), ]
for (j in seq_along(aln_subset)) set(aln_subset, which(is.na(aln_subset[[j]])), j, 0)
invariant_cols <- names(aln_subset)[vapply(aln_subset, function(x) all(x == x[1L]), logical(1))]
aln_subset[, (invariant_cols) := NULL]
X <- as.matrix(aln_subset[,-1])

cat("Computing kinship matrix...\n")
X_centered <- scale(X, center = TRUE, scale = FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)

y <- processed_data_subset$pheno_Topt_site_p50 * 0.01

# ─── data validation ─────────────────────────────────────────────────────────
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))
stopifnot(nrow(Z) == length(y))
cat("Data dimensions: Z =", dim(Z), ", K =", dim(K), ", y length =", length(y), "\n")

# ─── define grid ─────────────────────────────────────────────────────────────
seeds <- c(1,2,3,4,5)
grid <- data.frame(seed = seeds, task = seq_along(seeds))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  seed <- row$seed
  set.seed(seed)
  
  combo_name <- paste0("Embeddings_Kinship_Seed-", seed)
  cat("\n⏳ Running:", combo_name, "\n")
  
  ETA <- list(
    Embeddings = list(X = Z, model = "BayesC"),
    Kinship = list(K = K, model = "RKHS")
  )
  
  fit <- BGLR(
    y       = y,
    ETA     = ETA,
    nIter   = 100000,
    burnIn  = 25000,
    verbose = TRUE,
    saveAt  = file.path(output_dir, "models", paste0("fit_", combo_name, "_"))
  )
  
  saveRDS(fit, file.path(output_dir, "models", paste0("fit_", combo_name, ".rds")))
  
  preds <- fit$yHat
  truth <- y
  
  write.csv(
    data.frame(ID = processed_data_subset$ID, y_true = truth, y_pred = preds),
    file.path(output_dir, "predictions", paste0(combo_name, "_preds.csv")),
    row.names = FALSE
  )
  
  metrics <- data.frame(
    combo     = combo_name,
    seed      = seed,
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

write.csv(summary_df, file.path(output_dir, "embeddings_grid_results_summary.csv"), row.names = FALSE)
cat("\n Grid search complete! Summary at embeddings_grid_results_summary.csv\n")


library(coda)
mu1 <-  "results/embeds/models/fit_Embeddings_Kinship_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE)[1:17000])
mu2 <- "results/embeds/models/fit_Embeddings_Kinship_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE)[1:17000])
mu3 <-  "results/embeds/models/fit_Embeddings_Kinship_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE)[1:17000])
mu4 <-  "results/embeds/models/fit_Embeddings_Kinship_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE)[1:17000])
mu5 <-  "results/embeds/models/fit_Embeddings_Kinship_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE)[1:17000])

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt Embeddings + Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt Embeddings + Kinship")
