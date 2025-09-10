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
stopifnot(nrow(X) == length(y))
cat("Data dimensions: Z =", dim(Z), ", X =", dim(X), ", K =", dim(K), ", y length =", length(y), "\n")

# ─── define grid ─────────────────────────────────────────────────────────────
seeds <- c(1,2,3,4,5)
data_types <- c("Embeddings", "Markers")
models <- c("BayesC", "BL")

grid <- expand.grid(
  seed = seeds,
  data_type = data_types,
  model = models,
  stringsAsFactors = FALSE
)
grid$task <- seq_len(nrow(grid))

# ─── parallel task runner ────────────────────────────────────────────────────
run_task <- function(row_index) {
  row <- grid[row_index,]
  seed <- row$seed
  data_type <- row$data_type
  model <- row$model
  
  set.seed(seed)
  
  combo_name <- paste0(data_type, "_", model, "_Seed-", seed)
  cat("\n⏳ Running:", combo_name, "\n")
  
  # Select data matrix based on data_type
  if (data_type == "Embeddings") {
    data_matrix <- Z
  } else if (data_type == "Markers") {
    data_matrix <- X
  }
  
  ETA <- list(
    Embeddings = list(X = data_matrix, model = model),
    Kinship = list(K = K, model = "RKHS")
  )
  
  fit <- BGLR(
    y       = y,
    ETA     = ETA,
    nIter   = 100000,
    burnIn  = 25000,
    verbose = FALSE,
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
    data_type = data_type,
    model     = model,
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

write.csv(summary_df, file.path(output_dir, "data_comparison_grid_results.csv"), row.names = FALSE)

# ─── summary statistics ──────────────────────────────────────────────────────
cat("\n=== Summary Statistics ===\n")
summary_stats <- summary_df[, .(
  mean_RMSE = mean(RMSE),
  sd_RMSE = sd(RMSE),
  mean_R2 = mean(R2),
  sd_R2 = sd(R2)
), by = .(data_type, model)]

print(summary_stats)

cat("\n Grid search complete! Summary at data_comparison_grid_results.csv\n")


#hypothesis here - the embeddings explain more additional variance than the kinship 
#so markers + kinship R^2 < embeddings + kinship
#test BL as a baseline,

res <- read.csv("results/alnVsEmbeds/data_comparison_grid_results.csv")
par(mfrow=c(1,1))
boxplot(R2 ~ data_type, data=res, main="TOpt ~ X + K + e; Model fit comparison",
        xlab="X")
boxplot(Spearman ~ data_type, data=res)
boxplot(MAE ~ data_type, data=res)

#but there is an odd stripe in some of these marker plots, they are jsut guessing the mean

preds <- read.csv("results/alnVsEmbeds/predictions/Markers_BL_Seed-3_preds.csv")
hist(preds$y_pred)
median(preds$y_pred)
sum(preds$y_pred==median(preds$y_pred))
range(preds[grep("NC_", preds$ID),"y_pred"])

length(grep("NC_", processed_data_subset$ID))

nc_indices <- grep("NC_", processed_data_subset$ID)
other_indices <- setdiff(1:nrow(processed_data_subset), nc_indices)

cat("NC_ samples:", length(nc_indices), "\n")
cat("Other samples:", length(other_indices), "\n")

# Check marker variance
X_nc <- X[nc_indices, ]
X_other <- X[other_indices, ]

cat("Marker variance for NC_:", mean(apply(X_nc, 2, var)), "\n")
cat("Marker variance for others:", mean(apply(X_other, 2, var)), "\n")

# Check if NC_ samples are identical/near-identical
nc_distances <- dist(X_nc)
cat("Max distance between NC_ samples:", max(nc_distances), "\n")
cat("Min distance between NC_ samples:", min(nc_distances), "\n")

dim(X_nc)
sum(X_nc)

# Check embedding variance comparison
Z_nc <- Z[nc_indices, ]
Z_other <- Z[other_indices, ]



cat("Embedding variance for NC_:", mean(apply(Z_nc, 2, var)), "\n")
cat("Embedding variance for others:", mean(apply(Z_other, 2, var)), "\n")

# Look at the alignment subsetting step
cat("Original alignment dimensions:", dim(aln), "\n")
cat("Subset alignment dimensions:", dim(aln_subset), "\n")
cat("Invariant columns removed:", length(invariant_cols), "\n")

# Check if NC_ samples were properly matched
nc_in_original <- grep("NC", aln[[1]])
cat("NC_ samples in original alignment:", length(nc_in_original), "\n")

# Check matching success
match_results <- match(processed_data_subset$ID, aln[[1]])
nc_match_success <- sum(!is.na(match_results[nc_indices]))
cat("NC_ samples successfully matched:", nc_match_success, "\n")






library(coda)
mu1 <-  "results/alnVsEmbeds/models/fit_Embeddings_BayesC_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE)[1:18000])
mu2 <- "results/alnVsEmbeds/models/fit_Embeddings_BayesC_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE)[1:18000])
mu3 <-  "results/alnVsEmbeds/models/fit_Embeddings_BayesC_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE)[1:18000])
mu4 <-  "results/alnVsEmbeds/models/fit_Embeddings_BayesC_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE)[1:18000])
mu5 <-  "results/alnVsEmbeds/models/fit_Embeddings_BayesC_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE)[1:18000])

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt Embeddings + Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt Embeddings + Kinship")

mu1 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE)[1:10000])
mu2 <- "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE)[1:10000])
mu3 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE)[1:10000])
mu4 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE)[1:10000])
mu5 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE)[1:10000])

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BL Topt Embeddings + Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BL Topt Embeddings + Kinship")

mu1 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE)[1:10000])
mu2 <- "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE)[1:10000])
mu3 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE)[1:10000])
mu4 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE)[1:10000])
mu5 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE)[1:10000])

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt Markers + Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt Markers + Kinship")

ETA_Kinship_var1 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-1_ETA_Kinship_varU.dat"
m1 <- as.mcmc(scan(ETA_Kinship_var1, quiet = TRUE)[1:10000])
ETA_Kinship_var2 <- "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-2_ETA_Kinship_varU.dat"
m2 <- as.mcmc(scan(ETA_Kinship_var2, quiet = TRUE)[1:10000])
ETA_Kinship_var3 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-3_ETA_Kinship_varU.dat"
m3 <- as.mcmc(scan(ETA_Kinship_var3, quiet = TRUE)[1:10000])
ETA_Kinship_var4 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-4_ETA_Kinship_varU.dat"
m4 <- as.mcmc(scan(ETA_Kinship_var4, quiet = TRUE)[1:10000])
ETA_Kinship_var5 <-  "results/alnVsEmbeds/models/fit_Markers_BayesC_Seed-5_ETA_Kinship_varU.dat"
m5 <- as.mcmc(scan(ETA_Kinship_var5, quiet = TRUE)[1:10000])

chains <- mcmc.list(m1, m2, m3,m4,m5) 
plot(chains, main="BayesC Topt Markers + Kinship")
gelman.diag(chains)
gelman.plot(chains, main="BayesC Topt Markers + Kinship")

mu1 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-1_mu.dat"
m1 <- as.mcmc(scan(mu1, quiet = TRUE)[1:10000])
mu2 <- "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-2_mu.dat"
m2 <- as.mcmc(scan(mu2, quiet = TRUE)[1:10000])
mu3 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-3_mu.dat"
m3 <- as.mcmc(scan(mu3, quiet = TRUE)[1:10000])
mu4 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-4_mu.dat"
m4 <- as.mcmc(scan(mu4, quiet = TRUE)[1:10000])
mu5 <-  "results/alnVsEmbeds/models/fit_Embeddings_BL_Seed-5_mu.dat"
m5 <- as.mcmc(scan(mu5, quiet = TRUE)[1:10000])




