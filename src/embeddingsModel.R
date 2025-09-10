library(data.table)
library(arrow)
library(Matrix)
library(BGLR)
library(stringr)


# --- read data
processed_data <-as.data.table(read_parquet("data/processed_data.parquet"))



#let's go get them; 

#write all the genomes to a fasta file 
fasta_file <- "data/selectedGenomes.fasta"
writeLines(
  paste0(">", data$ID, "\n", processed_data$Genome),
  fasta_file
)


#example embedding 
embed <- fread("data/genome_embeddings_csv/NC0015681Epifagusvirginiana_embeddings.csv",
               header=TRUE)


# Get all CSV files in the folder
embed_all <- fread("data/combined_embeddings.csv")
embed_all$LUI <- embed_all$V1
embed_all <- embed_all[,-1]

embed_all[, Species := sub(".*[0-9]{7}(.*)$", "\\1", LUI)]
embed_all$Species <- str_split_i(embed_all$Species,"_",1)

processed_data_subset <- merge(
  processed_data,
  embed_all,
  by.x = "mergeSpecies",
  by.y = "Species",
  all.x = FALSE,   # only keep rows with embeddings
  sort = FALSE
)

embedding_cols <- setdiff(names(processed_data_subset),
                          names(processed_data))  
embedding_cols <- embedding_cols[sapply(processed_data_subset[, ..embedding_cols], is.numeric)]

# Build embedding matrix
Z <- as.matrix(processed_data_subset[, ..embedding_cols])

aln <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
aln_subset <- aln[match(processed_data_subset$ID, aln[[1]]), ]
for (j in seq_along(aln_subset)) set(aln_subset, which(is.na(aln_subset[[j]])), j, 0)
invariant_cols <- names(aln_subset)[vapply(aln_subset, function(x) all(x == x[1L]), logical(1))]
aln_subset[, (invariant_cols) := NULL]

X <- as.matrix(aln_subset[,-1])
X_centered <- scale(X, center = TRUE, scale = FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)

y <- processed_data_subset$pheno_Topt_site_p50 * 0.01

summary(as.vector(K))
hist(as.vector(K))

stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

dim(Z)
dim(K)
length(y)

ETA <- list(
  Embeddings = list(X = Z, model = "BayesC"),  # fixed effect from embeddings
  Kinship = list(K = K, model = "RKHS")        # random effect via kinship
)

seeds <- c(1,2,3)

for (seed in seeds) {
  fit <- BGLR(
    y = y,
    ETA = ETA,
    nIter = 100000,
    burnIn = 25000,
    verbose = TRUE,
    saveAt = "results/embeds100k/"
  )
  saveRDS(fit, "results/embeds100k/fit.RDS")
}

quit()
# --- Step 5: Extract results

#evaluate convergence



# Predicted values
y_hat <- fit$yHat

# Posterior means for embedding effects
beta_hat <- fit$ETA$Embeddings$b

# Variance components
var_components <- fit$varE


# --- construct K 

aln <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(aln)) set(aln, which(is.na(aln[[j]])), j, 0)
invariant_cols <- names(aln)[vapply(aln, function(x) all(x == x[1L]), logical(1))]
aln[, (invariant_cols) := NULL]
X <- as.matrix(aln[,-1])

cat("Computing kinship matrix K...\n")
X_centered <- scale(X, center=TRUE, scale=FALSE)
K <- tcrossprod(X_centered) / ncol(X_centered)
stopifnot(nrow(K) == length(y))
stopifnot(ncol(K) == length(y))

#merge K on processed_data
y <- processed_data$pheno_Topt_site_p50


pca <- prcomp(Z, rank.=4)

pvar <- pca$sdev^2 * 100 / sum(pca$sdev^2) 

par(mfrow=c(1,2))
plot(pvar[1:20], main="PCA on 6061 embeddings")
plot(pca$x[,2],pca$x[,3])
points(pca$x[grep("Poales", processed_data_subset$Taxonomy),2],
       pca$x[grep("Poales", processed_data_subset$Taxonomy),3],
       col="seagreen")
points(pca$x[grep("Asparagales", processed_data_subset$Taxonomy),2],
       pca$x[grep("Asparagales", processed_data_subset$Taxonomy),3],
       col="purple")

points(pca$x[grep("Fabales", processed_data_subset$Taxonomy),2],
       pca$x[grep("Fabales", processed_data_subset$Taxonomy),3],
       col="yellow")

plot(pca$x[,4],pca$x[,3])
points(pca$x[grep("Poales", processed_data_subset$Taxonomy),4],
       pca$x[grep("Poales", processed_data_subset$Taxonomy),3],
       col="seagreen")
points(pca$x[grep("Asparagales", processed_data_subset$Taxonomy),4],
       pca$x[grep("Asparagales", processed_data_subset$Taxonomy),3],
       col="purple")

points(pca$x[grep("Fabales", processed_data_subset$Taxonomy),4],
       pca$x[grep("Fabales", processed_data_subset$Taxonomy),3],
       col="yellow")

points(pca$x[grep("Pinales", processed_data_subset$Taxonomy),4],
       pca$x[grep("Pinales", processed_data_subset$Taxonomy),3],
       col="blue")

points(pca$x[grep("Brassicales", processed_data_subset$Taxonomy),4],
       pca$x[grep("Brassicales", processed_data_subset$Taxonomy),3],
       col="red")


mod <- lm(processed_data_subset$pheno_Topt_site_p50 * 0.01 ~
            pca$x[,1] + pca$x[,2] + pca$x[,3] + pca$x[,4])
summary(mod)  

plot(processed_data_subset$pheno_Topt_site_p50 * 0.01,mod$fitted.values)

mod <- lm(pheno_Topt_site_p50~
            M + T + A + I + L + Y + E + R + S + P + F + V + W + K + N + C + D + G + H + Q, data=processed_data_subset)
summary(mod)  

plot(processed_data_subset$pheno_Topt_site_p50 * 0.01,mod$fitted.values)
  