library(ape)
library(adephylo)
library(phylobase)
library(igraph)
library(arrow)
library(data.table)
data <- read_parquet("data/processed_data.parquet")
tree <- read.tree("data/Fasttree/Fasttree.nwk")

aln_data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))
for (j in seq_along(aln_data)) set(aln_data, which(is.na(aln_data[[j]])), j, 0)
invariant_cols <- names(aln_data)[vapply(aln_data, function(x) all(x == x[1L]), logical(1))]
aln_data[, (invariant_cols) := NULL]
aln_mat <- as.matrix(aln_data[,-1])
dim(aln_mat)

#pca <- prcomp(aln_mat, center = T, scale. = T, rank.=10)
#saveRDS(pca, "data/tmp/pca.rds")
pca <- readRDS("data/tmp/pca.rds")
barplot(pca$sdev^2 / sum(pca$sdev^2))
plot(pca$x[,5],pca$x[,6])
cat("pcadone!")
scores <- pca$x[,1:10]
#scores <- as.data.frame(matrix(rnorm(10857*5), nrow=10857, ncol=5))
#rownames(scores) <- data$ID

poales_data <- data[grep("Spermatophyt", data$Taxonomy),]
scores <- scores[grep("Spermatophyt", data$Taxonomy),]
poales_tree <- keep.tip(tree, poales_data$ID)
poales_data <- poales_data[match(poales_tree$tip.label, poales_data$ID), ]
scores <- scores[match(poales_tree$tip.label, poales_data$ID),]
complete_cases <- complete.cases(poales_data$pheno_Topt_site_p50)
poales_data <- poales_data[complete_cases, ]
poales_tree <- keep.tip(poales_tree, poales_data$ID)

tree <- poales_tree

g <- graph_from_edgelist(tree$edge, directed = FALSE)
topo.dist <- distances(g, v = 1:length(tree$tip.label), to = 1:length(tree$tip.label))
classes <- 1:max(topo.dist)

moran.results <- vector("list", length(classes))
names(moran.results) <- paste0("class_", classes)

for(k in classes){
  W <- matrix(0, nrow = nrow(topo.dist), ncol = ncol(topo.dist),
              dimnames = dimnames(topo.dist))
  idx <- which(topo.dist == k, arr.ind = TRUE)
  if(nrow(idx) > 0){
    for(r in seq_len(nrow(idx))){
      W[idx[r,1], idx[r,2]] <- 1
      W[idx[r,2], idx[r,1]] <- 1
    }
  }
  if(all(W == 0)){ moran.results[[k]] <- NA; next }
  sum(is.na(scores))
  moran.results[[k]] <- abouheif.moran(scores, W=W, nrepet=99, alter="two-sided")
  cat(sprintf("%s: %f\n", k, moran.results[[k]]$obs[1]))
}

saveRDS(moran.results, "results/multi_pc_moran_results.rds")

moran.results <- readRDS("results/multi_pc_moran_results.rds")
# extract values
I <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$obs[1]
})

pval <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$pvalue[1]
})

plot(I)
