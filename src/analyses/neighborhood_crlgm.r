library(ape)
library(adephylo)
library(phylobase)
library(igraph)
library(arrow)

data <- read_parquet("data/processed_data.parquet")
tree <- read.tree("data/Fasttree/Fasttree.nwk")

aln_data <- read_parquet("data/tmp/cds_supermatrix.parquet")
for (j in seq_along(aln_data)) set(aln_data, which(is.na(aln_data[[j]])), j, 0)
invariant_cols <- names(aln_data)[vapply(aln_data, function(x) all(x == x[1L]), logical(1))]
aln_data[, (invariant_cols) := NULL]
aln_mat <- as.matrix(aln_data[,-1])
dim(aln_mat)
pca <- prcomp(aln_mat, center = T, scale. = T, rank.=5)
saveRDS(pca, "data/tmp/pca.rds")
pca < - readRDS("data/tmp/pca.rds")
scores <- pca$x[,1;5]

poales_data <- data[grep(
  "Poales", data$Taxonomy
),]

# Match tree to data
poales_tree <- keep.tip(tree, poales_data$ID)

# Ensure data order matches tree
poales_data <- poales_data[match(poales_tree$tip.label, poales_data$ID), ]

# Remove NAs
complete_cases <- complete.cases(poales_data$pheno_Topt_site_p50)
poales_data <- poales_data[complete_cases, ]
poales_tree <- keep.tip(poales_tree, poales_data$ID)

poales_tree_fixed <- poales_tree

tree <- poales_tree_fixed
trait <- poales_data$pheno_Topt_site_p50
stopifnot(all(names(trait) == tree$tip.label))

# 1. compute node distances (topological, not patristic)
# dist.nodes gives distances among all nodes; subset to tips
g <- graph_from_edgelist(tree$edge, directed = FALSE)
topo.dist <- distances(g, v = 1:length(tree$tip.label), 
                       to = 1:length(tree$tip.label))
hist(topo.dist, main="Spermatophyte topological distance distribution |Fasttree|n=10857|")
sum(topo.dist == 4)

# 2. define classes by integer steps of node distance
max.dist <- max(topo.dist)
classes <- 1:max.dist

moran.results <- vector("list", length(classes))
names(moran.results) <- paste0("class_", classes)

for(k in classes){
  cat("Processing neighborhood class", k, "\n")
  W <- matrix(0, nrow = nrow(topo.dist), ncol = ncol(topo.dist),
              dimnames = dimnames(topo.dist))
  idx <- which(topo.dist == k, arr.ind = TRUE)
  if(nrow(idx) > 0){
    for(r in seq_len(nrow(idx))){
      W[idx[r,1], idx[r,2]] <- 1
      W[idx[r,2], idx[r,1]] <- 1
    }
  }
  if(all(W == 0)){
    moran.results[[k]] <- NA
    next
  }
  df <- as.data.frame(trait)
  rownames(df) <- names(trait)
  moran.results[[k]] <- abouheif.moran(df, W = W, nrepet = 99, alter = "two-sided")
}

saveRDS(moran.results, file = "results/spermatophyta_neighborhood_moran_results.rds")
moran.results <- readRDS("results/spermatophyta_neighborhood_moran_results.rds")
# extract values
I <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$obs[1]
})

pval <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$pvalue[1]
})

png("plots/spermatophyta_neighborhood_correlogram.png",
    width = 40, height = 10, units = "in", res = 300)
plot(classes, I, type = "b", pch = 16,
     main="Spermatophyta phylogenetic signal |T_opt~|topological distance|Fasttree|",
     cex.main = 0.9,font.main = 2,
     xlab = "topological distance (nodes apart)", ylab = "Moran's I (Abouheif)")
abline(h = 0, lty = 2)
points(classes, I, col = ifelse(pval < 0.05, "red", "black"))
legend("topright", legend = c("p < 0.05"), pch = 16, col = "red", bty = "n")
dev.off()

classSizes <- table(factor(topo.dist, levels=classes))
