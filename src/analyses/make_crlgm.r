library(ape)
library(data.table)
library(arrow)
library(dplyr)
library(ape)
library(adephylo)
library(geiger)
library(phylobase)


data <- as.data.frame(read_parquet("data/processed_data.parquet"))
tree <- read.tree("data/Fasttree/Fasttree.nwk")

poales_data <- data[grep(
    "Spermatophyta", data$Taxonomy
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
zero_branches <- poales_tree_fixed$edge.length == 0
very_small <- poales_tree_fixed$edge.length < 1e-6

cat("Fixing", sum(zero_branches), "zero branches and", 
    sum(very_small), "very small branches\n")

# Set minimum branch length
min_branch <- 1e-6
poales_tree_fixed$edge.length[poales_tree_fixed$edge.length < min_branch] <- min_branch
poales_tree_fixed$node.label <- NULL

rownames(poales_data) <- poales_data$ID

# convert to phylo4d object
#p4d <- phylo4d(poales_tree_fixed, poales_data$pheno_Topt_site_p50)
#barplot(p4d)
tree <- poales_tree_fixed
trait <- poales_data$pheno_Topt_site_p50
stopifnot(all(names(trait) == tree$tip.label) | all(tree$tip.label %in% names(trait)))

# 1. compute patristic distances
pat <- cophenetic.phylo(tree)         # symmetric matrix of tip distances

# 2. set up distance classes (choose sensible bins)
# e.g. 8-12 bins; adjust depending on tree depth and n tips
n.classes <- 100
d.vec <- as.numeric(pat[upper.tri(pat)])
breaks <- quantile(d.vec, probs = seq(0, 1, length.out = n.classes + 1))
# avoid zero-width bins by jittering if necessary
breaks <- unique(breaks)
if(length(breaks) < 3) stop("too few unique distances; reduce n.classes")

# 3. compute Moran's I (Abouheif's test) for each class
moran.results <- vector("list", length = length(breaks)-1)
names(moran.results) <- paste0("class_", seq_len(length(breaks)-1))

for(i in seq_len(length(breaks)-1)){
  cat("Processing distance class", i, "\n")
  lo <- breaks[i]
  hi <- breaks[i+1]
  # W: proximity matrix with positive entries for pairs whose patristic distance falls in [lo,hi)
  W <- matrix(0, nrow = nrow(pat), ncol = ncol(pat), dimnames = dimnames(pat))
  idx <- which(pat >= lo & pat < hi, arr.ind = TRUE)
  if(nrow(idx) > 0){
    for(r in seq_len(nrow(idx))){
      W[idx[r,1], idx[r,2]] <- 1
      W[idx[r,2], idx[r,1]] <- 1
    }
  }
  # if W is all zeros for a class, skip
  if(all(W == 0)){
    moran.results[[i]] <- NA
    next
  }
  # abouheif.moran expects a data.frame of traits; rownames must match W
  df <- as.data.frame(trait)
  rownames(df) <- names(trait)
  # run test (adjust nrepet as you like)
  moran.results[[i]] <- abouheif.moran(df, W = W, nrepet = 99, alter = "two-sided")
}
saveRDS(moran.results, file = "results/spermatophyta_moran_results.rds")
moran.results <- readRDS("results/spermatophyta_moran_results.rds")
# 4. extract Moran's I and p-values, plot correlogram
I <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$obs[1]
})

pval <- sapply(moran.results, function(x) {
  if (length(x) == 1 && is.na(x)) return(NA)
  x$pvalue[1]
})
#save all these objects 

# midpoints for plotting
mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
png("plots/spermatophyta_correlogram.png", width = 10, height = 10, units = "in", res = 300)
plot(mid, I, type = "b", pch = 16,
     main="Spermatophyte phylogenetic signal",
     xlab = "patristic distance (class midpoint)", ylab = "Moran's I (Abouheif)")
abline(h = 0, lty = 2)
points(mid, I, col = ifelse(pval < 0.05, "red", "black"))
legend("topright", legend = c("p < 0.05"), pch = 16, col = "red", bty = "n")
dev.off()
  
