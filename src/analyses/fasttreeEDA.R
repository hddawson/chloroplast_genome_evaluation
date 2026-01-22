library(ape)
library(phytools)
library(arrow)
library(phangorn)
library(Biostrings)

# load tree
tr <- read.tree("iqtree_input/fast.tree")
data <- read_parquet("data/processed_data.parquet")

## 1. Basic tree stats
cat("Tips:", Ntip(tr), "\n")
cat("Nodes:", tr$Nnode, "\n")
cat("Tree height:", max(nodeHeights(tr)), "\n")
cat("Total branch length:", sum(tr$edge.length), "\n")
is.rooted(tr)
head(tr$tip.label)
data$ID[grep("Pinales",data$Order)]


is.monophyletic(tr,data$ID[grep("Pinales",data$Order)])
pine_mrca <- getMRCA(tr,data$ID[grep("Pinales",data$Order)])
rooted_tree <- root(tr, node=pine_mrca)

pine_ids <- data$ID[data$Order=="Pinales"]
pine_labels <- setNames(data$Organism, data$ID)

pine_tree <- extract.clade(tr, pine_mrca)
pine_tree$tip.label <- pine_labels[pine_tree$tip.label]

plot(pine_tree, cex=0.6)

data$ID[grep("Pinales",data$Order)]

#aln = read.FASTA("iqtree_input/supergene.faa",type = "AA")
#aln = read.phyDat("iqtree_input/supergene.faa",format = "fasta",type="AA")
aln = readAAStringSet("iqtree_input/supergene.faa",format = "fasta")
phyAln <- as.phyDat(aln)
head(rooted_tree$tip.label)
head(names(phyAln))
setdiff(rooted_tree$tip.label, names(phyAln))
setdiff(names(phyAln), rooted_tree$tip.label)
any(duplicated(rooted_tree$tip.label))

any(is.na(rooted_tree$edge.length))
any(rooted_tree$edge.length <= 0)
hist(rooted_tree$edge.length)
rooted_tree$edge.length[rooted_tree$edge.length <= 0] <- 1e-8

is.binary(rooted_tree)
rooted_tree <- multi2di(rooted_tree)
rooted_tree$edge.length[rooted_tree$edge.length <= 0] <- 1e-8

standard_aa <- c("A","R","N","D","C","Q","E","G","H","I",
                 "L","K","M","F","P","S","T","W","Y","V","-")

# Check each sequence
has_nonstandard <- sapply(names(phyAln), function(tip) {
  seq_chars <- as.character(phyAln[[tip]])
  any(!seq_chars %in% standard_aa)
})

cat("Sequences with nonstandard AAs:", sum(has_nonstandard), "\n")

aln_matrix <- as.matrix(aln)

# Define standard AAs
standard_aa <- c("A","R","N","D","C","Q","E","G","H","I",
                 "L","K","M","F","P","S","T","W","Y","V","-")

# Replace nonstandard AAs with gaps
aln_matrix[!aln_matrix %in% standard_aa] <- "-"

# Check what was replaced
cat("Unique characters after cleaning:", 
    paste(sort(unique(as.vector(aln_matrix))), collapse=","), "\n")

# Convert to phyDat
phyAln_clean <- phyDat(aln_matrix, type="AA")

# Match with tree
shared_tips <- intersect(rooted_tree$tip.label, names(phyAln_clean))
cat("Shared tips:", length(shared_tips), "\n")
stopifnot(length(shared_tips) > 0)

rooted_tree_clean <- keep.tip(rooted_tree, shared_tips)
phyAln_clean <- subset(phyAln_clean, subset=shared_tips)

stopifnot(all(rooted_tree_clean$tip.label %in% names(phyAln_clean)))


set.seed(123)
sample_tips <- sample(shared_tips, min(9000, length(shared_tips)))
cat("Subsampled to:", length(sample_tips), "tips\n")

rooted_tree_sub <- keep.tip(rooted_tree_clean, sample_tips)
phyAln_sub <- subset(phyAln_clean, subset=sample_tips)

stopifnot(all(rooted_tree_sub$tip.label %in% names(phyAln_sub)))
stopifnot(Ntip(rooted_tree_sub) == length(phyAln_sub))

# Try pml
fit <- pml(rooted_tree_sub, phyAln_sub)
#got the subsample!
fit
cat("FITTED SMALL!")
# Try pml
fit <- pml(rooted_tree_clean, phyAln_clean)

cat("FITTING tree!")

#fit <- pml(rooted_tree, data = phyAln)
saveRDS("data/tmp/fit.rds")
# Optimize under a GTR+Γ+I model (common for DNA)
fit_opt <- optim.pml(
  fit,
  model = "LG",
  optInv = TRUE,
  optEdge = TRUE
)
saveRDS("data/tmp/fit_opt.rds")



## 2. Branch length distribution
bl <- tr$edge.length
summary(bl)

hist(bl, main = "Branch length distribution",
     xlab = "Branch length")

## 3. Terminal vs internal branches
is_tip_edge <- tr$edge[,2] <= Ntip(tr)
boxplot(
  bl[is_tip_edge],
  bl[!is_tip_edge],
  names = c("Terminal", "Internal"),
  main = "Terminal vs internal branch lengths",
  ylab = "Branch length"
)

## 4. Root-to-tip distances (clock signal check)
rtt <- nodeHeights(tr)
tip_rtt <- rtt[rtt[,2] <= Ntip(tr), 2]

hist(tip_rtt,
     main = "Root-to-tip distances",
     xlab = "Distance")

cat("RTT CV:", sd(tip_rtt) / mean(tip_rtt), "\n")

## 5. Extremely long branches (outliers)
q <- quantile(bl, 0.99)
long_edges <- which(bl > q)
cat("Edges >99th percentile:", length(long_edges), "\n")

long_edges <- which(bl > quantile(bl, 0.995))
long_tips <- tr$edge[long_edges, 2]
long_tips <- long_tips[long_tips <= Ntip(tr)]

head(tr$tip.label[long_tips])
length(long_tips)
boxplot(tr$edge.length, main="Edge length distribution")

library(arrow)

data <- read_parquet("data/processed_data.parquet") 

data$Taxonomy[match(tr$tip.label[long_tips], data$ID)]
data$Organism[match(tr$tip.label[long_tips], data$ID)]


# ---- make sure tips match ----

tr <- read.tree("iqtree_input/fast.tree")
al <- read.FASTA("iqtree_input/supergene.faa")

# tips in tree not in alignment
setdiff(tr$tip.label, names(al))

# tips in alignment not in tree
setdiff(names(al), tr$tip.label)


