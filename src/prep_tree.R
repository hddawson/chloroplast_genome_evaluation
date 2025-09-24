library(ape)
library(data.table)
library(arrow)
data <- read_parquet("data/processed_data.parquet")
barplot(table(data$Order), lax=0.45)
tree <- read.tree("data/2_global_order_level.tre")

grep("Cuscuta",data$Taxonomy)
grep("Zea",data$Taxonomy)
#the tree is order level
#GBIF.org (15 September 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.3rkm2d


tree$tip.label
table(data$Order)
setDT(data)
fasta_ids <- fread("data/tmp/cds_supermatrix.fasta", header=FALSE)
fasta_ids <- fasta_ids[grepl("^>", V1)]
fasta_ids[, ID := sub("^>", "", V1)]
map <- merge(fasta_ids[, .(ID)], data[, .(ID, Order)], by="ID")
fwrite(map[, .(ID, Order)], "data/id2order.tsv", sep="\t", col.names=FALSE)

x <- read.tree("data/2_global_order_level.tre")
x$edge.length[is.na(x$edge.length)] <- 1
write.tree(x, "data/2_global_order_level.clean.tre")

#/programs/raxml-ng_v1.2.0/raxml-ng --parse --msa data/tmp/cds_supermatrix.fasta --model GTR+G --tree data/2_global_order_level.clean.tre --force

backbone <- read.tree("data/2_global_order_level.clean.tre")
groups <- split(data$ID, data$Order)
backbone <- drop.tip(backbone, setdiff(backbone$tip.label, names(groups)))

for(ord in names(groups)){
  sub <- stree(length(groups[[ord]]), tip.label=groups[[ord]])
  backbone <- bind.tree(backbone, sub, where=which(backbone$tip.label==ord))
}

#plot(backbone, "f")
length(backbone$tip.label)             # should be ~ nrow(data)
length(unique(backbone$tip.label))     # should equal length(t$tip.label)

length(unique(data$Order))      # should be ~33
table(data$Order)               # number of samples per order

table(backbone$tip.label %in% data$ID) # all tips should be TRUE

# each order should form a clade
all(sapply(split(data$ID,data$Order),
           function(x) is.monophyletic(backbone, x)))

write.tree(backbone, "data/constraint_full.tre")

#/programs/raxml-ng_v1.2.0/raxml-ng --all --msa data/tmp/cds_supermatrix.fasta --model GTR+G --tree-constraint data/constraint_full.tre --bs-trees 100 --threads 40


t <- read.tree("data/constraint_full.tre")
aln_ids <- sub("^>", "", readLines("data/tmp/cds_supermatrix.fasta"))#[seq(1, by=2, length.out=length(t$tip.label))])
aln_ids <- aln_ids[seq(1,by=2, length.out=length(aln_ids)/2)]
length(aln_ids)
length(unique(aln_ids))
length(t$tip.label)
length(unique(t$tip.label))

setdiff(aln_ids, t$tip.label)
setdiff(t$tip.label, aln_ids)

#/programs/raxml-ng_v1.2.0/raxml-ng --check --msa ../data/tmp/cds_supermatrix.fasta --model GTR+G --prefix T1
#saves the new tree to  /local/workdir/hdd29/chloroplast_genome_evaluation/tree/T1.raxml.reduced.phy
#/programs/raxml-ng_v1.2.0/raxml-ng --parse --msa T1.raxml.reduced.phy --model GTR+G --prefix T2
# /programs/raxml-ng_v1.2.0/raxml-ng --all --msa T2.raxml.rba --tree-constraint ../data/constraint_full.tre --bs-trees 10 --threads 15
#ERROR: Following 250 taxa present in the constraint tree can not be found in the alignment:
library(phangorn)
tree <- read.tree("data/constraint_full.tre")
aln <- read.phyDat("tree/T1.raxml.reduced.phy", format="phylip") #SLOWOWOWOWOWOW
tips <- names(aln)
pruned <- drop.tip(tree, setdiff(tree$tip.label, tips))
write.tree(pruned, "tree/constraint_pruned.tre")


#/programs/raxml-ng_v1.2.0/raxml-ng --redo --msa data/tmp/cds_supermatrix.fasta --model GTR+G