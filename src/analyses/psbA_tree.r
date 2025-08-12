library(ape)
library(stringr)
library(ggtree)
tree <- read.tree("/workdir/hdd29/chloroplast_genome_evaluation/data/psbA.fa.tree")
Ntip(tree)
ggtree(tree) + geom_tiplab(aes(label=label), size=1)

tax_data <- read.csv("/workdir/hdd29/chloroplast_genome_evaluation/data/taxonomy_info.csv")

str_split_i(tree$tip.label,"_",2)[1:5]

tax_data$isPlant <- str_detect(tax_data$Taxonomy,"Spermatophyta")

acc <- str_split_i(tree$tip.label, "_", 2)
tip_tax <- tax_data[match(acc, tax_data$FileBasename), ]

cols <- ifelse(tip_tax$isPlant, "forestgreen", "orange")

clades <- list(
  plants = which(tip_tax$isPlant),
  nonplants = which(!tip_tax$isPlant)
)

tree_collapsed <- tree
for (grp in clades) {
  mrca_node <- getMRCA(tree_collapsed, tree_collapsed$tip.label[grp])
  if (!is.null(mrca_node)) tree_collapsed <- collapse.singles(drop.tip(tree_collapsed, setdiff(1:Ntip(tree_collapsed), grp)))
}

pdf("tree_isPlant.pdf", width=12, height=12)
plot(tree_collapsed, tip.color=cols, cex=0.3, no.margin=TRUE)
dev.off()

collapse_by_trait <- function(tr, trait) {
  repeat {
    collapsed <- FALSE
    for (node in (Ntip(tr)+1):max(tr$edge)) {
      tips <- extract.clade(tr, node)$tip.label
      idx <- match(tips, tr$tip.label)
      if (length(unique(trait[idx])) == 1) {
        keep_tip <- idx[1]
        drop_idx <- setdiff(idx, keep_tip)
        tr <- drop.tip(tr, drop_idx)
        trait <- trait[-drop_idx]
        collapsed <- TRUE
        break
      }
    }
    if (!collapsed) break
  }
  list(tree=tr, trait=trait)
}

res <- collapse_by_trait(tree, tip_tax$isPlant)

pdf("tree_isPlant_collapsed.pdf", width=12, height=12)
plot(res$tree, tip.color=ifelse(res$trait, "forestgreen", "orange"), cex=0.3, no.margin=TRUE)
dev.off()


