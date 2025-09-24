library(ape)
library(arrow)
tree <- read.tree("tree/T2.raxml.rba.raxml.lastTree.TMP")
plot(tree)

data <- read_parquet("data/processed_data.parquet")
