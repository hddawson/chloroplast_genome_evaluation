library(arrow)
library(data.table)

data <- as.data.table(read_parquet("data/tmp/majMinor_aln.pq"))

aln_cols <- grep("cds_supermatrix", colnames(data))
data$V1[1:10]

pca <- prcomp(data[,..aln_cols], scale.=TRUE, rank.=100)
saveRDS(pca, "data/tmp/aln_majMinor_pca.rds")