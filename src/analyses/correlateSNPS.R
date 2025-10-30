library(data.table)
library(arrow)
library(ggplot2)
aln <- read_parquet("data/tmp/majMinor_aln.pq")
setDT(aln)

stopifnot(is.character(aln[[1]]))
aln[is.na(aln)] <- 0
#drop all zero var columns 
invariant_cols <- names(aln)[sapply(aln, function(x) length(unique(x)) == 1)]
cat("Dropping", length(invariant_cols), "invariant columns\n")
aln <- aln[, !..invariant_cols]

pca <- prcomp(data[,-1], rank.=100, scale. = TRUE)
saveRDS(pca, "data/tmp/majMinor_aln_pca.rds")
cat("PCA done!")
quit()
data <- read_parquet("data/processed_data.parquet")
pca <- readRDS("data/tmp/majMinor_aln_pca.rds")
summary(pca)
pvar <- pca$sdev^2 / sum(pca$sdev^2)
pvar <- round(pvar*100,2)

barplot(pvar)
plot(pca$x[,1],pca$x[,2])
scores <- as.data.frame(pca$x)
df <- cbind(data,scores)

snp_cols <- grep("cds_supermatrix", colnames(aln))
cors <- sapply(aln[,..snp_cols], function(x) cor(as.numeric(x), data$pheno_wc2.1_2.5m_bio_8_p50, use="pairwise.complete.obs"))
hist(cors, main="SNP correlations with pheno_Topt_site_p50", xlab="Correlation")


which(cors==min(cors))
dev.off()
boxplot(aln$cds_supermatrix_34414, data$pheno_wc2.1_2.5m_bio_8_p50,
     pch=16, col=rgb(0,0,0,0.1), cex=0.6)
table(aln$cds_supermatrix_34414)

library(ggplot2)
ggplot(data, aes(x=factor(aln$cds_supermatrix_34414), y=pheno_wc2.1_2.5m_bio_8_p50)) +
  geom_violin(fill="gray", alpha=0.6) +
  geom_boxplot(width=0.1, outlier.size=0.5) +
  xlab("psbJ site 47 A-G") + ylab("pheno_wc2.1_2.5m_bio_8_p50")

ggplot(df, aes(x=Order, y=aln$cds_supermatrix_34414, fill=Order)) +
  geom_boxplot() +
  xlab("Order") + ylab("SNP state") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

tab <- as.data.table(table(df$Order, aln$cds_supermatrix_34414))
setnames(tab, c("Order","SNP","N"))
tab[, prop := N/sum(N), by=Order]

ggplot(tab, aes(x=Order, y=prop, fill=SNP)) +
  geom_bar(stat="identity") +
  xlab("Order") + ylab("Proportion") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggplot(data, aes(x=factor(aln$cds_supermatrix_34414), y=pheno_wc2.1_2.5m_bio_8_p50)) +
  geom_boxplot(width=0.1, outlier.size=0.5) +
  facet_wrap(~Order) +
  xlab("psbJ site 47 A-G") + ylab("pheno_wc2.1_2.5m_bio_8_p50")

#https://pmc.ncbi.nlm.nih.gov/articles/PMC10855486/
grep("Medicago sativa", df$Organism)
aln[grep("Medicago sativa", df$Organism), cds_supermatrix_34414]
df[grep("Medicago sativa", df$Organism), "Taxonomy"]


library(data.table)
library(ggplot2)

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

df <- read_parquet("data/processed_data.parquet")
setDT(df)

merged <- merge(clean_embeds[, c("ID","Gene",embed_cols), with=FALSE],
                df[, .(ID, pheno=pheno_wc2.1_2.5m_bio_8_p50)], by="ID")

gene_cors <- merged[, {
  em <- as.matrix(.SD[, ..embed_cols])
  cors <- cor(em, pheno, use="complete.obs")[,1]
  list(cor=as.vector(cors))
}, by=Gene]

hist(gene_cors$cor)

ggplot(gene_cors, aes(x=cor)) +
  geom_histogram(bins=50, fill="gray", color="black") +
  xlab("Correlation with pheno_wc2.1_2.5m_bio_8_p50") +
  ylab("Count") +
  ggtitle("Gene-level embedding–temperature correlations")



