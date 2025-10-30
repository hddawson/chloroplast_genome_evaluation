library(arrow)
library(data.table)
library(stringr)

bn_embeds <- read_parquet("data/psbJ_embeddings.parquet")

bn_embeds <- as.data.table(bn_embeds)
bn_embeds[, c("ID","Gene","Taxonomy") := tstrsplit(id, "\\|")]
bn_embeds[, Gene := sub("Gene_","",Gene)]

setDT(data)
merged <- merge(bn_embeds, data[, .(ID, pheno=data$pheno_wc2.1_2.5m_bio_5_p50)], by="ID")

merged[, lapply(.SD, function(x) sum(!grepl("^-?[0-9.]+$", x))), 
       .SDcols=!(names(merged) %in% c("ID","Gene","Taxonomy","pheno"))]

emb_mat <- merged[, lapply(.SD, function(x) as.numeric(as.character(x))), 
                  .SDcols=!(names(merged) %in% c("ID","id","Gene","Taxonomy","pheno"))]
plot(colSums(is.na(emb_mat))) # is 0 

cors <- apply(as.matrix(emb_mat), 2, function(x) cor(x, merged$pheno, use="pairwise.complete.obs"))


hist(cors, main="BottleNeck-Embedding–Phenotype Correlations", xlab="Correlation")

pca <- prcomp(bn_embeds[,-1], scale=TRUE)
barplot(pca$sdev^2 / sum(pca$sdev^2))

plot(pca$x[,1],pca$x[,2])

mean_embeds <- readRDS("data/tmp/embeds_with_mds.rds")

clean_embeds <- mean_embeds[ManualOutlier==FALSE]
psbJ_mean_embeds <- clean_embeds[Gene=="psbJ"]
emb_cols <- grep("embedding", colnames(psbJ_mean_embeds))

merged <- merge(psbJ_mean_embeds, data[, .(ID, pheno=pheno_wc2.1_2.5m_bio_5_p50)], by="ID")
emb_mat <- as.matrix(merged[, ..emb_cols])
cors <- apply(emb_mat, 2, function(x) cor(x, merged$pheno, use="pairwise.complete.obs"))
hist(cors, breaks=50, main="Mean Embedding–Phenotype Correlations", xlab="Correlation")



bn <- as.data.table(read_parquet("data/psbJ_embeddings.parquet"))
bn[, c("ID","Gene","Taxonomy") := tstrsplit(id, "\\|")]
bn[, Gene := sub("Gene_","",Gene)]
bn_m <- merge(bn, data[, .(ID, pheno=pheno_wc2.1_2.5m_bio_5_p50)], by="ID")
bn_mat <- bn_m[, !c("ID","id","Gene","Taxonomy","pheno"), with=FALSE]
bn_mat <- bn_mat[, lapply(.SD, as.numeric)]
bn_cors <- apply(as.matrix(bn_mat), 2, function(x) cor(x, bn_m$pheno, use="pairwise.complete.obs"))

mean_embeds <- readRDS("data/tmp/embeds_with_mds.rds")
psbJ_mean <- mean_embeds[ManualOutlier==FALSE & Gene=="psbJ"]
emb_cols <- grep("embedding", colnames(psbJ_mean), value=TRUE)
mean_m <- merge(psbJ_mean, data[, .(ID, pheno=pheno_wc2.1_2.5m_bio_5_p50)], by="ID")
mean_mat <- as.matrix(mean_m[, ..emb_cols])
mean_cors <- apply(mean_mat, 2, function(x) cor(x, mean_m$pheno, use="pairwise.complete.obs"))

hist(mean_cors, breaks=50, col=rgb(1,0,0,0.4), xlab="Correlation", main="Bottleneck vs Mean Embedding Correlations")
hist(bn_cors, breaks=50, col=rgb(0,0,1,0.4), add=TRUE)
legend("topright", legend=c("Mean","Bottleneck"), fill=c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)), bty="n")

X_dt <- readRDS("data/tmp/x_dt_imputed.rds")
df <- read_parquet("data/processed_data.parquet"); setDT(df)

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
pred_cols <- setdiff(colnames(X_dt), c("ID","Taxonomy"))
for (c in pred_cols) if (any(is.na(X_dt[[c]]))) X_dt[is.na(get(c)), (c) := median(X_dt[[c]], na.rm=TRUE)]
y_all <- df[match(X_dt$ID, df$ID), get(pheno_col)]

gene <- "psbJ"
cols <- grep(paste0("^", gene, "__"), pred_cols, value=TRUE)
mat <- as.matrix(X_dt[, ..cols])
ok <- complete.cases(mat) & !is.na(y_all)
mat <- mat[ok,,drop=FALSE]; y <- y_all[ok]

pr <- prcomp(mat, center=TRUE, scale.=TRUE)
cor_pcs <- apply(pr$x, 2, function(x) cor(x,y))
cor_embeds <- apply(mat, 2, function(x) cor(x,y))
d <- rbind(data.frame(Type="PC",R=cor_pcs), data.frame(Type="Embedding",R=cor_embeds))

p <- ggplot(d,aes(R,fill=Type))+geom_histogram(alpha=0.6,position="identity",bins=30)+
  scale_fill_viridis_d()+labs(title=gene,x="Pearson r",y="Count")+theme_minimal()
p

xd_cors <- apply(mat, 2, function(x) cor(x,y))

hist(xd_cors, col=rgb(1,0,0,0.4), xlab="Correlation", main="Mean vs Bottleneck Correlations (psbJ)")
hist(bn_cors, col=rgb(0,0,1,0.4), add=TRUE)
legend("topright", legend=c("Mean","Bottleneck"), fill=c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)), bty="n")


bn <- merge(bn, df[, .(ID, Order)], by="ID")
pdat <- data.table(PC1=pca$x[,3], PC2=pca$x[,2], Order=bn$Order)
ggplot(pdat, aes(PC1, PC2, color=Order))+
  geom_point(size=1)+
  theme_minimal()+
  labs(title="psbJ bottleneck PCA", x="PC1", y="PC2")
                   
y <- df[match(bn$ID, df$ID), pheno_wc2.1_2.5m_bio_5_p50]
cors_pcs <- apply(pca$x, 2, function(x) cor(x, y, use="pairwise.complete.obs"))
hist(cors_pcs, breaks=30, main="PC–Temperature Correlations", xlab="Correlation")


