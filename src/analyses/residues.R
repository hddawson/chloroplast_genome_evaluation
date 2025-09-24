library(data.table)
library(parallel)
library(arrow)
setwd("/workdir/hdd29/chloroplast_genome_evaluation")
df        <- as.data.table(read_parquet("data/processed_data.parquet"))
hist(df$Total_Amino_Acids, main="Total AA per plastome")
hist(df$geno_genomeLength)

y <- df$pheno_Topt_site_p50
X <- as.data.table(read_parquet("data/rbcL_embeddings.parquet"))

#merge on ID in df and sample_id in X, fit a lm


#the outliers I found by manual inspection are
outliers <- c("JN032131.1", "MZ475300.1", "MH325132.1", "MK645903.1", "PV424006.1")

merged <- merge(df[, .(ID, pheno_Topt_site_p50)], X, by.x="ID", by.y="sample_id")
fit <- lm(pheno_Topt_site_p50 ~ . -ID, data=merged)
summary(fit)
plot(fit$fitted.values, merged$pheno_Topt_site_p50)

aa_cols <- intersect(c("A","R","N","D","C","Q","E","G","H","I",
                       "L","K","M","F","P","S","T","W","Y","V"),
                     names(df))

form <- as.formula(paste("pheno_Topt_site_p50 ~", paste(aa_cols, collapse="+")))
fit2 <- lm(form, data=df)
summary(fit2)
plot(fit2$fitted.values, df$pheno_Topt_site_p50)

par(mfrow=c(1,2), mar=c(4,4,3,1))

plot(fit$fitted.values, merged$pheno_Topt_site_p50,
     xlab="Fitted T_opt (rbcL embeddings)",
     ylab="Observed T_opt",
     main="rbcL Embeddings Model")

plot(fit2$fitted.values, df$pheno_Topt_site_p50,
     xlab="Fitted T_opt (AA proportions)",
     ylab="Observed T_opt",
     main="Amino Acid Composition Model")

y <- merged$pheno_Topt_site_p50
cors <- sapply(merged[, !c("ID","pheno_Topt_site_p50"), with=FALSE], function(x) cor(x, y))
par(mfrow=c(1,1))
hist(cors, main="correlations of embeddings with T_opt_site")

high_cor <- cors[abs(cors) > 0.7]
high_cor

pca <- prcomp(X[,-1])
pvar <- pca$sdev^2 * 100 / sum(pca$sdev^2) 
barplot(pvar[1:10])
plot(pca$x[,1],pca$x[,2])
plot(pca$x[,3],pca$x[,4])
plot(pca$x[,4],pca$x[,5])
plot(pca$x[,6],pca$x[,7])








dup_counts <- X[, .N, by=names(X)][N > 1]
nrow(dup_counts)        # number of distinct duplicated rows
sum(dup_counts$N)  

plot(fit$fitted.values, merged$pheno_Topt_site_p50)

train <- merge(df[Order != "Poales", .(ID, pheno_Topt_site_p50)], X, by.x="ID", by.y="sample_id")
test  <- merge(df[Order == "Poales", .(ID, pheno_Topt_site_p50)], X, by.x="ID", by.y="sample_id")

train_lm <- train[, !c("ID"), with=FALSE]
test_lm  <- test[, !c("ID"), with=FALSE]

fit <- lm(pheno_Topt_site_p50 ~ ., data=train_lm)
s   <- summary(fit)$coefficients

effects <- data.table(
  predictor = rownames(s),
  estimate  = s[, "Estimate"],
  pval      = s[, "Pr(>|t|)"]
)[predictor != "(Intercept)"]

effects_by_size <- effects[order(-abs(estimate))]
effects_by_pval <- effects[order(pval)]

head(effects_by_size)
head(effects_by_pval)

pred <- predict(fit, newdata=test_lm)
cor(pred, test$pheno_Topt_site_p50,method="spearman")
plot(pred, test$pheno_Topt_site_p50)






aln <- as.data.table(read_parquet("data/tmp/aa_supermatrix.parquet"))
aln[1:10,1:10]
unique_counts <- sapply(aln[, -1], function(x) length(unique(x[!is.na(x) & x != ""])))
par(mfrow=c(1,1))
hist(unique_counts)

maj_counts <- sapply(aln[, -1], function(x) {
  clean_x <- x[!is.na(x) & x != ""]  # Remove missing data
  if(length(clean_x) == 0) return(0)
  
  # Find modal (most common) residue
  modal_residue <- names(sort(table(clean_x), decreasing=TRUE))[1]
  
  # Count non-modal residues
  sum(clean_x == modal_residue)
})

# View results
sum(maj_counts==10961)

sum(alternate_counts > 4) / length(aln)
boxplot(alternate_counts)
summary(maj_counts)



ggplot(df, aes(x = Order, y = Total_Amino_Acids)) +
  geom_boxplot() +
  xlab("Order") +
  title("Total Amino Acids per plastome") +
  ylab("Total AA") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

plot(df$geno_genomeLength, df$Total_Amino_Acids)
points(
  df[grep("Fabales", df$Order), ][["geno_genomeLength"]],
  df[grep("Fabales", df$Order), ][["Total_Amino_Acids"]],col="yellow"
)
points(
  df[grep("Poales", df$Order), ][["geno_genomeLength"]],
  df[grep("Poales", df$Order), ][["Total_Amino_Acids"]],col="green"
)
points(
  df[grep("Pinales", df$Order), ][["geno_genomeLength"]],
  df[grep("Pinales", df$Order), ][["Total_Amino_Acids"]],col="blue"
)
points(
  df[grep("Asparagales", df$Order), ][["geno_genomeLength"]],
  df[grep("Asparagales", df$Order), ][["Total_Amino_Acids"]],col="purple"
)
