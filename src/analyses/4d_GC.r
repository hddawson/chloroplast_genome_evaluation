library(data.table)
library(ggplot2)
library(stringr)
library(matrixStats)

data <- fread("data/selected_genomes.csv")

n_pos <- sum(str_detect(data$Taxonomy, "Poaceae"))

gc_cols <- grep("degenerate", colnames(data), value=TRUE)
sum(is.na(df[, ..gc_cols]))
valid_gc_cols <- sapply(data[, ..gc_cols], function(x) !all(is.na(x)))
gc_cols <- gc_cols[valid_gc_cols]


for (col in gc_cols) {
  nas <- which(is.na(data[[col]]))
  if(length(nas) > 0) {
    set(data, i=nas, j=col, value=median(data[[col]], na.rm=TRUE))
  }
}

gcs <- data[, ..gc_cols]
sum(is.na(gcs))
library(pheatmap)
pheatmap(cor(gcs))

plants <- data[which(str_detect(data$Taxonomy, "Tracheophyta"))]
gcs <- plants[, ..gc_cols]
cor_mat <- cor(gcs)

high_cor_cluster <- c("psaA", "rpoC2","psbB", "psbD", "rbcL", "rpoB", "atpA", "psaB")
high_cor_cols <- grep(paste(high_cor_cluster, collapse="|"), gc_cols, value=TRUE)

pairs(gcs[,..high_cor_cols])

pos <- data[which(str_detect(data$Taxonomy, "Poales"))]
sum(is.na(pos))
sum(is.na(pos[,..gc_cols]))

vars <- colVars(as.matrix(pos[, ..gc_cols]))

# Keep columns with variance > 0
non_const_cols <- gc_cols[vars > 0]

# Compute correlation and plot heatmap
cor_mat <- cor(pos[, ..non_const_cols])
pheatmap(cor_mat)


data <- data[, !duplicated(colnames(data)), with=FALSE]
gc_cols <- grep("_degenerate_GC$", colnames(data), value=TRUE)
trna_cols <- grep("^trn.{5}$", colnames(data), value=TRUE)

results <- list()
for (gc in gc_cols) {
  safe_trna_cols <- paste0("`", trna_cols, "`")
  formula <- as.formula(paste(gc,"~", paste(safe_trna_cols, collapse = "+")))
  model <- lm(formula, data=data)
  coefs <- summary(model)$coefficients
  # Skip intercept row, keep trna effects
  coefs_trna <- coefs[rownames(coefs) %in% safe_trna_cols, , drop=FALSE]
  if (nrow(coefs_trna) > 0) {
    res_dt <- data.table(
      gene = gc,
      trna = rownames(coefs_trna),
      estimate = coefs_trna[, "Estimate"],
      pvalue = coefs_trna[, "Pr(>|t|)"]
    )
    results[[gc]] <- res_dt
  }
}
final_res <- rbindlist(results)
# Adjust p-values (BH FDR)
final_res[, p_adj := p.adjust(pvalue, method="BH")]
final_res <- final_res[order(p_adj)]

sig_res <- final_res[p_adj < 0.05]

# Summarize counts of significant tRNAs per gene
table(sig_res$gene)

library(dplyr)
sig_res %>% 
  group_by(gene) %>% 
  arrange(p_adj) %>% 
  slice_head(n=5) %>% 
  print()

hist(-log10(sig_res$p_adj))

summary_df <- final_res %>%
  filter(p_adj < 0.05) %>%
  group_by(trna) %>%
  summarize(
    n_genes = n(),
    prop_positive = mean(estimate > 0),
    prop_negative = mean(estimate < 0)
  ) %>%
  arrange(desc(n_genes))

print(summary_df)

pca <- prcomp(gcs, scale.=TRUE)
scores <- as.data.frame(pca$x)
scores$Poaceae <- str_detect(plants$Taxonomy, "Poaceae")
plot(scores$PC1, scores$PC2)
var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

ggplot(scores, aes(x=PC1, y=PC2, color=Poaceae)) +
  geom_point(alpha=0.7) +
  theme_minimal() +
  labs(title="PCA of gcs data", color="Poaceae", x=x_lab, y=y_lab)



ggplot(scores, aes(x=PC1, y=PC2, color=Poaceae, label=plants$Organism)) +
  geom_point(alpha=0.7) +
  geom_text(check_overlap = TRUE, vjust = -1, size = 3) +
  coord_cartesian(xlim = c(2.5, 3.5), ylim = c(-5, -4.5)) +
  theme_minimal() +
  labs(title="PCA of gcs data", color="Poaceae")

loadings <- as.data.frame(pca$rotation)
head(loadings)


pc_trn <- prcomp(data[,..trna_cols], scale=T)
summary(pc_trn)
plot(pc_trn$x[,1],pc_trn$x[,2])

pc_trn <- prcomp(plants[,..trna_cols], scale=T)
summary(pc_trn)
plot(pc_trn$x[,1],pc_trn$x[,2])

scores <- as.data.frame(pc_trn$x)
scores$Poaceae <- str_detect(plants$Taxonomy, "Poaceae")

ggplot(scores, aes(x=PC1, y=PC2, color=Poaceae)) +
  geom_point(alpha=0.7) +
  theme_minimal() +
  labs(title="PCA of trn gain loss data", color="Poaceae")

binary_gene_cols <- colnames(data)[sapply(data, function(x) all(na.omit(x) %in% c(0,1)))]
binary_gene_cols <- binary_gene_cols[1:120]


pc_trn <- prcomp(plants[,..binary_gene_cols], scale=T)
summary(pc_trn)
plot(pc_trn$x[,1],pc_trn$x[,2])

scores <- as.data.frame(pc_trn$x)
scores$Poaceae <- str_detect(plants$Taxonomy, "Poaceae")

var_exp <- pc_trn$sdev^2 / sum(pc_trn$sdev^2) * 100
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

ggplot(scores, aes(x=PC1, y=PC2, color=Poaceae, label=plants$Organism)) +
  geom_point(alpha=0.7) +
  geom_text(check_overlap=TRUE, vjust=-1, size=3) +
  theme_minimal() +
  labs(title="PCA of gene gain loss data", color="Poaceae",x=x_lab, y=y_lab)


loadings <- as.data.frame(pca$rotation)
head(loadings)

sum(plants$petD)

nitrogix <- plants[plants$FixationLabel!=""]
safe_gc_cols <- paste0("`", gc_cols, "`")

formula <- as.formula(paste("FixationLabel ~", paste(safe_gc_cols, collapse = "+")))
nitrogix
# Fit logistic regression (FixationLabel should be factor with 2 levels)
nitrogix[, FixationLabel := factor(FixationLabel)]
nitrogix_scaled <- copy(nitrogix)
nitrogix_scaled[, (gc_cols) := lapply(.SD, scale), .SDcols = gc_cols]

# optionally remove row 1645
nitrogix_scaled <- nitrogix_scaled[c(1645,324,1511,1793)]$Taxonomy

model <- glm(formula, datga=nitrogix_scaled, family=binomial)

summary(model)
plot(model)
confint(model)                # confidence intervals for coefficients  
anova(model, test="Chisq")   # likelihood ratio tests for terms  
library(car)  
vif(model)  

# For each predictor, create a contingency table with FixationLabel
perfect_preds <- sapply(gc_cols, function(col) {
  vals <- nitrogix[[col]]
  label <- nitrogix$FixationLabel
  groups <- vals > median(vals, na.rm=TRUE)
  tbl <- table(groups, label)
  any(apply(tbl, 1, function(x) (x[1] == 0 && x[2] > 0) || (x[2] == 0 && x[1] > 0)))
})

names(perfect_preds[perfect_preds])

# Check for predictors with near-zero variance
near_zero_var <- sapply(gc_cols, function(col) {
  var(nitrogix[[col]], na.rm=TRUE) < 1e-6
})
names(near_zero_var[near_zero_var])

# Check correlation with response if it's numeric (if factor, use point-biserial)
library(Hmisc)
rcorr_matrix <- rcorr(as.matrix(nitrogix[, ..gc_cols]), as.numeric(nitrogix$FixationLabel))
head(rcorr_matrix$r)

# Or plot predictors against outcome to see weird patterns
#gc_cols <- grep("_degenerate_GC$", colnames(nitrogix), value=TRUE)
cols <- c("blue", "red")

for (col in gc_cols) {
  data1 <- nitrogix[FixationLabel == levels(factor(FixationLabel))[1], ][[col]]
  data2 <- nitrogix[FixationLabel == levels(factor(FixationLabel))[2], ][[col]]
  
  hist(data1, col=adjustcolor(cols[1], alpha.f=0.5), main=col, xlab="GC Content", ylab="Frequency", 
       breaks=20, xlim=range(nitrogix[[col]], na.rm=TRUE), 
       ylim=c(0, max(hist(data1, plot=FALSE)$counts, hist(data2, plot=FALSE)$counts)))
  hist(data2, col=adjustcolor(cols[2], alpha.f=0.5), add=TRUE, breaks=20)
  legend("topright", legend=levels(factor(nitrogix$FixationLabel)), fill=adjustcolor(cols, alpha.f=0.5))
}

boxplot(nitrogix$genomicGC~nitrogix$FixationLabel)



######

dir()
df <- fread("data/selected_genomes.csv")

fourFoldgc_cols <- grep("degenerate", colnames(df), value=TRUE)

sum(is.na(df[, ..gc_cols]))

valid_gc_cols <- sapply(data[, ..gc_cols], function(x) !all(is.na(x)))
gc_cols <- gc_cols[valid_gc_cols]

#plot histogram of proportion NA per column

nas <- c()
for (gc_col in gc_cols) {
  nones <- sum(is.na(df[, ..gc_col])) / nrow(df)
  nas <- append(nas, nones)
}

goods <- nas < 0.1
sum(goods) / length(nas)

good_cols <- gc_cols[goods]

df2 <- copy(df)
for (col in good_cols) {
  nas <- which(is.na(df2[[col]]))
  if(length(nas) > 0) {
    set(df2, i=nas, j=col, value=median(df2[[col]], na.rm=TRUE))
  }
}

pca <- prcomp(df2[,..good_cols])
scores <- data.frame(pca$x)
summary(pca)
plot(scores$PC1, scores$PC2)
points(scores[grep("Poaceae", df2$Taxonomy), "PC1"],
       scores[grep("Poaceae", df2$Taxonomy), "PC2"],
       col="coral")

points(scores[grep("C4", df2$photosyntheticPathway), "PC1"],
       scores[grep("C4", df2$photosyntheticPathway), "PC2"],
       col="lightgreen")

#add legend

text(scores[grep("Zea mays", df2$Organism), "PC1"],
     scores[grep("Zea mays", df2$Organism), "PC2"], labels="Zea mays", 
     col="cornflowerblue")

missing_genes <- c("accD", "ycf1", "yc2", "rpl22", "infA")

good_cols

good_genes <- sapply(good_cols, function(a) str_split_i(a, "_", 1))

pca <- prcomp(df2[,..good_genes])
scores <- data.frame(pca$x)
summary(pca)
plot(scores$PC1, scores$PC2)
points(scores[grep("Poaceae", df2$Taxonomy), "PC1"],
       scores[grep("Poaceae", df2$Taxonomy), "PC2"],
       col="coral")

points(scores[grep("Viridiplantae", df2$Taxonomy), "PC1"],
       scores[grep("Viridiplantae", df2$Taxonomy), "PC2"],
       col="green")

loadings <- pca$rotation
print(loadings)

#1026 other samples than the 18044 oplants - but they are alllll over in the middle 
plants <- df2[grep("Spermatophyta", df2$Taxonomy),]
grep("Zea", plants$Taxonomy)
pca <- prcomp(plants[,..good_cols])
scores <- data.frame(pca$x)
summary(pca)
plot(scores$PC1, scores$PC2)
points(scores[grep("Poaceae", df2$Taxonomy), "PC1"],
       scores[grep("Poaceae", df2$Taxonomy), "PC2"],
       col="coral")

length(grep("Spermatophyta", df2$Taxonomy))


pca <- prcomp(plants[,..good_cols])
scores <- data.frame(pca$x)
summary(pca)
plot(scores$PC1, scores$PC2)
points(scores[grep("Poaceae", plants$Taxonomy), "PC1"],
       scores[grep("Poaceae", plants$Taxonomy), "PC2"],
       col="coral")

points(scores[grep("eudicotyledons", plants$Taxonomy), "PC1"],
       scores[grep("eudicotyledons", plants$Taxonomy), "PC2"],
       col="pink")

points(scores[grep("BOP clade", plants$Taxonomy), "PC1"],
       scores[grep("BOP clade", plants$Taxonomy), "PC2"],
       col="darkgreen")

points(scores[grep("PACMAD", plants$Taxonomy), "PC1"],
       scores[grep("PACMAD", plants$Taxonomy), "PC2"],
       col="lightgreen")

points(scores[grep("Triticodae", plants$Taxonomy), "PC1"],
       scores[grep("Triticodae", plants$Taxonomy), "PC2"],
       col="red")

points(scores[grep("Yes", plants$FixationLabel), "PC1"],
       scores[grep("Yes", plants$FixationLabel), "PC2"],
       col="yellow")

points(scores[grep("No", plants$FixationLabel), "PC1"],
       scores[grep("No", plants$FixationLabel), "PC2"],
       col="red")

df2[grep("Zea mays", plants$Organism), "Taxonomy"]
points(scores[grep("Zea mays", plants$Organism), "PC1"],
     scores[grep("Zea mays", plants$Organism), "PC2"],
     col="cornflowerblue")

points(scores[grep("C4", plants$photosyntheticPathway), "PC1"],
       scores[grep("C4", plants$photosyntheticPathway), "PC2"],
       col="lightgreen")

points(scores[grep("C4", plants$photosyntheticPathway), "PC1"],
       scores[grep("C4", plants$photosyntheticPathway), "PC2"],
       col="lightgreen")


randDF <- copy(df)
for (col in good_cols) {
  nas <- which(is.na(randDF[[col]]))
  if(length(nas) > 0) {
    observed <- randDF[[col]][!is.na(randDF[[col]])]
    set(randDF, i=nas, j=col, value=sample(observed, length(nas), replace=TRUE))
  }
}

plants <- randDF[grep("Spermatophyta", randDF$Taxonomy),]
grep("Zea", plants$Taxonomy)
pca <- prcomp(plants[,..good_cols])
scores <- data.frame(pca$x)

plot(scores$PC1, scores$PC2)
points(scores[grep("Poaceae", plants$Taxonomy), "PC1"],
       scores[grep("Poaceae", plants$Taxonomy), "PC2"],
       col="coral")

points(scores[grep("Triticodae", plants$Taxonomy), "PC1"],
       scores[grep("Triticodae", plants$Taxonomy), "PC2"],
       col="cornflowerblue")

points(scores[grep("Zea mays", plants$Organism), "PC1"],
       scores[grep("Zea mays", plants$Organism), "PC2"],
       col="yellow")

plants <- randDF[grep("Spermatophyta", randDF$Taxonomy),]
pca <- prcomp(plants[,..good_cols], scale.=TRUE)
scores <- data.frame(pca$x)
ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

plot(scores$PC1, scores$PC2,
     xlab=paste0("PC1 (", ve[1], "% variance)"),
     ylab=paste0("PC2 (", ve[2], "% variance)"),
     main="PCA of Degenerate GC in Genes with >90% Conservation in Spermatophyta",
     col="black",cex=0.7)

points(scores[grep("Poales", plants$Taxonomy), c("PC1","PC2")], col="coral", pch=16, cex=1.4)
points(scores[grep("Poaceae", plants$Taxonomy), c("PC1","PC2")], col="tomato", pch=16)
points(scores[grep("Fabaceae", plants$Taxonomy), c("PC1","PC2")], col="forestgreen", pch=16)

points(scores[grep("Triticodae", plants$Taxonomy), c("PC1","PC2")], col="cornflowerblue", pch=16)
points(scores[grep("Zea mays", plants$Organism), c("PC1","PC2")], col="yellow", pch=16)

legend("bottomleft",
       legend=c("Poales","Poaceae","Triticodae","Zea mays","Fabaceae"),
       col=c("coral","tomato","cornflowerblue","yellow","forestgreen"),
       pch=c(16,16,16,16,16),
       pt.cex=c(1,1,1,1,1))

par(xpd=TRUE, mar=c(5, 4, 4, 8))
plot(scores$PC1, scores$PC2,
     xlab=paste0("PC1 (", ve[1], "% variance)"),
     ylab=paste0("PC2 (", ve[2], "% variance)"),
     main="PCA of Degenerate GC in Genes with >90% Conservation in Spermatophyta",
     col="black", cex=0.7)

points(scores[grep("Poales", plants$Taxonomy), c("PC1","PC2")], col="coral", pch=16, cex=1.4)
points(scores[grep("Poaceae", plants$Taxonomy), c("PC1","PC2")], col="tomato", pch=16)
points(scores[grep("Fabaceae", plants$Taxonomy), c("PC1","PC2")], col="forestgreen", pch=16)
points(scores[grep("Triticodae", plants$Taxonomy), c("PC1","PC2")], col="cornflowerblue", pch=16)
points(scores[grep("Zea mays", plants$Organism), c("PC1","PC2")], col="yellow", pch=16)

legend("topright", inset=c(-0.25,0),
       legend=c("Poales","Poaceae","Triticodae","Zea mays","Fabaceae"),
       col=c("coral","tomato","cornflowerblue","yellow","forestgreen"),
       pch=16, pt.cex=1.2)




extreme_ids <- unique(c(
  which.max(scores$PC1),
  which.min(scores$PC1),
  which.max(scores$PC2),
  which.min(scores$PC2)
))

plants[extreme_ids, c("Organism", "Taxonomy")]
