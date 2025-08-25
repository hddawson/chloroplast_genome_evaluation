library(data.table)
library(ggplot2)
library(stringr)
library(matrixStats)

data <- fread("data/tmp/rbcL_aln/merged_aa_counts.csv")
sset <- c("rbcL_nuc_A", "rbcL_nuc_C", "rbcL_nuc_T", "rbcL_nuc_G",
          "F", "T", "S", "W", "R", "G", "K", "P", "M", "N", "D", "L",
          "C", "Q", "H", "V", "E", "I", "Y", "A")
sset <- data[,..sset]

par(mfrow=c(2,2))
hist(sset$rbcL_nuc_A,col="coral")
hist(sset$rbcL_nuc_T,col="pink")
hist(sset$rbcL_nuc_G,col="seagreen")
hist(sset$rbcL_nuc_C,col="blue")

# nitrogen atoms per base
N_nuc <- c(
  rbcL_nuc_A = 5,
  rbcL_nuc_C = 3,
  rbcL_nuc_G = 5,
  rbcL_nuc_T = 2
)

# nitrogen atoms per amino acid (side chain + backbone)
# (canonical free AAs; ignoring peptide bond effects)
N_aa <- c(
  F=1, T=1, S=1, W=2, R=4, G=1, K=2, P=1, M=1, N=2, D=1, L=1,
  C=1, Q=2, H=3, V=1, E=1, I=1, Y=1, A=1
)

# combine
N_map <- c(N_nuc, N_aa)

# multiply each column by its N count
N_usage <- as.data.frame(sweep(sset, 2, N_map[colnames(sset)], `*`))
colnames(N_usage) <- paste0(colnames(N_usage),"_N_usage")
# total N per sample
N_usage$total_N <- rowSums(N_usage, na.rm=TRUE)



# plot histograms of total N and components
par(mfrow=c(1,1))
hist(N_usage$total_N)

df <- cbind(N_usage, data)
df <- df[df$total_N > 1100,]
plot(df$GC_content_99_consensus,df$genomicGC)

par(mfrow=c(1,1))
boxplot(~rbcL_nuc_A+rbcL_nuc_T+rbcL_nuc_C+rbcL_nuc_G,data=df)
boxplot(df[, c("rbcL_nuc_A","rbcL_nuc_T","rbcL_nuc_C","rbcL_nuc_G")],
        main="fourfold nucleotide counts in rbcL",
        ylab="Count",
        col=c("coral","pink","blue","seagreen"))

df$FixationLabel
df$photosyntheticPathway

par(mfrow=c(2,2))
hist(df$rbcL_nuc_A,col="coral")
hist(df$rbcL_nuc_T,col="pink")
hist(df$rbcL_nuc_G,col="seagreen")
hist(df$rbcL_nuc_C,col="blue")
# colnames(data)[478:497] OG AA counts
monomer_counts <- c(colnames(df)[503:522], "rbcL_nuc_C",  "rbcL_nuc_G",  "rbcL_nuc_A", "rbcL_nuc_T")

plot(df$rbcL_nuc_A, df$total_N)
points(df[grep("Poaceae", df$Taxonomy),"rbcL_nuc_A"],
       df[grep("Poaceae", df$Taxonomy),"total_N"], col="seagreen", pch=19)
points()

a <- df[,monomer_counts]
cm <- cor(a)
library(pheatmap)
pheatmap(cor(a))

pca <- prcomp(df[,monomer_counts], scale.=TRUE)
summary(pca)
scores <- as.data.frame(pca$x)
df <- cbind(df, scores)

var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

scale_factor <- 5
loadings <- pca$rotation[,1:2] * scale_factor

par(mfrow=c(1,1))

# PCA scores
plot(df$PC1, df$PC2,
     xlab=x_lab, ylab=y_lab,
     main="PCA on AA counts + ACTG counts at 4d sites in rbcL", pch=19, col="grey")
points(df[grep("Poaceae", df$Taxonomy),"PC1"],
       df[grep("Poaceae", df$Taxonomy),"PC2"], col="seagreen", pch=19)
points(df[grep("Montiaceae", df$Taxonomy),"PC1"],
       df[grep("Montiaceae", df$Taxonomy),"PC2"], col="gold", pch=19)
points(df[grep("Cactaceae", df$Taxonomy),"PC1"],
       df[grep("Cactaceae", df$Taxonomy),"PC2"], col="red", pch=19)
par(xpd=NA)
legend("topright", inset=c(-0.2,-.1), legend=c("Other","Poaceae","Fabaceae"),
       col=c("grey","seagreen","gold"), pch=19, bty="n")
par(xpd=FALSE)

# Loadings
plot(0,0, type="n", xlim=range(loadings[,1])*1.2, ylim=range(loadings[,2])*1.2,
     xlab=x_lab, ylab=y_lab, main="")
arrows(0,0, loadings[,1], loadings[,2], length=0.1, col="red")
text(loadings[,1], loadings[,2], labels=rownames(loadings),
     pos=3, cex=0.7, col="black")
abline(h=0,v=0,lty=2,col="lightgrey")

#great, is the result same when we just look at AA's

monomer_counts <- c(colnames(df)[503:522])
pca <- prcomp(df[,monomer_counts], scale.=TRUE)
summary(pca)
scores <- as.data.frame(pca$x)
df <- cbind(df, scores)
a <- df[,c("Organism","total_N", "PC1")] # sort by 

var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

scale_factor <- 5
loadings <- pca$rotation[,1:2] * scale_factor

par(mfrow=c(1,2))

# PCA scores
plot(df$PC1, df$PC2,
     xlab=x_lab, ylab=y_lab,
     main="PCA on AA counts in rbcL", pch=19, col="grey")
points(df[grep("Poaceae", df$Taxonomy),"PC1"],
       df[grep("Poaceae", df$Taxonomy),"PC2"], col="seagreen", pch=19)
points(df[grep("Fabaceae", df$Taxonomy),"PC1"],
       df[grep("Fabaceae", df$Taxonomy),"PC2"], col="gold", pch=19)
par(xpd=NA)
legend("topright", inset=c(-0.2,-.1), legend=c("Other","Poaceae","Fabaceae"),
       col=c("grey","seagreen","gold"), pch=19, bty="n")
par(xpd=FALSE)

# Loadings
plot(0,0, type="n", xlim=range(loadings[,1])*1.2, ylim=range(loadings[,2])*1.2,
     xlab=x_lab, ylab=y_lab, main="")
arrows(0,0, loadings[,1], loadings[,2], length=0.1, col="red")
text(loadings[,1], loadings[,2], labels=rownames(loadings),
     pos=3, cex=0.7, col="black")
abline(h=0,v=0,lty=2,col="lightgrey")

monomer_counts <- c("rbcL_nuc_C",  "rbcL_nuc_G",  "rbcL_nuc_A", "rbcL_nuc_T")
pca <- prcomp(df[,monomer_counts], scale.=TRUE)
summary(pca)
scores <- as.data.frame(pca$x)
df <- cbind(df, scores)

var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

scale_factor <- 5
loadings <- pca$rotation[,1:2] * scale_factor

par(mfrow=c(1,2))

# PCA scores
plot(df$PC1, df$PC2,
     xlab=x_lab, ylab=y_lab,
     main="PCA on 4d nucleotide counts in rbcL", pch=19, col="grey")
points(df[grep("Poaceae", df$Taxonomy),"PC1"],
       df[grep("Poaceae", df$Taxonomy),"PC2"], col="seagreen", pch=19)
points(df[grep("Fabaceae", df$Taxonomy),"PC1"],
       df[grep("Fabaceae", df$Taxonomy),"PC2"], col="gold", pch=19)
par(xpd=NA)
legend("topright", inset=c(-0.2,-.1), legend=c("Other","Poaceae","Fabaceae"),
       col=c("grey","seagreen","gold"), pch=19, bty="n")
par(xpd=FALSE)

# Loadings
plot(0,0, type="n", xlim=range(loadings[,1])*1.2, ylim=range(loadings[,2])*1.2,
     xlab=x_lab, ylab=y_lab, main="")
arrows(0,0, loadings[,1], loadings[,2], length=0.1, col="red")
text(loadings[,1], loadings[,2], labels=rownames(loadings),
     pos=3, cex=0.7, col="black")
abline(h=0,v=0,lty=2,col="lightgrey")




plot(pca$x[,1],pca$x[,2], main="PCA on AA counts + ACTG counts at 4d sites in rbcL",
     xlab=x_lab, ylab=y_lab)
points(df[grep("Poaceae", df$Taxonomy),"PC1"],df[grep("Poaceae", df$Taxonomy),"PC2"],
       col="seagreen")
points(df[grep("Fabaceae", df$Taxonomy),"PC1"],df[grep("Fabaceae", df$Taxonomy),"PC2"],
       col="yellow")
loadings <- pca$rotation[,1:2]
arrows(0, 0, loadings[,1], loadings[,2], length=0.1, col="red")
text(loadings[,1], loadings[,2], labels=rownames(loadings), col="red", cex=0.7)

# create group column
df$group <- fifelse(grepl("Poaceae", df$Taxonomy), "Poaceae",
                    fifelse(grepl("Fabaceae", df$Taxonomy), "Fabaceae", "Other"))

# melt the dataframe for ggplot
library(reshape2)
monomer_counts <- c(monomer_counts, "total_N")
plot_df <- melt(df, id.vars=c("Organism","group"), measure.vars=monomer_counts,
                variable.name="Monomer", value.name="Count")

# save to PDF
pdf("plots/monomer_boxplots.pdf", width=12, height=8)

ggplot(plot_df, aes(x=group, y=Count, fill=group)) +
  geom_boxplot(outlier.shape=NA) +
  #geom_jitter(width=0.2, size=1, alpha=0.3) +
  facet_wrap(~Monomer, scales="free_y", ncol=5) +
  scale_fill_manual(values=c("Poaceae"="seagreen", "Fabaceae"="gold", "Other"="grey")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none") +
  labs(x="Group", y="Count", title="Monomer Counts by Group")

dev.off()


#let's do some mapping!!

genomic_cols <- c("rbcL_nuc_A","rbcL_nuc_C","rbcL_nuc_T","rbcL_nuc_G",
                  "F","T","S","W","R","G","K","P","M","N","D","L",
                  "C","Q","H","V","E","I","Y","A","genomicGC","genomeLength","total_N")

# phenotype columns
pheno_cols <- grep("^pheno_", colnames(df), value=TRUE)

# initialize empty list
all_results <- list()

for(pcol in pheno_cols){
  fmla <- as.formula(paste(pcol, "~", paste(genomic_cols, collapse="+")))
  lm_fit <- lm(fmla, data=df)
  sum_fit <- summary(lm_fit)
  
  coefs <- as.data.frame(sum_fit$coefficients)
  coefs <- cbind(Predictor = rownames(coefs), coefs)
  rownames(coefs) <- NULL
  
  # rename safely
  colnames(coefs)[2:5] <- c("Estimate","StdError","tValue","pValue")
  coefs$Phenotype <- pcol
  coefs$R2 <- sum_fit$r.squared
  
  all_results[[pcol]] <- coefs
}

results_df <- bind_rows(all_results)
head(results_df)

hist(-log10(results_df$pValue))


df$FixationLabel <- ifelse(df$FixationLabel == "Yes", "Yes", "No")
df$FixationLabel <- factor(df$FixationLabel, levels=c("No","Yes"))
table(df$FixationLabel)
# 3. Build formula for logistic regression
fmla <- as.formula(paste("FixationLabel ~", paste(genomic_cols, collapse="+")))

# 4. Fit logistic regression
fit <- glm(fmla, data=df, family=binomial)

# 5. Summary
summary(fit)


#pseudo R^2
1 - fit$deviance/fit$null.deviance 

df$pred_prob <- predict(fit, type="response")
y_binary <- as.integer(df$FixationLabel == "Yes")
# scatter plot of predicted probability vs actual
plot(df$pred_prob, y_binary + runif(length(y_binary), -0.05, 0.05),
     xlab = "Predicted probability",
     ylab = "Observed (1=Yes, 0=No)")

library(pROC)
roc_obj <- roc(df$FixationLabel, df$pred_prob)
plot(roc_obj)
auc(roc_obj)

df_yes <- df %>% filter(FixationLabel == "Yes")
df_no  <- df %>% filter(FixationLabel == "No")

# downsample No to match Yes
set.seed(123)

library(caret)
library(dplyr)
library(stringr)

taxonomy_lists <- str_split(df$Taxonomy, ";")
all_terms <- unlist(taxonomy_lists)
all_terms <- str_trim(all_terms)
all_terms <-unique(all_terms)
Orders <- all_terms[grep("ales", all_terms)]

get_orders <- function(taxonomy, orders) {
  terms_found <- orders[str_detect(taxonomy, orders)]
  if(length(terms_found) == 0) return(NA_character_)
  return(terms_found)
}

df$Order_found <- sapply(df$Taxonomy, get_orders, orders=Orders)
df$nOrders <- sapply(df$Order_found, length)
table(df$nOrders)

df$Taxonomy[which(df$nOrders>1)]

df$Order <- sapply(df$Order_found, function(x) if(length(x) > 0) x[1] else NA_character_)
table(df$Order)
#make a balanced dataset
df$strata <- interaction(df$FixationLabel, df$Order)

set.seed(123)
train_index <- createDataPartition(df$strata, p = 0.5, list = FALSE, times = 5)

yes_samples <- df %>% filter(FixationLabel == "Yes")
no_samples <- df %>% filter(FixationLabel == "No")

# count number of unique Orders
orders <- unique(no_samples$Order)
n_orders <- length(orders)

# target number per Order for "No" to get roughly uniform sampling
n_per_order <- floor(nrow(yes_samples) / n_orders)

# sample "No" for each Order
no_balanced <- no_samples %>%
  group_by(Order) %>%
  group_modify(~ {
    n_target <- min(n_per_order, nrow(.x))
    sample_n(.x, n_target)
  }) %>%
  ungroup()

# combine with all "Yes" samples
balanced_df <- bind_rows(yes_samples, no_balanced)

# check distribution
table(balanced_df$FixationLabel)
table(balanced_df$Order, balanced_df$FixationLabel)

fmla <- as.formula(paste("FixationLabel ~", paste(genomic_cols, collapse="+")))

# 4. Fit logistic regression
fit <- glm(fmla, data=balanced_df, family=binomial)

# 5. Summary
summary(fit)

#pseudo R^2
1 - fit$deviance/fit$null.deviance 

balanced_df$pred_prob <- predict(fit, type="response")
y_binary <- as.integer(balanced_df$FixationLabel == "Yes")
# scatter plot of predicted probability vs actual
plot(balanced_df$pred_prob, y_binary + runif(length(y_binary), -0.05, 0.05),
     xlab = "Predicted probability",
     ylab = "Observed (1=Yes, 0=No)")

roc_obj <- roc(balanced_df$FixationLabel, balanced_df$pred_prob)
plot(roc_obj)
auc(roc_obj)

######################################################
######################################################
######################################################
######################################################
######################################################
######################################################


yes_samples <- df %>% filter(FixationLabel == "Yes")
no_samples  <- df %>% filter(FixationLabel == "No")
yes_holdout_idx <- sample(nrow(yes_samples), size = floor(0.3 * nrow(yes_samples)))
yes_holdout <- yes_samples[yes_holdout_idx, ]
yes_train   <- yes_samples[-yes_holdout_idx, ]

# 2. Sample a balanced number of No samples for training
n_orders <- length(unique(no_samples$Order))
n_per_order <- floor(nrow(yes_train) / n_orders)

no_train <- no_samples %>%
  group_by(Order) %>%
  group_modify(~ {
    n_target <- min(n_per_order, nrow(.x))
    sample_n(.x, n_target)
  }) %>%
  ungroup()

# Combine training data
train_df <- bind_rows(yes_train, no_train)

# Fit logistic regression
fmla <- as.formula(paste("FixationLabel ~", paste(genomic_cols, collapse="+")))
fit <- glm(fmla, data=train_df, family=binomial)

# Training summary
summary(fit)
pseudo_R2 <- 1 - fit$deviance/fit$null.deviance
pseudo_R2

# Predict on holdout
# include all held-out Yes samples and an equal-sized No set sampled from remaining Nos
n_no_holdout <- nrow(yes_holdout)
no_holdout <- no_samples %>%
  filter(!row_number() %in% rownames(no_train)) %>% 
  sample_n(n_no_holdout)

test_df <- bind_rows(yes_holdout, no_holdout)
test_df$pred_prob <- predict(fit, newdata=test_df, type="response")
y_binary <- as.integer(test_df$FixationLabel == "Yes")

# Plot predicted probability vs observed
plot(test_df$pred_prob, y_binary + runif(length(y_binary), -0.05, 0.05),
     xlab="Predicted probability", ylab="Observed (1=Yes, 0=No)",
     main="Holdout Predictions", pch=16, col="grey")

fab_idx <- grep("Fabaceae", test_df$Taxonomy)
points(test_df$pred_prob[fab_idx],
       y_binary[fab_idx] + runif(length(fab_idx), -0.05, 0.05),
       col="forestgreen", pch=16)

# ROC and AUC
roc_obj <- roc(test_df$FixationLabel, test_df$pred_prob)
plot(roc_obj, main="ROC Curve on Holdout")
auc(roc_obj)

fab_idx  <- grep("Fabaceae", test_df$Taxonomy)
nonfab_idx <- setdiff(which(test_df$FixationLabel == "Yes"), fab_idx)

# Compute accuracy for all positives
all_pos_idx <- which(test_df$FixationLabel == "Yes")
all_pos_pred <- ifelse(test_df$pred_prob[all_pos_idx] > 0.5, "Yes", "No")
all_pos_acc <- mean(all_pos_pred == test_df$FixationLabel[all_pos_idx])

# Accuracy for Fabaceae positives
fab_pred <- ifelse(test_df$pred_prob[fab_idx] > 0.5, "Yes", "No")
fab_acc <- mean(fab_pred == test_df$FixationLabel[fab_idx])

# Accuracy for non-Fabaceae positives
nonfab_pred <- ifelse(test_df$pred_prob[nonfab_idx] > 0.5, "Yes", "No")
nonfab_acc <- mean(nonfab_pred == test_df$FixationLabel[nonfab_idx])

cat("Accuracy on all positives:", all_pos_acc, "\n")
cat("Accuracy on Fabaceae positives:", fab_acc, "\n")
cat("Accuracy on non-Fabaceae positives:", nonfab_acc, "\n")

run_once <- function(df, genomic_cols) {
  set.seed(NULL)  # for randomness each run
  
  # Split Yes into train/holdout
  yes_samples <- df %>% filter(FixationLabel == "Yes")
  no_samples  <- df %>% filter(FixationLabel == "No")
  
  yes_holdout_idx <- sample(nrow(yes_samples), size = floor(0.3 * nrow(yes_samples)))
  yes_holdout <- yes_samples[yes_holdout_idx, ]
  yes_train   <- yes_samples[-yes_holdout_idx, ]
  
  # Balance No training samples across Orders
  n_orders <- length(unique(no_samples$Order))
  n_per_order <- floor(nrow(yes_train) / n_orders)
  
  no_train <- no_samples %>%
    group_by(Order) %>%
    group_modify(~ sample_n(.x, min(n_per_order, nrow(.x)))) %>%
    ungroup()
  
  train_df <- bind_rows(yes_train, no_train)
  
  # Fit logistic regression
  fmla <- as.formula(paste("FixationLabel ~", paste(genomic_cols, collapse="+")))
  fit <- glm(fmla, data=train_df, family=binomial)
  pseudo_R2 <- 1 - fit$deviance/fit$null.deviance
  
  # Build holdout set
  n_no_holdout <- nrow(yes_holdout)
  no_holdout <- no_samples %>% 
    filter(!row_number() %in% rownames(no_train)) %>% 
    sample_n(n_no_holdout)
  
  test_df <- bind_rows(yes_holdout, no_holdout)
  test_df$pred_prob <- predict(fit, newdata=test_df, type="response")
  
  # ROC and AUC
  roc_obj <- roc(test_df$FixationLabel, test_df$pred_prob, quiet=TRUE)
  auroc <- auc(roc_obj)
  
  # Accuracy breakdown
  fab_idx  <- grep("Fabaceae", test_df$Taxonomy)
  nonfab_idx <- setdiff(which(test_df$FixationLabel == "Yes"), fab_idx)
  all_pos_idx <- which(test_df$FixationLabel == "Yes")
  
  pred <- ifelse(test_df$pred_prob > 0.5, "Yes", "No")
  
  all_pos_acc   <- mean(pred[all_pos_idx]   == test_df$FixationLabel[all_pos_idx])
  fab_acc       <- ifelse(length(fab_idx) > 0, mean(pred[fab_idx] == test_df$FixationLabel[fab_idx]), NA)
  nonfab_acc    <- ifelse(length(nonfab_idx) > 0, mean(pred[nonfab_idx] == test_df$FixationLabel[nonfab_idx]), NA)
  
  return(data.frame(
    pseudo_R2 = pseudo_R2,
    AUROC = as.numeric(auroc),
    Acc_all_pos = all_pos_acc,
    Acc_fab = fab_acc,
    Acc_nonfab = nonfab_acc
  ))
}

# Example: run 1000 replicates
res <- bind_rows(replicate(1000, run_once(df, genomic_cols), simplify=FALSE))
saveRDS(res, "data/res.RDS")
summary(res)
hist(res$Acc_nonfab)
hist(res$Acc_fab)
cor(res$Acc_nonfab, res$Acc_fab)


calib_plot <- function(df, group_name) {
  df %>%
    mutate(bin = ntile(pred_prob, 10)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(pred_prob),
      obs_rate  = mean(FixationLabel == "Yes"),
      n = n()
    ) %>%
    ggplot(aes(x = mean_pred, y = obs_rate, size = n)) +
    geom_point() +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Calibration -", group_name),
         x = "Predicted probability", y = "Observed probability") +
    theme_minimal()
}

fab_df    <- test_df %>% filter(grepl("Fabaceae", Taxonomy))
nonfab_df <- test_df %>% filter(!grepl("Fabaceae", Taxonomy))

calib_plot(fab_df, "Fabaceae")
calib_plot(nonfab_df, "non-Fabaceae")


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
