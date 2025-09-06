library(data.table)

df        <- fread(file.path(data_dir, "process"))
y         <- df$pheno_Topt_site_p50
hold      <- df$Order %in% c("Fabales", "Poales")
yNA       <- y
yNA[hold] <- NA

data <- as.data.table(read_parquet("data/processed_data.parquet"))
df <- data[!is.na(data$pheno_Topt_site_p50),]
sum(!is.na(data$pheno_Topt_site_p50))
sum(is.na(data$pheno_wc))
pca <- readRDS("data/tmp/onehot_aln_pca.rds")
pca_scores <- pca$x

table(df$Order)

df <-cbind(df, pca_scores)

mod <- lm(pheno_Topt_site_p50 ~ PC1 + PC2 + PC3 + PC4 + PC5
          + PC6 + PC7 + PC8 + PC9 + PC10 + 
            PC11 + PC12 + PC13 + PC14 + PC15
          + PC16 + PC17 + PC18 + PC19 + PC20 + geno_genomicGC, data=df)
summary(mod)
plot(mod)
plot(mod$fitted.values, df$pheno_Topt_site_p50)

mod <- lm(pheno_Topt_site_p50 ~ PC1 + PC2 + PC3 + PC4 + PC5
          + PC6 + PC7 + PC8 + PC9 + PC10 + 
            PC11 + PC12 + PC13 + PC14 + PC15
          + PC16 + PC17 + PC18 + PC19 + PC20 + geno_genomicGC,
          data=df[grep("Poaceae", df$Taxonomy),])
summary(mod)
plot(mod$fitted.values, df[grep("Poaceae", df$Taxonomy),]$pheno_Topt_site_p50)

mod <- lm(pheno_Topt_site_p50 ~ PC1 + PC2 + PC3 + PC4 + PC5
          + PC6 + PC7 + PC8 + PC9 + PC10 + 
            PC11 + PC12 + PC13 + PC14 + PC15
          + PC16 + PC17 + PC18 + PC19 + PC20 + geno_genomicGC,
          data=df[grep("Poaceae", df$Taxonomy),])
plot(mod$fitted.values, df[grep("Poaceae", df$Taxonomy),]$pheno_Topt_site_p50)




plot(mod)

mod <- lm(pheno_Topt_site_p50 ~ PC1 + PC2 + PC3 + PC4 + PC5
          + PC6 + PC7 + PC8 + PC9 + PC10 + 
            PC11 + PC12 + PC13 + PC14 + PC15
          + PC16 + PC17 + PC18 + PC19 + PC20 + geno_genomicGC,
          data=df[grep("Fabaceae", df$Taxonomy),])
summary(mod)

mod <- lm(pheno_Topt_site_p50 ~ PC1 + PC2 + PC3 + PC4 + PC5
          + PC6 + PC7 + PC8 + PC9 + PC10 + 
            PC11 + PC12 + PC13 + PC14 + PC15
          + PC16 + PC17 + PC18 + PC19 + PC20 + geno_genomicGC,
          data=df[grep("Brassicales", df$Taxonomy),])
summary(mod)


dinuc_cols <- grep("^geno_dinucleotide_", names(df), value = TRUE)
pc_cols <-grep("^PC", names(df), value = TRUE)
aa_codes <- c("A", "R", "N", "D", "C", 
              "Q", "E", "G", "H", "I", 
              "L", "K", "M", "F", "P", 
              "S", "T", "W", "Y", "V")
# Create formula: response ~ all matching predictors
f <- reformulate(termlabels = c(dinuc_cols, aa_codes, pc_cols[1:20]), response = "pheno_Topt_site_p50")
f <- reformulate(termlabels = pc_cols, response = "pheno_Topt_site_p50")

# Fit model
model <- lm(f, data = df)
summary(model)
#plot(model)
plot(model$fitted.values, df$pheno_Topt_site_p50)
abline(a=0,b=1, col="red")

model <- lm(f, data = df[grep("Poaceae", df$Taxonomy),])
summary(model)
#plot(model)
plot(model$fitted.values, df[grep("Poaceae", df$Taxonomy),]$pheno_Topt_site_p50)
abline(a=0,b=1, col="red")



