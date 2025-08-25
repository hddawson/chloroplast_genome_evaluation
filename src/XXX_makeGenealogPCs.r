library(data.table)

data <- fread("data/tmp/aa_supermatrix_expanded_3pcs.csv")

stopifnot(is.character(data[[1]]))
#drop all zero var columns 
invariant_cols <- names(data)[sapply(data, function(x) length(unique(x)) == 1)]
cat("Dropping", length(invariant_cols), "invariant columns\n")
data <- data[, !..invariant_cols]

pca <- prcomp(data[,-1], rank.=100, scale. = TRUE)
saveRDS(pca, "data/tmp/aa_supermatrix_expanded_10pcs_pca.rds")
cat("PCA done!")

pca <- readRDS("data/tmp/aa_supermatrix_expanded_10pcs_pca.rds")
summary(pca)
plot(pca)
plot(pca$x[,1],pca$x[,2])
scores <- as.data.frame(pca$x)
df <- cbind(data,scores)
df$ID <- df$V1
tax_data <- fread("data/tmp/rbcL_aln/merged_aa_counts.csv")
head(df$ID)
head(tax_data$ID)
setDT(tax_data)
setDT(df)
df <- df[tax_data, on = "ID"] 

with(df[grep("Poaceae", Taxonomy), ],
     points(PC1, PC2, col="seagreen"))

df[order(PC1), .(Organism, PC1)]
#I get wierd results for 
#Organism        PC1
#<char>      <num>
#  1:   Pontederia cordata var. lancifolia -1278.9312
#2:                       Solanum aureum -1070.0965
#3: Eremophila phyllopoda subsp. obliqua  -876.8482
#4:                   Solanum glutinosum  -872.7635. - also nice floweing plant, I can't determine the distinction
#5:                         Senna siamea  -795.1061
#---  
#Pontderia cordata is an aquatic plant - stems underwater
#Eremophila phyllopoda - soime sort of shrub, " grey, pendulous leaves" 
plot(df$PC1, df$total_N)
loadings <- pca$rotation

top_pc1 <- sort(loadings[, "PC1"], decreasing = TRUE)[1:10]
top_pc1_neg <- sort(loadings[, "PC1"], decreasing = FALSE)[1:10]

top_pc1
top_pc1_neg


plot(pca$x[,3],pca$x[,4])
with(df[df[["rps18"]]==0, ],
     points(PC1, PC2, col="red"))

plot(df$PC1,df$PC2)
with(df[grep("Poaceae", Taxonomy), ],
     points(PC3, PC4, col="seagreen"))
with(df[grep("Fabales", Taxonomy), ],
     points(PC3, PC4, col="yellow"))


plot(pca$x[,1],pca$x[,2])
with(df[grep("Poaceae", Taxonomy), ],
     points(PC2, PC3, col="seagreen"))

var_exp <- pca$sdev^2
var_exp_pct <- round(100 * var_exp / sum(var_exp), 1)

# create a group column for plotting
df <- df %>%
  mutate(Order_group = case_when(
    grepl("Poaceae", Taxonomy) ~ "Poaceae",
    grepl("Fabales", Taxonomy) ~ "Fabales",
    TRUE ~ "Other"
  ))

# plot
df$Order_group <- factor(df$Order_group, levels = c("Other", "Fabales", "Poaceae"))

# plot
ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(data = df %>% filter(Order_group == "Other"),
             color = "black", alpha = 0.6, size = 1.5) +
  geom_point(data = df %>% filter(Order_group != "Other"),
             aes(color = Order_group), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Poaceae" = "seagreen",
                                "Fabales" = "yellow")) +
  labs(x = paste0("PC1 (", var_exp_pct[1], "%)"),
       y = paste0("PC2 (", var_exp_pct[2], "%)"),
       title="PCA of AA alignment, 3pc expansion",
       color = "Order") +
  theme_minimal() +
  theme(legend.position = "right")

with(df[df[["psbA"]]==0, ],
     points(PC1, PC2, col="red"))

numeric_cols <- names(df)[sapply(df, is.numeric)]
numeric_cols <- numeric_cols[!grepl("aa_supermatrix", numeric_cols)]

cors <- sapply(df[ , ..numeric_cols], function(x) cor(x, df$PC2, use="pairwise"))
head(sort(abs(cors), decreasing = TRUE), 20)

cors <- sapply(df[ , ..numeric_cols], function(x) cor(x, df$PC1, use="pairwise"))
head(sort(abs(cors), decreasing = TRUE), 20)

plot(df$psbZ, df$PC2)

sum(df$X > 0)

