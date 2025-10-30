library(data.table)
library(arrow)
library(ggplot2)
aln <- read_parquet("data/tmp/majMinor_aln.pq")
data <- setDT(aln)

stopifnot(is.character(data[[1]]))
for (j in seq_along(data)) set(data, which(is.na(data[[j]])), j, median(data[[j]], na.rm=TRUE))#drop all zero var columns 
invariant_cols <- names(data)[sapply(data, function(x) length(unique(x)) == 1)]
cat("Dropping", length(invariant_cols), "invariant columns\n")
data <- data[, !..invariant_cols]

pca <- prcomp(data[,-1], rank.=1000, scale. = TRUE)
saveRDS(pca, "data/tmp/majMinor_aln_pca.rds")
cat("PCA done!")
quit()
data <- read_parquet("data/processed_data.parquet")
pca <- readRDS("data/tmp/majMinor_aln_pca.rds")
summary(pca)
pvar <- pca$sdev^2 / sum(pca$sdev^2)
pvar <- round(pvar*100,2)

barplot(pvar)
dev.off()
plot(cumsum(pvar)[1:1000],
     main="Cumulative percent variance explained by Principal Components",
     xlab="Number of PCs",
     ylab="Cumulative Percent Variance Explained")
plot(pca$x[,1],pca$x[,2])
scores <- as.data.frame(pca$x)
df <- cbind(data,scores)

plot(pca$x[,2],pca$x[,3])

table(df$Order)

orders <- unique(df$Order)
n_orders <- length(orders)

colors <- c("red", "blue", "seagreen", "orange", "purple", "brown", "black", "darkgreen")
shapes <- c(3,4,11,16, 17, 18, 15, 19)  # circle, triangle, diamond, square, filled circle

# Create color and shape vectors (recycle as needed)
order_colors <- rep(colors, length.out = n_orders)
order_shapes <- rep(shapes, length.out = n_orders)

pairs <- expand.grid(color=colors, shape=shapes)
pairs <- pairs[sample(nrow(pairs), n_orders),]

order_colors <- pairs$color
order_shapes <- pairs$shape

plot_pca_orders <- function(A, B, df, pvar, orders, order_colors, order_shapes, main_title = NULL) {
  stopifnot(A %in% seq_along(pvar), B %in% seq_along(pvar))
  
  # Default title if not provided
  if (is.null(main_title)) {
    main_title <- sprintf("PCA by Order (PC%d vs PC%d)", A, B)
  }
  
  dev.off()  # reset graphics device
  
  # Layout: main plot + legend panel
  layout(matrix(c(1,2), 1, 2), widths = c(3, 2))
  
  # Main plot
  plot(df[[paste0("PC", A)]], df[[paste0("PC", B)]],
       type = "n",
       xlab = sprintf("PC%d (%.1f%%)", A, pvar[A]),
       ylab = sprintf("PC%d (%.1f%%)", B, pvar[B]),
       main = main_title)
  
  # Add points for each order
  for (i in seq_along(orders)) {
    order_subset <- df$Order == orders[i]
    points(df[[paste0("PC", A)]][order_subset],
           df[[paste0("PC", B)]][order_subset],
           col = order_colors[i],
           pch = order_shapes[i],
           cex = 0.8)
  }
  
  # Side panel for legend
  par(mar = c(5, 0, 4, 2))
  plot.new()
  legend("topleft", legend = orders,
         col = order_colors, pch = order_shapes,
         cex = 0.7, pt.cex = 0.7,
         bty = "n", ncol = 2, y.intersp = 0.6)
  
  # Reset layout
  par(mfrow = c(1,1))
}
plot_pca_orders(1, 2, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(3, 4, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(5, 6, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(7, 8, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(9, 10, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(11, 12, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(13, 14, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(15, 16, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(17, 18, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(19, 20, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(21, 22, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2

library(ggplot2)
library(gridExtra)

pairs <- split(1:20, rep(1:10, each=2))
plots <- lapply(seq_along(pairs), function(i){
  A <- pairs[[i]][1]; B <- pairs[[i]][2]
  ggplot(df, aes_string(x=paste0("PC",A), y=paste0("PC",B), color="Order", shape="Order")) +
    geom_point(size=1, alpha=0.7) +
    labs(x=sprintf("PC%d (%.1f%%)",A,pvar[A]),
         y=sprintf("PC%d (%.1f%%)",B,pvar[B]),
         title=sprintf("PC%d vs PC%d",A,B)) +
    theme_bw() +
    theme(legend.position="none",
          plot.title=element_text(size=10),
          axis.title=element_text(size=8),
          axis.text=element_text(size=7))
})

legend_plot <- ggplot(df, aes(PC1, PC2, color=Order, shape=Order)) +
  geom_point() + theme_void() +
  theme(legend.position="bottom", legend.title=element_blank()) +
  guides(color=guide_legend(ncol=4), shape=guide_legend(ncol=4))

g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

legend_grob <- g_legend(legend_plot)
title_grob <- grid::textGrob(sprintf("Cumulative variance explained: %.1f%% (PC1–PC20)", sum(pvar[1:20])),
                             gp=grid::gpar(fontsize=12, fontface="bold"))

gridExtra::grid.arrange(
  do.call(gridExtra::arrangeGrob, c(plots, ncol=5, nrow=2)),
  legend_grob,
  title_grob,
  ncol=1,
  heights=c(10,1,0.5)
)





plot(pca$x[,9],pca$x[,10], col="red")
text(pca$x[,9],pca$x[,10],labels=df$Organism)
head(order(df$PC10))

# flag outliers (top/right extremes)
thr9 <- 100
thr10 <- 100
outliers <- df[df$PC9 > thr9 & df$PC10 > thr10,]
outliers[, "ID"]

embeds <- readRDS("data/tmp/embeds_with_mds.rds")
embed_cols <- grep("embedding", colnames(embeds), value=TRUE)
clean_embeds <- embeds[ManualOutlier == FALSE]

out_ids <- outliers$ID
psa_out_ids <- embeds[ManualOutlier & Gene=="psaA", ID]
embeds[ManualOutlier & Gene=="psaA", MDS_zscore]

plot(embeds[Gene=="psaA", MDS1],embeds[Gene=="psaA", MDS2])
points(embeds[Gene=="psaA" & ID %in% out_ids, MDS1],
       embeds[Gene=="psaA" & ID %in% out_ids, MDS2],
       col="red")

list(
  overlap = intersect(out_ids, psa_out_ids),
  only_in_pc_outliers = setdiff(out_ids, psa_out_ids),
  only_in_psa_outliers = setdiff(psa_out_ids, out_ids)
)

# loadings
vars <- colnames(data)[-1]
loadings <- as.data.table(pca$rotation, keep.rownames="variable")
loadings[, variable := vars]

pc9_top <- loadings[order(-abs(PC9))][1:200, .(variable, PC9)]
pc10_top <- loadings[order(-abs(PC10))][1:200, .(variable, PC10)]

list(pc9_top=pc9_top, pc10_top=pc10_top)

table(data$cds_supermatrix_14817)

with(df[grep("Poales", df$Taxonomy), ],
     points(PC2, PC3, col="seagreen"))

with(df[grep("Fabales", df$Taxonomy), ],
     points(PC2, PC3, col="yellow"))

with(df[grep("Pinales", df$Taxonomy), ],
     points(PC2, PC3, col="blue"))




# Create the plot
plot(df$PC2, df$PC3, type="n", 
     xlab="PC2", ylab="PC3", 
     main="PCA by Order")

# Plot each order
for(i in 1:n_orders) {
  order_subset <- df$Order == orders[i]
  points(df$PC2[order_subset], df$PC3[order_subset], 
         col = order_colors[i], 
         pch = order_shapes[i],
         cex = 0.8)
}

legend("topright", legend = orders, 
       col = order_colors, pch = order_shapes, 
       cex = 0.6, ncol = 2)

dev.off()
# Set up layout with main plot and side panel for legend
layout(matrix(c(1,2), 1, 2), widths = c(3, 2))
# Main plot
plot(df$PC2, df$PC3, type="n", 
     xlab=paste0("PC2 ",pvar[2],"%"),
     ylab=paste0("PC3 ",pvar[3],"%"), 
     main="PCA by Order (onehot aln encoding)")

# Plot each order
for(i in 1:n_orders) {
  order_subset <- df$Order == orders[i]
  points(df$PC2[order_subset], df$PC3[order_subset], 
         col = order_colors[i], 
         pch = order_shapes[i],
         cex = 0.8)
}

# Side panel for legend
par(mar = c(5, 0, 4, 2))
plot.new()
legend("left", legend = orders,
       col = order_colors, pch = order_shapes,
       cex = 0.7, pt.cex = 0.7,
       bty = "n", ncol = 2, y.intersp = 0.6)
par(mfrow = c(1,1))

plot_pca_orders <- function(A, B, df, pvar, orders, order_colors, order_shapes, main_title = NULL) {
  stopifnot(A %in% seq_along(pvar), B %in% seq_along(pvar))
  
  # Default title if not provided
  if (is.null(main_title)) {
    main_title <- sprintf("PCA by Order (PC%d vs PC%d)", A, B)
  }
  
  dev.off()  # reset graphics device
  
  # Layout: main plot + legend panel
  layout(matrix(c(1,2), 1, 2), widths = c(3, 2))
  
  # Main plot
  plot(df[[paste0("PC", A)]], df[[paste0("PC", B)]],
       type = "n",
       xlab = sprintf("PC%d (%.1f%%)", A, pvar[A]),
       ylab = sprintf("PC%d (%.1f%%)", B, pvar[B]),
       main = main_title)
  
  # Add points for each order
  for (i in seq_along(orders)) {
    order_subset <- df$Order == orders[i]
    points(df[[paste0("PC", A)]][order_subset],
           df[[paste0("PC", B)]][order_subset],
           col = order_colors[i],
           pch = order_shapes[i],
           cex = 0.8)
  }
  
  # Side panel for legend
  par(mar = c(5, 0, 4, 2))
  plot.new()
  legend("topleft", legend = orders,
         col = order_colors, pch = order_shapes,
         cex = 0.7, pt.cex = 0.7,
         bty = "n", ncol = 2, y.intersp = 0.6)
  
  # Reset layout
  par(mfrow = c(1,1))
}
plot_pca_orders(1, 2, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(2, 3, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(3, 4, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(5, 6, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(7, 8, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(9, 10, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(11, 12, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2

plot_pca_orders(13, 12, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(20, 21, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2





sort(abs(pca$rotation[,1]),decreasing=TRUE)[1:10]

map <- read.csv("data/tmp/cds_supermatrix.fasta_mapping.csv")
map$site <- seq_along(1:nrow(map))
head(map)
loadings <- as.data.frame(pca$rotation)
site_nums <- as.numeric(gsub("cds_supermatrix_([0-9]+)_[ATCG]", "\\1", rownames(loadings)))
gene_mapping <- map$gene[match(site_nums, map$site)]
gene_contrib <- data.frame(
  site = site_nums,
  gene = gene_mapping,
  loading = abs(loadings$PC1)
)


top_genes <- aggregate(loading ~ gene, gene_contrib, sum)
top_genes[order(-top_genes$loading),]


# Get top contributing genes programmatically
top_genes <- top_genes[order(-top_genes$loading),]
top_gene_names <- head(top_genes$gene, nrow(top_genes))  # Top 10 genes

# Calculate gene presence stats
gene_stats <- data.frame(
 gene = top_gene_names,
 pc1_loading = top_genes$loading,
 prop_present = sapply(df[,..top_gene_names], mean),
 variance = sapply(df[,..top_gene_names], var)
)

ggplot(gene_stats, aes(x = prop_present, y = pc1_loading)) +
  geom_point(aes(size = variance, color = variance), alpha = 0.7) +
  geom_text(aes(label = gene), vjust = -0.5, hjust = 0.5, size = 3) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    x = "Proportion Present (0 = always absent, 1 = always present)",
    y = "PC1 Loading (Absolute Sum)", 
    title = "PC1 Explained by Gene Presence Variation",
    subtitle = "Size and color = variance in presence/absence"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Quick correlation check
cor(gene_stats$variance, gene_stats$pc1_loading)
cor(gene_stats$prop_present, gene_stats$pc1_loading)

hist(df$pheno_Topt_site_p50)
hist(df$pheno_wc)


pca <- readRDS("data/tmp/majMinor_aln_pca.rds")
summary(pca)
pvar <- pca$sdev^2 / sum(pca$sdev^2)
pvar <- round(pvar*100, 3)
plot(pvar)
plot(pca$x[,1],pca$x[,2])
scores <- as.data.frame(pca$x)
df <- cbind(data,scores)

plot_pca_orders(1, 4, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(6, 5, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(8, 7, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(8, 9, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2
plot_pca_orders(30, 31, df, pvar, orders, order_colors, order_shapes)   # PC1 vs PC2


with(df[grep("Poaceae", Taxonomy), ],
     points(PC1, PC2, col="seagreen"))

with(df[grep("Fabales", Taxonomy), ],
     points(PC1, PC2, col="yellow"))

with(df[grep("Pinales", Taxonomy), ],
     points(PC1, PC2, col="blue"))



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

