#gc 

library(data.table)
library(stringr)

par(mfrow=c(1,1))
gc_data <- fread("/workdir/hdd29/chloroplast_genome_evaluation/data/chloroplast_genes_gc_content.csv")
gc_data$LUI <- str_split_i(gc_data$sample,"_",1)


hist(gc_data[which(gc_data$LUI=="NC0016662Zeamays")]$gc_content, main= "Zea mays plastid Genes GC")
hist(gc_data[which(gc_data$LUI=="NC0313331Oryzasativa")]$gc_content, main= "Oryzasativa Genes GC")

plot(gc_data[which(gc_data$LUI=="NC0016662Zeamays")]$gc_content, gc_data[which(gc_data$LUI=="NC0313331Oryzasativa")]$gc_content)

genes <- unique(gc_data$gene)

for (gene_name in genes) {
  sset <- gc_data[which(gc_data$gene==gene_name),]
  hist(sset$gc_content, main= paste0(gene_name, " GC, n=", nrow(sset)))
}

library(dplyr)

gene_gc <- gc_data %>%
  group_by(gene) %>%
  summarise(
    n=n(),
    median_gc = median(gc_content),
    gc_sd = sd(gc_content)
  )

hist(gene_gc$median_gc)
gene_gc[which(str_detect(gene_gc$gene, "rn")),]


h1 <- hist(gene_gc[which(str_detect(gene_gc$gene, "trn")),]$median_gc, plot = FALSE)
h2 <- hist(gene_gc[which(!str_detect(gene_gc$gene, "trn")),]$median_gc, plot = FALSE)
xlim <- range(c(h1$breaks, h2$breaks))
ylim <- range(c(h1$counts, h2$counts))
plot(h1, col = rgb(0.2, 0.4, 0.8, 0.5), xlim = xlim, ylim = ylim,
     main = "tRNA & non tRNA contributions to GC", xlab = "Value", ylab = "Frequency")
plot(h2, col = rgb(0.2, 0.8, 0.4, 0.5), add = TRUE)
legend("topright", legend = c("trn", "others"),
       fill = c(rgb(0.2, 0.4, 0.8, 0.5), rgb(0.2, 0.8, 0.4, 0.5)))

h1 <- hist(gene_gc[which(str_detect(gene_gc$gene, "rn")),]$median_gc, plot = FALSE)
h2 <- hist(gene_gc[which(!str_detect(gene_gc$gene, "rn")),]$median_gc, plot = FALSE)
xlim <- range(c(h1$breaks, h2$breaks))
ylim <- range(c(h1$counts, h2$counts))
plot(h1, col = rgb(0.2, 0.4, 0.8, 0.5), xlim = xlim, ylim = ylim,
     main = "RNA gene & non RNA contributions to GC", xlab = "Value", ylab = "Frequency")
plot(h2, col = rgb(0.2, 0.8, 0.4, 0.5), add = TRUE)
legend("topright", legend = c("rn", "others"),
       fill = c(rgb(0.2, 0.4, 0.8, 0.5), rgb(0.2, 0.8, 0.4, 0.5)))


hist(gene_gc$gc_sd)
plot(gene_gc$median_gc, gene_gc$gc_sd)

cor(gene_gc)
