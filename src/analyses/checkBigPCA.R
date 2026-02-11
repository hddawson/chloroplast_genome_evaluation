library(arrow)
library(data.table)
library(ggplot2)

pca <- read_parquet("results/embeddings_full_pca.parquet")
library(dplyr)
library(Peptides)

res_keep <- aaList()  # common residues only
data <- read_parquet("data/processed_data.parquet")

zmID <- data$ID[grep("Zea", data$Organism)]

zm_dat <- pca %>% filter(ID == zmID) 
zm_dat %>% 
  ggplot(aes(PC1, PC2)) +
  stat_binhex(bins = 120) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() + ggtitle("Zm residues in ESMC space")

tdID <- data$ID[grep("rabidopsis th", data$Organism)]

at_dat <- pca %>% filter(ID == tdID) 
at_dat %>% 
  ggplot(aes(PC1, PC2)) +
  stat_binhex(bins = 120) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() + ggtitle("At residues in ESMC space")

pinusID <- data$ID[grep("Pinus", data$Organism)]

pca %>%
  filter(ID %in% pinusID) %>%
  ggplot(aes(PC1, PC2)) +
  stat_binhex(bins = 120) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() + ggtitle("Pinus residues in ESMC space (n=26)")

genes <- unique(pca$Gene)

data(AAdata)

AAdata$kideraFactors
table(pca$Residue[1:10000])

library(ggplot2)
library(dplyr)
library(viridis)

# turn Kidera factors into a long lookup
kf_lookup <- bind_rows(
  lapply(names(AAdata$kideraFactors), function(kf) {
    tibble(
      Residue = names(AAdata$kideraFactors[[kf]]),
      value   = as.numeric(AAdata$kideraFactors[[kf]]),
      KF      = kf
    )
  })
)


kf_vec <- AAdata$kideraFactors[["KF1"]]

ggplot(pca,
       aes(PC3, PC4)) +
  stat_summary_hex(
    aes(z = kf_vec[Residue], fill = after_stat(value)),
    fun = mean,
    bins = 120
  ) +
  scale_fill_viridis_c()

plots <- lapply(names(AAdata$kideraFactors), function(kf) {
  kf_vec <- AAdata$kideraFactors[[kf]]
  
  ggplot(pca,
         aes(PC3, PC4)) +
    stat_summary_hex(
      aes(z = kf_vec[Residue], fill = after_stat(value)),
      fun = mean,
      bins = 120
    ) +
    scale_fill_viridis_c() +
    labs(title = kf)
})
#gridrange them '


pdf("kidera_hexplots_pc3_pc4.pdf", width = 14, height = 10)
do.call(grid.arrange, plots)
dev.off()


kf_vec <- AAdata$kideraFactors[[kf]]  # named by residue

ggplot(pca %>% filter(Gene %in% genes),
       aes(PC4, PC4)) +
  stat_summary_hex(
    aes(
      z = kf_vec[Residue],
      fill = after_stat(value)
    ),
    fun = mean,
    bins = 120
  ) +
  facet_wrap(~Residue) +
  scale_fill_viridis_c()


ggplot(pca %>% filter(Gene %in% genes),
       aes(PC4, PC4)) +
  stat_summary_hex(
    aes(
      z = kf_vec[Residue],
      fill = after_stat(value)
    ),
    fun = mean,
    bins = 120
  ) +
  facet_wrap(~Residue) +
  scale_fill_viridis_c()

plots <- lapply(unique(base$KF), function(kf) {
  base %>%
    filter(KF == kf) %>%
    ggplot(aes(PC3, PC4)) +
    stat_summary_hex(
      aes(z = value, fill = after_stat(value)),
      fun = mean,
      bins = 120
    ) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = kf, fill = "Mean Kidera")
})

# e.g. print KF1
plots[[1]]

library(gridExtra)
do.call(grid.arrange, plots)

pca %>%
  ggplot(aes(PC4, PC4)) +
  stat_binhex(bins = 120) +
  facet_wrap(~Residue) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal()

pca %>%
  mutate(pos_bin = cut_number(Residue_Index, 20)) %>%  # ~5% gene chunks
  ggplot(aes(PC1, PC2)) +
  stat_binhex(bins = 150) +
  facet_wrap(~pos_bin) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal()

barplot(c(0.28295055, 0.08695123, 0.048718642, 0.028368318, 0.022306975), main="Variance Explained")

library(ape)
tree <- read.tree("results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree")
data <- as.data.table(read_parquet("data/processed_data.parquet"))

# Drop gymnosperms
pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree$tip.label)
tree <- root(tree, outgroup = pinales_in_tree)
is.rooted(tree)

pc1_tip <- pca %>% group_by(ID) %>% summarise(PC1 = mean(PC1, na.rm=TRUE))
tree$tip.label <- intersect(tree$tip.label, pc1_tip$ID)
library(ape); library(dplyr); library(phytools)

plot(contMap(tree, setNames(pc1_tip$PC1, pc1_tip$ID), plot=FALSE), ftype="off")

head(data$pheno_wc2.1_2.5m_bio_8_p50)
head(data$ID)

pca %>%
  left_join(data[, c("ID","pheno_wc2.1_2.5m_bio_8_p50")], by="ID") %>%
  ggplot(aes(PC3, PC5, z = pheno_wc2.1_2.5m_bio_8_p50)) +
  stat_summary_hex(fun = mean, bins = 150) +
  scale_fill_viridis_c(name = "Temp (BIO8 p50)") +
  theme_minimal() +
  ggtitle("PCA colored by temperature")
pc_cols <- grep("^PC[0-9]+$", names(pca), value = TRUE)

pc_cor <- pca %>%
  left_join(data[, c("ID","pheno_wc2.1_2.5m_bio_8_p50")], by="ID") %>%
  summarise(across(all_of(pc_cols),
                   ~ cor(.x, pheno_wc2.1_2.5m_bio_8_p50,
                         use="complete.obs", method="spearman")))
