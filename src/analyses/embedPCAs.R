library(data.table)
library(arrow)
library(ggplot2)
library(patchwork)
library(viridisLite)


X_dt <- readRDS("data/tmp/x_dt_imputed.rds")
df <- read_parquet("data/processed_data.parquet")
setDT(df)

pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
embed_pred_cols <- setdiff(colnames(X_dt), c("ID","Taxonomy"))

aln <- read_parquet("data/tmp/majMinor_aln.pq")
setDT(aln)
for (j in seq_along(aln)) set(aln, which(is.na(aln[[j]])), j, median(aln[[j]], na.rm=TRUE))
#drop all zero var columns 
invariant_cols <- names(aln)[sapply(aln, function(x) length(unique(x)) == 1)]
aln <- aln[, !..invariant_cols]
aln_cols <- grep("cds", colnames(aln))

pheno <- df[[pheno_col]]

aln_cors <- sapply(aln_cols, function(idx)
  cor(aln[[idx]], pheno, use="pairwise.complete.obs"))

sum(X_dt$ID==df$ID)

setDT(df)
X_dt <- merge(
  X_dt,
  df[, .(ID, pheno = get(pheno_col))],
  by = "ID"
)

e_pheno <- X_dt$pheno
embed_cors <- sapply(embed_pred_cols, function(cn)
  cor(X_dt[[cn]], e_pheno, use="pairwise.complete.obs"))


hist(aln_cors)
hist(embed_cors)

pdf("figures/embed_aln_cors.pdf", width=6, height=4)
hist(embed_cors, breaks=50, col=rgb(0,0,1,0.4), border=NA,
     xlab="Correlation with phenotype", main="Embedding vs Alignment correlations")
hist(aln_cors, breaks=50, col=rgb(1,0,0,0.4), border=NA, add=TRUE)
legend("topright", legend=c("Embedding","Alignment"),
       fill=c(rgb(0,0,1,0.4),rgb(1,0,0,0.4)), border=NA)
dev.off()

for (c in pred_cols) if (any(is.na(X_dt[[c]]))) X_dt[is.na(get(c)), (c) := median(X_dt[[c]], na.rm=TRUE)]
y_all <- df[match(X_dt$ID, df$ID), get(pheno_col)]
genes <- unique(vapply(strsplit(pred_cols, "__"), `[`, 1, FUN.VALUE=character(1)))
dir.create("plots/genePcas", showWarnings=FALSE)





for (gene in genes) {
  cols <- grep(paste0("^", gene, "__"), pred_cols, value=TRUE)
  if (length(cols) < 2) next
  mat <- as.matrix(X_dt[, ..cols])
  ok <- complete.cases(mat) & !is.na(y_all)
  if (sum(ok) < 10) next
  mat2 <- mat[ok, , drop=FALSE]
  y <- y_all[ok]
  pr <- prcomp(mat2, center=TRUE, scale.=TRUE)
  
  cor_pcs <- sapply(seq_len(ncol(pr$x)), function(i) cor(pr$x[,i], y))
  cor_embeds <- sapply(seq_len(ncol(mat2)), function(i) cor(mat2[,i], y))
  df_comp <- rbind(
    data.frame(Type="PC", R=cor_pcs),
    data.frame(Type="Embedding", R=cor_embeds)
  )
  
  p_comp <- ggplot(df_comp, aes(x=R, fill=Type)) +
    geom_histogram(alpha=0.6, position="identity", bins=30) +
    scale_fill_viridis_d() +
    labs(title=paste(gene,"Temp correlations"), x="Pearson r", y="Count") +
    theme_minimal()
  
  ggsave(sprintf("plots/genePcas/%s_pc_vs_embed_corrhist.png", gene), p_comp, width=6, height=4, dpi=300)
}

for (gene in genes) {
  cols <- grep(paste0("^", gene, "__"), pred_cols, value=TRUE)
  if (length(cols) < 2) next
  mat <- as.matrix(X_dt[, ..cols])
  ok <- complete.cases(mat) & !is.na(y_all)
  if (sum(ok) < 10) next
  mat2 <- mat[ok, , drop=FALSE]
  y <- y_all[ok]
  pr <- prcomp(mat2, center=TRUE, scale.=TRUE)
  varp <- pr$sdev^2 / sum(pr$sdev^2) * 100
  kshow <- min(length(varp), 9)
  df_var <- data.frame(PC = paste0("PC", seq_along(varp))[1:kshow], Var = varp[1:kshow])
  pcs <- as.data.frame(pr$x)
  pcs$Temp <- y
  
  p_var <- ggplot(df_var, aes(PC, Var)) + geom_col() + labs(y="% variance", x=NULL, title=paste(gene, "variance")) + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1))
  
  p12 <- ggplot(pcs, aes(PC1, PC2, color=Temp)) + geom_point(size=1.2) + labs(x=paste0("PC1 (", round(varp[1],1),"%)"), y=paste0("PC2 (", round(varp[2],1),"%)"), title="PC1 vs PC2") + scale_color_viridis_c() + theme_minimal()
  if (ncol(pcs) >= 4) p34 <- ggplot(pcs, aes(PC3, PC4, color=Temp)) + geom_point(size=1.2) + labs(x=paste0("PC3 (", round(varp[3],1),"%)"), y=paste0("PC4 (", round(varp[4],1),"%)"), title="PC3 vs PC4") + scale_color_viridis_c() + theme_minimal() else p34 <- ggplot() + geom_blank() + labs(title="PC3 vs PC4 (NA)")
  
  cors <- sapply(seq_len(ncol(pr$x)), function(i) cor.test(pr$x[,i], y, method="pearson")$estimate)
  ps <- sapply(seq_len(ncol(pr$x)), function(i) cor.test(pr$x[,i], y, method="pearson")$p.value)
  df_cor <- data.frame(PC = paste0("PC", seq_along(cors)), R = cors, P = ps)
  df_cor$sig <- ifelse(df_cor$P < .001, "***", ifelse(df_cor$P < .01, "**", ifelse(df_cor$P < .05, "*", "")))
  p_cor <- ggplot(df_cor, aes(PC, R, fill = -log10(P))) + geom_col() + geom_text(aes(label=sig), vjust=-0.5) + labs(y="Pearson r", x=NULL, title="PC vs Temperature correlations") + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1))
  
  out <- (p_var | (p12 / p34)) / p_cor + plot_layout(heights = c(1,1))
  ggsave(sprintf("plots/genePcas/%s_pca_summary.png", gene), out, width=10, height=9, dpi=300)
}
