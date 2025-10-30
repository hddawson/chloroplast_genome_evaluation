library(arrow)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(Biostrings)


df <- read_parquet("data/all_residue_embeddings.parquet")

mat <- as.matrix(df[,grep("embedding", colnames(df))])
cat("made mat!")
#pca <- prcomp(mat,rank.=10, scale.=T, center=T)

#saveRDS(pca, "data/tmp/residuespca.rds")
pca <- readRDS("data/tmp/residuespca.rds")

data <- read_parquet("data/processed_data.parquet")

# specify phenotype column
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln <- read_parquet("data/tmp/majMinor_aln.pq")
pcs_IDS <- aln$index
scores <- as.data.frame(ev_pcs$x)
scores <- cbind(ID = pcs_IDS, scores)
n_pcs <- 100
pc_names <- paste0("PC", seq_len(n_pcs))
colnames(scores)[-1] <- paste0("PC", seq_len(ncol(scores)-1))
scores <- scores %>% select(ID, all_of(pc_names))

fasta <- "data/tmp/alignedGenes/psbA_AA_aligned.fasta"
aln <- readAAStringSet(fasta)

# Get sequence names (should match your IDs or Gene+ID pattern)
names(aln) <- sub("\\|.*", "", names(aln))

aln_mat <- as.matrix(aln)
seq_ids <- rownames(aln_mat)
aligned_len <- ncol(aln_mat)

colnames(aln_mat) <- seq_len(ncol(aln_mat))
aln_df <- as.data.frame(aln_mat) %>%
  tibble::rownames_to_column("ID") %>%
  tidyr::pivot_longer(-ID, names_to = "Aligned_Position", values_to = "Residue") %>%
  mutate(Aligned_Position = as.integer(Aligned_Position)) %>%
  group_by(ID) %>%
  mutate(Residue_Index = ifelse(Residue == "-", NA_integer_, cumsum(Residue != "-"))) %>%
  ungroup() %>%
  filter(!is.na(Residue_Index)) %>%
  select(ID, Residue_Index, Aligned_Position)



common_ids <- intersect(df$ID, pcs_IDS)
df <- df %>% filter(ID %in% common_ids)
scores <- scores %>% filter(ID %in% common_ids)

df_joined <- df %>%
  inner_join(data %>% select(ID, !!sym(pheno_col)), by="ID") %>%
  rename(pheno = !!sym(pheno_col)) %>%
  left_join(scores, by="ID")

df_joined <- df_joined %>%
  mutate(ID = sub("\\|.*", "", ID)) %>%
  left_join(aln_df, by = c("ID", "Residue_Index"))

emb_cols <- grep("^embedding_", colnames(df_joined), value = TRUE)
df_joined <- df_joined %>% mutate(GroupID = paste(Gene, Residue_Index, sep = "_"))




groups <- unique(df_joined$GroupID)
results_list <- vector("list", length(groups))


require(doMC)
registerDoMC(cores = 10)

# loop over residues
for (i in seq_along(groups)) {
  gid <- groups[i]
  sub <- df_joined %>% filter(GroupID == gid)
  if (nrow(sub) < 1000) next
  y <- sub$pheno
  X_emb <- scale(as.matrix(sub[, emb_cols]))
  X_pcs <- scale(as.matrix(sub[, pc_names]))
  X_all <- cbind(X_emb, X_pcs)
  penalty.factor <- c(rep(1, ncol(X_emb)), rep(0, ncol(X_pcs)))
  
  fit <- tryCatch(cv.glmnet(X_all, y, alpha=0, penalty.factor=penalty.factor, parallel = T), error=function(e) NULL)
  if (is.null(fit)) next
  best <- glmnet(X_all, y, alpha=0, lambda=fit$lambda.min, penalty.factor=penalty.factor)
  coefs <- as.data.frame(as.matrix(coef(best)))
  coefs$Predictor <- rownames(coefs)
  colnames(coefs)[1] <- "Estimate"
  coefs <- coefs %>%
    filter(Predictor != "(Intercept)", Estimate != 0) %>%
    mutate(Gene=sub$Gene[1],
           Residue_Index=sub$Residue_Index[1],
           Residue=sub$Residue[1],
           Lambda=fit$lambda.min,
           N=nrow(sub)) %>%
    select(Gene, Residue_Index, Residue, Predictor, Estimate, Lambda, N)
  
  yhat <- predict(best, X_all)
  r2 <- cor(y, yhat)^2
  coefs <- coefs %>% mutate(R2 = r2)
  results_list[[i]] <- coefs
  if (i %% 50 == 0) message("Processed ", i, "/", length(groups), " residues")
}
saveRDS(results_list, "results/results_list.rds")
results_df <- do.call(rbind, results_list)
saveRDS(results_df, "results/residue_predictor_coefs.rds")

results_df <- readRDS("results/residue_predictor_coefs.rds")

r1 <- df_joined[df_joined$Residue_Index==234,]
table(df_joined[df_joined$Residue_Index==234,"Residue"])
table(df_joined[df_joined$Residue_Index==236,"Residue"])
table(df_joined[df_joined$Residue_Index==238,"Residue"])

df_plot <- results_df %>%
  filter(grepl("embedding", Predictor)) %>%
  mutate(Predictor=factor(Predictor, levels=unique(Predictor)),
         Residue_Index=as.numeric(Residue_Index))

ggplot(df_plot, aes(x=Residue_Index, y=Predictor, fill=Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Residue Index", y="Predictor", fill="Estimate")

ggplot(df_plot, aes(x=Residue_Index, y=Predictor, fill=Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low="#313695", mid="white", high="#A50026",
                       midpoint=0, limits=c(-max(abs(df_plot$Estimate)), max(abs(df_plot$Estimate))),
                       oob=scales::squish) +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x="Residue Index", y=NULL, fill="Estimate")

ggplot(df_plot %>% filter(Residue_Index >= 200, Residue_Index <= 300),
       aes(x=Residue_Index, y=Predictor, fill=Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low="#313695", mid="white", high="#A50026",
                       midpoint=0, limits=c(-max(abs(df_plot$Estimate)), max(abs(df_plot$Estimate))),
                       oob=scales::squish) +
  facet_wrap(~Gene, scales="free_x") +
  scale_x_continuous(breaks=seq(200, 300, by=1)) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x="Residue Index", y=NULL, fill="Estimate")

df_bar <- df_plot %>%
  group_by(Gene, Residue_Index) %>%
  summarise(TotalAbsEffect=sum(abs(Estimate), na.rm=TRUE), .groups="drop")

ggplot(df_bar, aes(x=Residue_Index, y=TotalAbsEffect)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Residue Index", y="Cumulative |Effect Size|") +
  scale_x_continuous(breaks=seq(min(df_bar$Residue_Index), max(df_bar$Residue_Index), by=10))


mat <- df_plot %>%
  filter(grepl("embedding", Predictor)) %>%
  select(Gene, Residue_Index, Predictor, Estimate) %>%
  pivot_wider(names_from=Predictor, values_from=Estimate, values_fill=0)

residue_meta <- mat %>% select(Gene, Residue_Index)
mat_num <- as.matrix(mat %>% select(-Gene, -Residue_Index))
pca <- prcomp(mat_num, scale.=TRUE, center=TRUE)

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
xlab <- sprintf("PC1 (%.1f%%)", var_expl[1]*100)
ylab <- sprintf("PC2 (%.1f%%)", var_expl[2]*100)

plot(pca$x[,1], pca$x[,2], col="white", xlab=xlab, ylab=ylab, main="Residue PCA")
text(pca$x[,1], pca$x[,2], labels=residue_meta$Residue_Index, col="red", cex=1)

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
xlab <- sprintf("PC3 (%.1f%%)", var_expl[3]*100)
ylab <- sprintf("PC4 (%.1f%%)", var_expl[4]*100)

plot(pca$x[,3], pca$x[,4], col="white", xlab=xlab, ylab=ylab, main="Residue PCA")
text(pca$x[,3], pca$x[,4], labels=residue_meta$Residue_Index, col="red", cex=1)

scores <- cbind(residue_meta, pca$x)


table(results_df[results_df$Residue_Index==234, "Residue"])
a <- results_df[results_df$Residue_Index==234,]

plot(-log10(results$Pvalue), main="significance of residue-wise bio8 ~ emb_1 + emb_2 + ... + emb_960",
     ylab="-log10(p)", xlab="Residue index Unaligned")
segments(x0=220,y0=280, x1=275, y1=280, col="red")
text(240, 282,"D-E loop", col="red")

plot(results$R2)
