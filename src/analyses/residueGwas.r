library(arrow)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(Biostrings)

#r
# Load data
df <- read_parquet("data/all_residue_embeddings.parquet")
data <- read_parquet("data/processed_data.parquet")
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

# Load PCs
ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln <- read_parquet("data/tmp/majMinor_aln.pq")
pcs_IDS <- aln$index
scores <- as.data.frame(ev_pcs$x)
scores <- cbind(ID = pcs_IDS, scores)
n_pcs <- 100
pc_names <- paste0("PC", seq_len(n_pcs))
colnames(scores)[-1] <- paste0("PC", seq_len(ncol(scores)-1))
scores <- scores %>% select(ID, all_of(pc_names))

# Load alignment and create lookup
fasta <- "data/tmp/alignedGenes/psbA_AA_aligned.fasta"
aln <- readAAStringSet(fasta)
names(aln) <- sub("\\|.*", "", names(aln))
aln_mat <- as.matrix(aln)
colnames(aln_mat) <- seq_len(ncol(aln_mat))


aln_lookup <- as.data.frame(aln_mat) %>%
  tibble::rownames_to_column("ID") %>%
  tidyr::pivot_longer(-ID, names_to = "Aligned_Position", values_to = "Residue") %>%
  mutate(Aligned_Position = as.integer(gsub("V", "", Aligned_Position))) %>%
  group_by(ID) %>%
  mutate(Residue_Index = ifelse(Residue == "-", NA_integer_, cumsum(Residue != "-"))) %>%
  ungroup() %>%
  filter(!is.na(Residue_Index)) %>%
  select(ID, Residue_Index, Aligned_Position)

head(aln_lookup)# Filter to common IDs
common_ids <- intersect(df$ID, pcs_IDS)
stopifnot(length(common_ids) > 0)

# Single join pipeline
df_joined <- df %>%
  mutate(ID = sub("\\|.*", "", ID)) %>%
  inner_join(aln_lookup, by = c("ID", "Residue_Index")) %>%
  inner_join(data %>% select(ID, pheno = !!pheno_col), by = "ID") %>%
  inner_join(scores, by = "ID") %>%
  filter(ID %in% common_ids)

stopifnot(nrow(df_joined) > 0)
stopifnot(!all(is.na(df_joined$Aligned_Position)))

# Group by aligned position
df_joined <- df_joined %>% 
  mutate(GroupID = paste(Gene, Aligned_Position, sep = "_"))

groups <- unique(df_joined$GroupID)

require(doMC)
registerDoMC(cores = 10)

results_list <- vector("list", length(groups))
emb_cols <- grep("^embedding_", colnames(df_joined), value = TRUE)

# Loop over aligned positions
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
           Aligned_Position=sub$Aligned_Position[1],
           Residue=sub$Residue[1],
           Lambda=fit$lambda.min,
           N=nrow(sub)) %>%
    select(Gene, Aligned_Position, Residue, Predictor, Estimate, Lambda, N)

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

a <- results_df %>%
  select(Aligned_Position, R2) %>%
  unique()
lines(a$Aligned_Position,a$R2)



df_plot <- results_df %>%
  filter(grepl("embedding", Predictor)) %>%
  mutate(Predictor=factor(Predictor, levels=unique(Predictor)),
         Aligned_Position=as.numeric(Aligned_Position))

ggplot(df_plot, aes(x=Aligned_Position, y=Predictor, fill=Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Aligned_Position", y="Predictor", fill="Estimate")

ggplot(df_plot, aes(x=Aligned_Position, y=Predictor, fill=Estimate)) +
  geom_tile() +
  scale_fill_gradient2(low="#313695", mid="white", high="#A50026",
                       midpoint=0, limits=c(-max(abs(df_plot$Estimate)), max(abs(df_plot$Estimate))),
                       oob=scales::squish) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x="Aligned_Position", y=NULL, fill="Estimate")

ggplot(df_plot %>% filter(Aligned_Position >= 200, Aligned_Position <= 300),
       aes(x=Aligned_Position, y=Predictor, fill=Estimate)) +
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
  labs(x="Aligned_Position", y=NULL, fill="Estimate")

df_bar <- df_plot %>%
  group_by(Gene, Aligned_Position) %>%
  summarise(TotalAbsEffect=sum(abs(Estimate), na.rm=TRUE), .groups="drop")

ggplot(df_bar, aes(x=Aligned_Position, y=TotalAbsEffect)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Aligned_Position", y="Cumulative |Effect Size|") +
  scale_x_continuous(breaks=seq(min(df_bar$Aligned_Position), max(df_bar$Aligned_Position), by=10))

data[grep("Arabidopsis thaliana", data$Organism), "ID"] #AP000423.1
at_ID <- "AP000423.1"

at_df <- df_joined[df_joined$ID==at_ID,] %>% select(Residue_Index, Aligned_Position)
str(at_df)
library(jsonlite)

feat_json <- 'data/psbA_features.json'
j <- fromJSON(feat_json) 
feats <- as.data.frame(j) 


#a$location$end$value
#a$ligand$name

lines(df_bar$TotalAbsEffect)

# assume df_bar$TotalAbsEffect is already plotted
plot(df_bar$TotalAbsEffect, type = "l", main = "Per-residue effect", xlab = "Residue Index", ylab = "Total Abs Effect")

# extract positions and labels
pos <- a$location$end$value
labs <- a$ligand$name

# draw vertical lines at each residue position
abline(v = pos, col = "peachpuff", lty = 2, lwd = 1.5)

# add text labels slightly above the plot
text(x = pos, 
     y = par("usr")[4] * 0.7,  # near top of plot
     labels = labs, 
     srt = 90, cex = 0.7, col = "brown3", adj = 0)

feats$loc <- feats$features.location$start$value
feats$labl <- feats$features.ligand$name
feat_df <- feats %>%
  left_join(at_df, by=c("loc"="Residue_Index")) %>%
  filter(!is.na(Aligned_Position))

p <- ggplot(df_bar, aes(x=Aligned_Position, y=TotalAbsEffect)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Aligned_Position", y="Cumulative |Effect Size|")

p + geom_vline(data=feat_df, aes(xintercept=Aligned_Position),
               color="peachpuff", linetype="dashed", linewidth=0.7) +
  geom_text(data=feat_df,
            aes(x=Aligned_Position, y=max(df_bar$TotalAbsEffect)*0.7,
                label=labl),
            angle=90, vjust=0, hjust=0, size=5, color="black")


a <- results_df %>%
  select(Aligned_Position, R2) %>%
  distinct()

feat_df <- feats %>%
  transmute(Residue_Index = features.location$start$value,
            Label = features.ligand$name) %>%
  left_join(at_df, by="Residue_Index") %>%
  filter(!is.na(Aligned_Position))

ggplot(a, aes(x=Aligned_Position, y=R2)) +
  geom_line(color="#2C7BB6", linewidth=0.8) +
  geom_vline(data=feat_df, aes(xintercept=Aligned_Position),
             color="orange", linetype="dashed", linewidth=0.6) +
  geom_text(data=feat_df,
            aes(x=Aligned_Position, y=0.38, label=Label),
            angle=90, vjust=0, hjust=0, size=3, color="black") +
  theme_minimal(base_size=12) +
  labs(x="Aligned Position", y=expression(R^2))

a$dist_to_bind <- sapply(a$Aligned_Position, function(x)
  min(abs(x - feat_df$Aligned_Position)))

ggplot(a, aes(x=dist_to_bind, y=R2)) +
  geom_point(alpha=0.6, color="#2C7BB6") +
  geom_smooth(method="loess", se=FALSE, color="black") +
  theme_minimal(base_size=12) +
  labs(x="Distance to Nearest Binding Site (aligned residues)", y=expression(R^2))

LSD::heatscatter(as.vector(a$dist_to_bind),as.vector(a$R2))
cor(as.vector(a$dist_to_bind),as.vector(a$R2), method="pearson")
