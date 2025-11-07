library(arrow)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(Biostrings)
library(tibble)

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
n_pcs <- 1000
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

results_list <- vector("list", length(groups))
emb_cols <- grep("^embedding_", colnames(df_joined), value = TRUE)

# Loop over aligned positions
# Loop over aligned positions
for (i in seq_along(groups)) {
  gid <- groups[i]
  sub <- df_joined %>% filter(GroupID == gid)
  if (nrow(sub) < 1000) next
  
  y <- sub$pheno
  X_emb <- scale(as.matrix(sub[, emb_cols]))
  X_pcs <- scale(as.matrix(sub[, pc_names]))
  
  #-----------------------------------
  # Base (reduced) model: PCs only
  #-----------------------------------
  df_reduced <- as.data.frame(cbind(y = y, X_pcs))
  fit_reduced <- tryCatch(
    lm(y ~ ., data = df_reduced),
    error = function(e) NULL
  )
  if (is.null(fit_reduced)) next
  
  #-----------------------------------
  # Full model: embeddings + PCs
  #-----------------------------------
  df_full <- as.data.frame(cbind(y = y, X_emb, X_pcs))
  fit_full <- tryCatch(
    lm(y ~ ., data = df_full),
    error = function(e) NULL
  )
  if (is.null(fit_full)) next
  
  # Extract metrics
  r2_reduced <- summary(fit_reduced)$r.squared
  r2_full <- summary(fit_full)$r.squared
  p_value <- anova(fit_reduced, fit_full)$`Pr(>F)`[2]
  loglik_ratio <- logLik(fit_full)[1] - logLik(fit_reduced)[1]
  
  #-----------------------------------
  # Extract coefficients
  #-----------------------------------
  coefs <- as.data.frame(coef(fit_full)) %>%
    tibble::rownames_to_column("Predictor")
  colnames(coefs)[2] <- "Estimate"
  coefs <- coefs %>%
    dplyr::filter(Predictor != "(Intercept)") %>%
    dplyr::mutate(
      Gene             = sub$Gene[1],
      Aligned_Position = sub$Aligned_Position[1],
      Residue          = sub$Residue[1],
      N                = nrow(sub),
      R2_full          = r2_full,
      r2_reduced       = r2_reduced,
      P_value          = p_value,
      LogLik_Ratio     = loglik_ratio
    )
  
  results_list[[i]] <- coefs
  if (i %% 50 == 0) message("Processed ", i, "/", length(groups))
}

saveRDS(results_list, "results/results_list_cve.rds")
results_list <- readRDS("results/results_list_cve.rds")
results_df <- do.call(rbind, results_list)
saveRDS(results_df, "results/residue_predictor_coefs_cve.rds")

hist(results_df$LogLik_Ratio)
hist(results_df$R2_full)
hist(results_df$r2_reduced)
line(results_df$r2_reduced, results_df$R2_full)

results_df$CVMRatio <- results_df$CVE_Reduced / results_df$CVE_Full
hist(results_df$CVMRatio)
hist(results_df$CVE_Full)
hist(results_df$CVE_Reduced)
results_df$Predictor <- as.factor(results_df$Predictor)
plot(results_df[results_df$Aligned_Position==238, ]$Estimate)
boxplot(Estimate ~ Predictor,
        data=results_df[grep(".embedding", results_df$Predictor),])


df_plot <- results_df %>%
  filter(grepl("embedding", Predictor)) %>%
  mutate(Predictor=factor(Predictor, levels=unique(Predictor)),
         Aligned_Position=as.numeric(Aligned_Position))

ggplot(df_plot, aes(x=Aligned_Position, y=Predictor, fill=CVMRatio)) +
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

df_bar <- df_plot %>%
  group_by(Gene, Aligned_Position) %>%
  summarise(
    TotalAbsEffect = sum(abs(Estimate), na.rm = TRUE),
    LogLik_Ratio = mean(LogLik_Ratio, na.rm = TRUE),
    .groups = "drop"
  )


plot(df_bar$TotalAbsEffect, df_bar$LogLik_Ratio)
text(df_bar$TotalAbsEffect, df_bar$LogLik_Ratio, df_bar$Aligned_Position)

ggplot(df_bar, aes(x=Aligned_Position, y=LogLik_Ratio)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Aligned_Position", y="CVMRatio") +
  scale_x_continuous(breaks=seq(min(df_bar$Aligned_Position), max(df_bar$Aligned_Position), by=10)) +
  coord_cartesian(ylim=range(df_bar$LogLik_Ratio))

ggplot(df_bar, aes(x=Aligned_Position, y=LogLikRatio)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Aligned_Position", y="LogLikRatio") +
  scale_x_continuous(breaks=seq(min(df_bar$Aligned_Position), max(df_bar$Aligned_Position), by=10))


data[grep("Arabidopsis thaliana", data$Organism), "ID"] #AP000423.1
at_ID <- "AP000423.1"

at_df <- df_joined[df_joined$ID==at_ID,] %>% select(Residue_Index, Aligned_Position)
str(at_df)
library(jsonlite)

feat_json <- 'data/psbA_features.json'
j <- fromJSON(feat_json) 
feats <- as.data.frame(j) 
feat_df <- feats %>%
  transmute(Residue_Index = features.location$start$value,
            Label = features.ligand$name) %>%
  left_join(at_df, by="Residue_Index") %>%
  filter(!is.na(Aligned_Position))


#a$location$end$value
#a$ligand$name

lines(df_bar$LogLikRatio)


feats$loc <- feats$features.location$start$value
feats$labl <- feats$features.ligand$name
feat_df <- feats %>%
  left_join(at_df, by=c("loc"="Residue_Index")) %>%
  filter(!is.na(Aligned_Position))

p <- ggplot(df_bar, aes(x=Aligned_Position, y=LogLik_Ratio)) +
  geom_col(fill="#2C7BB6") +
  facet_wrap(~Gene, scales="free_x") +
  theme_minimal(base_size=12) +
  labs(x="Aligned_Position", y="LogLikRatio")

p + geom_vline(data=feat_df, aes(xintercept=Aligned_Position),
               color="peachpuff", linetype="dashed", linewidth=0.7) +
  geom_text(data=feat_df,
            aes(x=Aligned_Position, y=max(df_bar$LogLik_Ratio)*0.7,
                label=labl),
            angle=90, vjust=0, hjust=0, size=5, color="red")


cor(df_bar$LogLikRatio, df_bar$TotalAbsEffect)

cor(results_df$R2, results_df$LogLikRatio)

a <- results_df %>%
  select(Aligned_Position, LogLik_Ratio) %>%
  distinct()

feat_df <- feats %>%
  transmute(Residue_Index = features.location$start$value,
            Label = features.ligand$name) %>%
  left_join(at_df, by="Residue_Index") %>%
  filter(!is.na(Aligned_Position))




ggplot(a, aes(x=Aligned_Position, y=LogLik_Ratio)) +
  geom_line(color="#2C7BB6", linewidth=0.8) +
  geom_vline(data=feat_df, aes(xintercept=Aligned_Position),
             color="orange", linetype="dashed", linewidth=0.6) +
  geom_text(data=feat_df,
            aes(x=Aligned_Position, y=0.8, label=Label),
            angle=90, vjust=0, hjust=0, size=3, color="black") +
  theme_minimal(base_size=12) +
  labs(x="Aligned Position", y=expression(R^2))

library(zoo)

# rolling average of LogLikRatio over 15-position window
a <- a %>%
  arrange(Aligned_Position) %>%
  mutate(LogLik_Ratio_smooth = rollapply(LogLik_Ratio, width = 15, FUN = mean, align = "center", fill = NA))

plot(df_bar$TotalAbsEffect, df_bar$LogLik_Ratio)

plot(a$Aligned_Position,a$LogLik_Ratio_smooth, col="white",
     main="Mean Cross Validated Error Ratio (Reduced / Full)",
     xlab="Aligned position", ylab="LogLik_Ratio_smooth")
lines(a$Aligned_Position,a$LogLik_Ratio_smooth, col="blue")
e_prior <- 0
for (i in 1:nrow(tm_df)) {
  s <- tm_df[i,]$start_aln
  e <- tm_df[i,]$end_aln
  abline(v=s, col="coral")
  abline(v=e, col="coral")
  sset <- a[a$Aligned_Position <= s &a$Aligned_Position >=e_prior,]
  m <- mean(sset$LogLik_Ratio)
  segments(e_prior, m, s, m, col="tomato")
  e_prior <-e
}
s <- e_prior
sset <- a[a$Aligned_Position <= max(a$Aligned_Position) & a$Aligned_Position >=s,]
m <- mean(sset$CVMRatio)
segments(e_prior, m, max(a$Aligned_Position), m, col="tomato")
text(305, 1000, "DE-loop", col="tomato")
for (i in 1:nrow(feat_df)) {
  abline(v = feat_df[i,]$Aligned_Position, col="green", lty=2, lwd=0.6)
  text(feat_df[i,]$Aligned_Position, 1.01, feat_df[i,]$Label, 
       srt=90, adj=c(0,0), cex=0.8, col="black")
}

summary(a$CVMRatio)



ggplot(a, aes(x = Aligned_Position, y = CVMRatio_smooth)) +
  geom_line(color = "#2C7BB6", linewidth = 0.8) +
  geom_vline(data = feat_df, aes(xintercept = Aligned_Position),
             color = "orange", linetype = "dashed", linewidth = 0.6) +
  geom_text(data = feat_df,
            aes(x = Aligned_Position, y = max(a$CVMRatio_smooth, na.rm=TRUE)*0.95, label = Label),
            angle = 90, vjust = 0, hjust = 0, size = 3, color = "black") +
  scale_x_continuous(limits = c(0, 400)) +   # <— constrain to 400 positions
  theme_minimal(base_size = 12) +
  labs(x = "Aligned Position", y = expression(R^2))

library(jsonlite)
library(dplyr)
library(ggplot2)
library(zoo)

# --- Load UniProt JSON ---
j <- fromJSON("data/psbA_TM_uniprot.json")

# Extract start, end for Transmembrane regions
tm_uniprot <- j$features %>%
  filter(type == "Transmembrane") %>%
  transmute(
    start = location$start$value,
    end   = location$end$value
  )
# --- Map to your alignment ---
at_ID <- "AP000423.1"
at_df <- df_joined %>%
  filter(ID == at_ID) %>%
  select(Residue_Index, Aligned_Position)

at_start <- dplyr::rename(at_df, start_aln = Aligned_Position)
at_end   <- dplyr::rename(at_df, end_aln = Aligned_Position)

tm_df <- tm_uniprot %>%
  left_join(at_start, by = c("start" = "Residue_Index")) %>%
  left_join(at_end,   by = c("end"   = "Residue_Index")) %>%
  filter(!is.na(start_aln), !is.na(end_aln)) %>%
  mutate(Label = "TM helix")

# --- Smooth your CVM ratio ---
a <- results_df %>%
  select(Aligned_Position, CVMRatio) %>%
  distinct() %>%
  arrange(Aligned_Position) %>%
  mutate(CVMRatio_smooth = rollapply(CVMRatio, width=15, FUN=mean, align="center", fill=NA))

# --- Plot with TM regions ---
ggplot(a, aes(Aligned_Position, CVMRatio_smooth)) +
  geom_line(color="#2C7BB6", linewidth=0.8) +
  geom_rect(
    data=tm_df,
    aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
    inherit.aes=FALSE,
    fill="orange", alpha=0.15
  ) +
  geom_text(
    data=tm_df,
    aes(x=(start+end)/2, y=max(a$CVMRatio_smooth, na.rm=TRUE)*0.95, label=Label),
    color="orange", angle=90, size=3, vjust=0
  ) +
  theme_minimal(base_size=12) +
  labs(x="Aligned Position", y=expression(R^2)) +
  scale_x_continuous(limits=c(0,400))









a2 <- a %>% mutate(logCVM = log2(CVMRatio))
ggplot(a2, aes(Aligned_Position, logCVM)) +
  geom_line(color="#2C7BB6") +
  labs(y="log2(CVMRatio)") +
  theme_minimal()

a$dist_to_bind <- sapply(a$Aligned_Position, function(x)
  min(abs(x - feat_df$Aligned_Position)))

ggplot(a, aes(x=dist_to_bind, y=R2)) +
  geom_point(alpha=0.6, color="#2C7BB6") +
  geom_smooth(method="loess", se=FALSE, color="black") +
  theme_minimal(base_size=12) +
  labs(x="Distance to Nearest Binding Site (aligned residues)", y=expression(R^2))

LSD::heatscatter(as.vector(a$dist_to_bind),as.vector(a$R2))
cor(as.vector(a$dist_to_bind),as.vector(a$R2), method="pearson")
