library(arrow)
library(data.table)
library(Biostrings)

# ---------------------------------------------------------------------
# LOAD COMMON DATA
# ---------------------------------------------------------------------

data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")

# Load embedding QC data
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))

# Create clean ID sets per gene (exclude ManualOutlier == TRUE)
setDT(embeds_with_mds)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

# Prepare PC table
pcs_IDS <- aln_index$index
scores <- as.data.table(ev_pcs$x)
setnames(scores, paste0("PC", seq_len(ncol(scores))))
scores[, ID := pcs_IDS]

n_pcs <- 1000L
pc_names <- paste0("PC", seq_len(n_pcs))
scores <- scores[, c("ID", pc_names), with = FALSE]

# Pre-join phenotype + PCs once
pheno_pcs <- data[, .(ID, pheno = get(pheno_col))][scores, on = "ID", nomatch = 0]
stopifnot(nrow(pheno_pcs) > 0)

# ---------------------------------------------------------------------
# PREPARE FILE LISTS
# ---------------------------------------------------------------------

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes_to_process <- get_gene(aln_files)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes: ", paste(genes_to_process, collapse = ", "))

# ---------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------

set.seed(123)

for (gene in genes_to_process) {
  message("\n=== Processing gene: ", gene, " ===")
  
  out_path <- paste0("results/residue_models_clean/", gene, "_effects.rds")
  
  if (file.exists(out_path)) {
    message("Skipping ", gene, " (output already exists)")
    next
  }
  
  # Get clean IDs for this gene
  clean_ids_gene <- clean_ids_by_gene[Gene == gene, ID]
  if (length(clean_ids_gene) == 0) {
    message("Skipping ", gene, ": no clean IDs after filtering ManualOutlier")
    next
  }
  
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # --- Read alignment ---
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln) || length(aln) == 0) {
    message("Skipping ", gene, ": could not read alignment")
    next
  }
  
  # --- Strip IDs ---
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  seq_ids <- names(aln)
  
  # Get common IDs: intersection of (alignment IDs, phenotype IDs, clean IDs for this gene)
  common_ids <- intersect(intersect(seq_ids, pheno_pcs$ID), clean_ids_gene)
  if (length(common_ids) < 8500L) {
    message("Skipping ", gene, ": insufficient overlap after filtering outliers (n=", 
            length(common_ids), ")")
    next
  }
  
  # Subset to common IDs
  aln_mat <- aln_mat[match(common_ids, seq_ids), , drop = FALSE]
  pheno_pcs_sub <- pheno_pcs[ID %in% common_ids]
  setkey(pheno_pcs_sub, ID)
  pheno_pcs_sub <- pheno_pcs_sub[match(common_ids, ID)]
  
  y <- pheno_pcs_sub$pheno
  X_pcs <- scale(as.matrix(pheno_pcs_sub[, ..pc_names]))
  
  # -------------------------------------------------------------------
  # LOOP OVER POSITIONS
  # -------------------------------------------------------------------
  
  positions <- seq_len(ncol(aln_mat))
  results_list <- vector("list", length(positions))
  
  for (pos in positions) {
    
    residues <- aln_mat[, pos]
    residue_table <- table(residues)
    
    # Skip gaps-only or monomorphic positions
    if (length(residue_table) < 2L || all(names(residue_table) == "-")) next
    
    # Skip positions with only gaps
    non_gap <- residues != "-"
    if (sum(non_gap) < 8500L) next
    
    # Check special case
    if (pos == 19 & gene == "psbJ") next
    
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    var0 <- which(apply(X_aa, 2, var) == 0)
    if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # -----------------------------------------------------------------
    # FIT MODELS
    # -----------------------------------------------------------------
    
    fit_reduced <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
    if (is.null(fit_reduced)) next
    
    fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_reduced <- summary(fit_reduced)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    p_res <- tryCatch(anova(fit_reduced, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    # Extract only residue effect sizes and SEs
    coef_summary <- summary(fit_full)$coefficients
    aa_coefs <- grep("^X_aa", rownames(coef_summary), value = TRUE)
    
    if (length(aa_coefs) > 0) {
      effects_dt <- data.table(
        Residue = sub("res_factor", "", aa_coefs),
        Effect = coef_summary[aa_coefs, "Estimate"],
        SE = coef_summary[aa_coefs, "Std. Error"],
        T_value = coef_summary[aa_coefs, "t value"],
        P_value = coef_summary[aa_coefs, "Pr(>|t|)"]
      )
    } else {
      effects_dt <- NULL
    }
    
    results_list[[pos]] <- list(
      Gene = gene,
      Aligned_Position = pos,
      N = sum(non_gap),
      R2_reduced = r2_reduced,
      R2_full = r2_full,
      P_res = p_res,
      residue_counts = as.list(residue_table),
      effects = effects_dt,
      IDs = common_ids
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) > 0) {
    dir.create("results/residue_models", showWarnings = FALSE, recursive = TRUE)
    saveRDS(results_list, out_path)
    message("\nCompleted ", gene, ": ", length(results_list), " positions with effects")
  } else {
    message("No positions retained for gene ", gene)
  }
}


# Load all model files
model_files <- list.files("results/residue_models_clean/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

f <- model_files[41]
models <- readRDS(f)
length(models)
i <- 333
models[i]

#for each site, I want to get the significance, effect size, and the amino acid proportions

cols <- c("P_res", "Effect", "Gene",
          "M", "T", "A", "I", "L", "E", "R", "S", "C", "K",
          "D", "N", "Y", "F", "W", "P", "V", "H", "Q", "G")
res <- data.frame(
  matrix(ncol=length(cols),
         nrow=length(models))
)
names(res) <- cols
aa_cols <- c(
          "M", "T", "A", "I", "L", "E", "R", "S", "C", "K",
          "D", "N", "Y", "F", "W", "P", "V", "H", "Q", "G")

f <- model_files[1]
models <- readRDS(f)
for (i in 1:length(models)) {
  mod <- models[i]
  mod <- mod[[1]]
  print(mod$effects, mod$residue_counts)
  #effect value of the significant terms, p < 0.001
  effs <- mod$effects
  sigs <- effs[effs$P_value==min(effs$P_value)]
  counts <- mod$residue_counts
  counts_vec <- unlist(counts)
  res[i,names(counts_vec)] <- counts_vec
  res[i,]$P_res <- mod$P_res
}
res[is.na(res)] <-0
hits <- res[res$P_res < 0.001,]

summary(hits)

all_res <- list()

for(f in model_files){
  models <- readRDS(f)
  res <- data.frame(matrix(ncol=length(cols), nrow=length(models)))
  names(res) <- cols
  print(f)
  
  for(i in seq_along(models)){
    mod <- models[[i]]
    pos <- mod$Aligned_Position
    effs <- mod$effects
    aa_cols <- cols[4:length(cols)]
    
    counts_vec <- unlist(mod$residue_counts)
    counts_vec <- counts_vec[names(counts_vec) %in% aa_cols]
    #hist(counts_vec, breaks = 10)
    
    res[i, aa_cols] <- 0
    res[i, names(counts_vec)] <- counts_vec
    res$P_res[i] <- mod$P_res
    effs <- mod$effects
    res$Effect[i] <- mean(effs$Effect)
    res$Gene[i] <- mod$Gene
  }
  print(sum(is.na(res[is.na(res)])))
  #res[is.na(res)] <- 0
  all_res[[f]] <- res
  #print(summary(res))
}

final_res <- do.call(rbind, all_res)
hits <- final_res[final_res$P_res < quantile(final_res$P_res,0.05),]
non_hits <- final_res[final_res$P_res > quantile(final_res$P_res,0.05),]
df <- final_res
df$logp <- -log10(df$P_res)
df$group <- ifelse(df$P_res < quantile(df$P_res, 0.05), "hit", "non_hit")

df <- df[order(df$P_res), ]

qq <- data.frame(
  expected = -log10(ppoints(nrow(df))),
  observed = df$logp,
  group = df$group
)

ggplot(qq, aes(expected, observed, color = group)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values=c("hit"="red","non_hit"="black")) +
  theme_bw()

df <- final_res
df$logp <- -log10(df$P_res)
df$group <- ifelse(df$P_res < quantile(df$P_res,0.05),"hit","non_hit")
df <- df[order(df$P_res),]

expected <- -log10(ppoints(nrow(df)))
observed <- df$logp

par(mfrow=c(1,2))
plot(expected, observed, pch=19, col=ifelse(df$group=="hit","red","black"),
     main="GWAS hit selection",
     xlab="Expected -log10(P)", ylab="Observed -log10(P)")
abline(0,1)
hist(final_res$P_res, main = "", xlab="P")
abline(v = quantile(final_res$P_res, 0.05), col = "red")
text(0.2, 1000, round(quantile(final_res$P_res, 0.05),5), col="red")
#and it sounds like there are some geometric reasonons for that 

plot(hits$Effect, -log10(hits$P_res))
summary(hits)


df_long <- hits |> 
  select(all_of(aa_cols)) |> 
  pivot_longer(everything(), names_to="AA", values_to="count") |>
  filter(count > 10)

ggplot(df_long, aes(AA, count)) +
  geom_boxplot() +
  #scale_y_log10() +
  geom_jitter() +
  ggtitle("Amino acid counts in top decile GWAS hits") +
  theme_bw()

LSD::heatboxplot(count ~ AA, df_long)

aa_list <- split(df_long$count, df_long$AA)
totals <- sapply(aa_list, sum)


layout(matrix(c(1,2), nrow=1), widths=c(3,1))
heatboxplot(
  x = aa_list,
  colpals = rep("heat", length(aa_list)),
  labels = names(aa_list),
  horizontal = T,
  ylim = c(min(df_long$count), max(df_long$count)),
  main = "AA counts (hits only)",
  axes = TRUE
)
barplot(totals, horiz=TRUE, las=1)

for (aa in aa_cols) {
  hist(df_long$count[df_long$AA==aa], main=aa)
}

#protein product functional enrichment 
hits$complex <- substr(hits$Gene, 1, 3)
non_hits$complex <- substr(non_hits$Gene, 1, 3)
df_long <- hits |>
  select(Gene, complex, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to="AA", values_to="count") |>
  filter(count > 10)

non_hits_long <-  non_hits |>
  select(Gene, complex, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to="AA", values_to="count") |>
  filter(count > 10)

by_complex <- split(df_long, df_long$complex)
#FROM WHAT AMINO ACIDS ARE WE GOING (major)
majs <- df_long[df_long$count > 4000,]
by_complex <- split(majs, majs$complex)


for(cmp in names(by_complex)){
  x <- split(by_complex[[cmp]]$count, by_complex[[cmp]]$AA)
  aa_counts <- sapply(x, length)
  barplot(aa_counts, las=2, main=paste(cmp, "majs"))
}

mins <- df_long[df_long$count < 4000,]
by_complex <- split(mins, mins$complex)


for(cmp in names(by_complex)){
  x <- split(by_complex[[cmp]]$count, by_complex[[cmp]]$AA)
  aa_counts <- sapply(x, length)
  barplot(aa_counts, las=2, main=paste(cmp, "mins"))
}

# ---- enrichment analysuis --- 
bg_tot <- colSums(non_hits[, aa_cols])
barplot(bg_tot)
cor(bg_tot,all_bg_tot)

aa_tot <- colSums(hits[, aa_cols])
barplot(aa_tot)
plot(aa_tot, bg_tot,
     xlab="AA counts gwas hits",
     ylab="AA counts all sites",col="white",
     main="AA counts in GWAS hits vs background")
text(aa_tot, bg_tot, aa_cols)

obs_prop <- aa_tot / sum(aa_tot)
bg_prop  <- bg_tot / sum(bg_tot)


enrich <- obs_prop / bg_prop
barplot(enrich)
hist(enrich)
sort(enrich, decreasing=TRUE)
#In general, there is min #I SHOULD BE WRITING THIS DOWN 
#froms and tos, baby, froms and tos
aa_tot <- colSums(hits[, aa_cols])
barplot(aa_tot)

hist(rowSums(hits[,aa_cols]!=0), main="Number of unique amino acids present at each site")
table(rowSums(hits[,aa_cols]!=0))

obs_prop <- aa_tot / sum(aa_tot)
bg_prop  <- bg_tot / sum(bg_tot)
enrich <- obs_prop / bg_prop
barplot(enrich)
hist(enrich)
sort(enrich, decreasing=TRUE)

clean_long <- df_long[df_long$AA!="-",]

FROM <- clean_long[clean_long$count > 4000,]
from_bg <- non_hits_long[non_hits_long$count > 4000,]

  
TO   <- clean_long[clean_long$count < 4000,]
to_bg <- non_hits_long[non_hits_long$count < 4000,]


from_counts <- table(FROM$AA)
from_bg_counts <- table(from_bg$AA)
from_bg_prop  <- from_bg_counts / sum(from_bg_counts)

to_counts   <- table(TO$AA)
to_bg_counts <- table(to_bg$AA)
to_bg_prop  <- to_bg_counts / sum(to_bg_counts)


bg_tot <- bg_tot[names(from_counts)]
chisq.test(from_counts, p = bg_tot / sum(bg_tot))
chisq.test(to_counts, p = bg_tot / sum(bg_tot))

aa <- "C"

k <- sum(TO$AA == aa)
n <- nrow(TO)
p0 <- to_bg_counts[aa] / sum(to_bg_counts)

binom.test(k, n, p = p0)

aa_list <- names(to_counts)

pvals <- sapply(aa_list, function(a){
  k <- to_counts[a]
  n <- sum(to_counts)
  p0 <- to_bg_counts[a] / sum(to_bg_counts)
  binom.test(k, n, p=p0)$p.value
})

to_prop <- to_counts / sum(to_counts)
to_enrich <- to_prop / to_bg_prop[names(to_prop)]

stars <- sapply(pvals, function(p){
  if(p < 0.001) "***"
  else if(p < 0.01) "**"
  else if(p < 0.05) "*"
  else ""
})

bp <- barplot(to_enrich, las=2,
              main="Minor allele enrichment in top 10%",ylim=c(0,1.5))
text(bp, to_enrich, labels=stars, pos=3, cex=0.9)

aa_list <- names(from_counts)

pvals <- sapply(aa_list, function(a){
  k <- from_counts[a]
  n <- sum(from_counts)
  p0 <- from_bg_counts[a] / sum(from_bg_counts)
  binom.test(k, n, p=p0)$p.value
})

from_prop <- from_counts / sum(from_counts)
from_enrich <- from_prop / from_bg_prop[names(from_bg_prop)]

stars <- sapply(pvals, function(p){
  if(p < 0.001) "***"
  else if(p < 0.01) "**"
  else if(p < 0.05) "*"
  else ""
})

bp <- barplot(from_enrich, las=2, ylim=c(0,1.7),
              main="Major allele enrichment in top 10%")
text(bp, from_enrich, labels=stars, pos=3, cex=0.9)


# ---- AA anova ----
# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Bind results and classify by GWAS significance
final_res <- do.call(rbind, all_res)
stopifnot(nrow(final_res) > 0)
stopifnot("P_res" %in% names(final_res))

threshold <- quantile(final_res$P_res, 0.1)
final_res$gwas_class <- ifelse(final_res$P_res < threshold, "hit", "non_hit")

# Create long format with allele class
df_processed <- final_res |>
  mutate(site = row_number()) |>
  select(site, gwas_class, Effect, all_of(aa_cols)) |>
  pivot_longer(cols = all_of(aa_cols), names_to = "AA", values_to = "count") |>
  filter(AA != "-", count > 10) |>
  mutate(allele_class = ifelse(count > 4000, "major", "minor")) |>
  group_by(site, gwas_class, Effect, allele_class, AA) |>
  summarise(count = sum(count), .groups = "drop") |>
  group_by(site, gwas_class, Effect, allele_class) |>
  mutate(proportion = count / sum(count)) |>
  ungroup() |>
  filter(!(allele_class == "major" & proportion < 0.9))

summary(df_processed)
summary(df_processed[df_processed$allele_class == "major",])
summary(df_processed[df_processed$allele_class == "minor",])
stopifnot(all(df_processed$proportion >= 0 & df_processed$proportion <= 1))
stopifnot(all(c("major", "minor") %in% df_processed$allele_class))

# Pivot to wide format for analysis
df_wide <- df_processed |>
  select(site, gwas_class, allele_class, AA, proportion) |>
  pivot_wider(
    names_from = AA,
    values_from = proportion,
    values_fill = 0
  )

# ANOVA for each amino acid by allele class and GWAS status
aa_anova_results <- lapply(aa_cols[aa_cols != "-"], function(amino_acid) {
  if (amino_acid %in% names(df_wide)) {
    formula_str <- paste(amino_acid, "~ gwas_class * allele_class")
    aov_result <- aov(as.formula(formula_str), data = df_wide)
    data.frame(
      AA = amino_acid,
      summary(aov_result)[[1]],
      term = rownames(summary(aov_result)[[1]])
    )
  }
})

anova_summary <- do.call(rbind, aa_anova_results)

hist(df_processed$proportion[df_processed$allele_class=="major"])

ggplot(df_processed, aes(x = AA, y = proportion, fill = gwas_class)) +
  geom_boxplot() +
  facet_grid(~allele_class, labeller = label_both) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "AA proportions by allele class and GWAS status",
    y = "Proportion"
  )

df_comparison <- df_processed |>
  filter(allele_class == "minor") |>
  group_by(AA, gwas_class) |>
  summarise(mean_prop = mean(proportion), .groups = "drop") |>
  pivot_wider(names_from = gwas_class, values_from = mean_prop)

stopifnot(all(c("hit", "non_hit") %in% names(df_comparison)))

ggplot(df_comparison, aes(x = hit, y = non_hit, label = AA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_text(hjust = -0.2, vjust = -0.2) +
  theme_bw() +
  labs(
    title = "Mean AA proportions: hits vs non-hits (minor alleles)",
    x = "Hit mean proportion",
    y = "Non-Hit mean proportion"
  ) +
  coord_fixed()


for (aa in unique(df_processed$AA)) {
  print(paste("-------------", aa, "-------------"))
  sset <- df_processed[df_processed$AA==aa & df_processed$allele_class == "minor",]
  print(summary(aov(proportion ~ gwas_class, sset)))
}
# Calculate ANOVA p-values for minor sites only
minor_anova <- df_processed |>
  filter(allele_class == "minor") |>
  group_by(AA) |>
  summarise(
    p_value = summary(aov(proportion ~ gwas_class))[[1]]["gwas_class", "Pr(>F)"],
    .groups = "drop"
  ) |>
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

stopifnot(all(!is.na(minor_anova$p_value)))

# Plot with annotations
df_processed |>
  filter(allele_class == "minor") |>
  ggplot(aes(x = proportion, y = AA, fill = gwas_class)) +
  geom_density_ridges(alpha = 0.6) +
  geom_text(data = minor_anova, aes(x = 0.95, y = AA, label = p_label), 
            inherit.aes = FALSE, hjust = 1, size = 5) +
  theme_ridges() +
  labs(
    title = "AA proportion distributions at minor allele sites",
    x = "Proportion",
    y = "Amino Acid"
  )
# ---- sites wide AA prop ----

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Bind results and classify by GWAS significance
final_res <- do.call(rbind, all_res)
stopifnot(nrow(final_res) > 0)
stopifnot("P_res" %in% names(final_res))

threshold <- quantile(final_res$P_res, 0.05)
final_res$gwas_class <- ifelse(final_res$P_res < threshold, "hit", "non_hit")

# Create long format with allele class
df_long <- final_res |>
  mutate(site = row_number()) |>
  select(site, gwas_class, Effect, all_of(aa_cols)) |>
  pivot_longer(
    cols = all_of(aa_cols),
    names_to = "AA",
    values_to = "count"
  ) |>
  filter(AA != "-", count > 10) |>
  mutate(
    allele_class = ifelse(count > 4000, "major", "minor"),
    effect_class = case_when(
      gwas_class == "non_hit" ~ "no_effect",
      Effect > 0 ~ "hot",
      Effect < 0 ~ "cold",
      TRUE ~ "no_effect"
    )
  )

table(df_long$effect_class)
# Filter to keep only high proportion major alleles
major_high <- df_long |>
  filter(allele_class == "major") |>
  group_by(site, gwas_class) |>
  mutate(site_total = sum(count),
         proportion = count / site_total) |>
  filter(proportion > 0.9) |>
  ungroup()

minor_all <- df_long |>
  filter(allele_class == "minor")

df_filtered <- bind_rows(major_high, minor_all)

# Calculate overall AA proportions by gwas_class and allele_class
aa_proportions <- df_filtered |>
  group_by(effect_class, allele_class, AA) |>
  summarise(total_count = sum(count), .groups = "drop") |>
  group_by(effect_class, allele_class) |>
  mutate(
    class_total = sum(total_count),
    proportion = total_count / class_total
  ) |>
  ungroup()

stopifnot(all(aa_proportions$proportion >= 0 & aa_proportions$proportion <= 1))

# Verify proportions sum to 1 within each group
check_sums <- aa_proportions |>
  group_by(effect_class, allele_class) |>
  summarise(prop_sum = sum(proportion), .groups = "drop")
stopifnot(all(abs(check_sums$prop_sum - 1) < 0.001))

# Boxplot by allele class and gwas status
ggplot(aa_proportions, aes(x = AA, y = proportion, fill = effect_class)) +
  geom_col(position = "dodge") +
  facet_grid(~allele_class, labeller = label_both) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Overall AA proportions by allele class and GWAS status",
    y = "Proportion of total AA count"
  )

ggplot(aa_proportions, aes(x = AA, y = proportion, fill = effect_class)) +
  geom_col(position = "dodge") +
  facet_grid(~allele_class, labeller = label_both) +
  scale_fill_manual(values = c(
    "cold" = "#3B82F6",      # blue
    "hot" = "#EF4444",       # red
    "no_effect" = "#9CA3AF"  # gray
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Overall AA proportions by allele class and effect",
    y = "Proportion of total AA count",
    fill = "Effect class"
  )


par(mfrow=c(1,3))
plot(
  aa_proportions$proportion[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="no_effect"& aa_proportions$allele_class=="minor"],
  main="AA proportions in cold adaptive sites vs background", 
  xlab="aa proportion across cold adaptive sites",
  ylab="aa proportion across background sites",col="white"
)
text(
  aa_proportions$proportion[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="no_effect"& aa_proportions$allele_class=="minor"],
  aa_proportions$AA[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"]
)
abline(0,1)

plot(
  aa_proportions$proportion[aa_proportions$effect_class=="hot" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="no_effect"& aa_proportions$allele_class=="minor"],
  main="AA proportions in hot adaptive sites vs background", 
  xlab="aa proportion across hot adaptive sites",
  ylab="aa proportion across background sites",col="white"
)
text(
  aa_proportions$proportion[aa_proportions$effect_class=="hot" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="no_effect"& aa_proportions$allele_class=="minor"],
  aa_proportions$AA[aa_proportions$effect_class=="hot" & aa_proportions$allele_class=="minor"]
)
abline(0,1)

plot(
  aa_proportions$proportion[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="hot"& aa_proportions$allele_class=="minor"],
  main="AA proportions in cold vs hot adaptive sites", 
  xlab="aa proportion across cold adaptive sites",
  ylab="aa proportion across hot adaptive sites",col="white"
)
text(
  aa_proportions$proportion[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$effect_class=="hot"& aa_proportions$allele_class=="minor"],
  aa_proportions$AA[aa_proportions$effect_class=="cold" & aa_proportions$allele_class=="minor"]
)
abline(0,1)



# Scatterplot comparing hit vs non-hit for minor alleles
aa_comparison <- aa_proportions |>
  filter(allele_class == "minor") |>
  select(AA, gwas_class, proportion) |>
  pivot_wider(names_from = gwas_class, values_from = proportion)

stopifnot(all(c("hit", "non_hit") %in% names(aa_comparison)))

ggplot(aa_comparison, aes(x = hit, y = non_hit, label = AA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_text(hjust = -0.2, vjust = -0.2) +
  theme_bw() +
  labs(
    title = "Overall AA proportions: hits vs non-hits (minor alleles)",
    x = "Hit proportion",
    y = "Non-hit proportion"
  ) +
  coord_fixed()

# Chi-square test for each allele class
for (ac in c("major", "minor")) {
  cat("\n========", ac, "alleles ========\n")
  
  contingency <- aa_proportions |>
    filter(allele_class == ac) |>
    select(AA, gwas_class, total_count) |>
    pivot_wider(names_from = gwas_class, values_from = total_count, values_fill = 0) |>
    column_to_rownames("AA")
  
  chi_result <- chisq.test(as.matrix(contingency))
  print(chi_result)
}
enrichment <- aa_comparison |>
  mutate(enrichment = hit / non_hit) |>
  arrange(desc(enrichment))

barplot(enrichment$enrichment, names.arg = enrichment$AA, las = 2,
        main = "AA enrichment in adaptive sites (minor alleles)",
        ylab = "Fold enrichment")
abline(h = 1, lty = 2)

# ---- adap vs background ---- 
# Bind results and classify by GWAS significance
final_res <- do.call(rbind, all_res)
stopifnot(nrow(final_res) > 0)
stopifnot("P_res" %in% names(final_res))

threshold <- quantile(final_res$P_res, 0.05)
final_res$gwas_class <- ifelse(final_res$P_res < threshold, "adaptive", "background")

# Create long format with allele class
df_long <- final_res |>
  mutate(site = row_number()) |>
  select(site, gwas_class, all_of(aa_cols)) |>
  pivot_longer(
    cols = all_of(aa_cols),
    names_to = "AA",
    values_to = "count"
  ) |>
  filter(AA != "-", count > 10) |>
  mutate(allele_class = ifelse(count > 4000, "major", "minor"))

# Filter to keep only high proportion major alleles
major_high <- df_long |>
  filter(allele_class == "major") |>
  group_by(site, gwas_class) |>
  mutate(site_total = sum(count),
         proportion = count / site_total) |>
  filter(proportion > 0.9) |>
  ungroup()

minor_all <- df_long |>
  filter(allele_class == "minor")

df_filtered <- bind_rows(major_high, minor_all)

# Calculate overall AA proportions by gwas_class and allele_class
aa_proportions <- df_filtered |>
  group_by(gwas_class, allele_class, AA) |>
  summarise(total_count = sum(count), .groups = "drop") |>
  group_by(gwas_class, allele_class) |>
  mutate(
    class_total = sum(total_count),
    proportion = total_count / class_total
  ) |>
  ungroup()

stopifnot(all(aa_proportions$proportion >= 0 & aa_proportions$proportion <= 1))

# Verify proportions sum to 1 within each group
check_sums <- aa_proportions |>
  group_by(gwas_class, allele_class) |>
  summarise(prop_sum = sum(proportion), .groups = "drop")
stopifnot(all(abs(check_sums$prop_sum - 1) < 0.001))

# Barplot by allele class and gwas status
ggplot(aa_proportions, aes(x = AA, y = proportion, fill = gwas_class)) +
  geom_col(position = "dodge") +
  facet_grid(~allele_class, labeller = label_both) +
  scale_fill_manual(values = c(
    "adaptive" = "#EF4444",      # red
    "background" = "#9CA3AF"     # gray
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Overall AA proportions by allele class and GWAS status",
    y = "Proportion of total AA count",
    fill = "GWAS class"
  )

# Scatterplot comparing adaptive vs background
par(mfrow=c(1,2))
plot(
  aa_proportions$proportion[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="major"],
  aa_proportions$proportion[aa_proportions$gwas_class=="background" & aa_proportions$allele_class=="major"],
  main="Major allele AA proportions in adaptive sites vs background", 
  xlab="AA proportion across adaptive sites",
  ylab="AA proportion across background sites",
  col="white"
)
text(
  aa_proportions$proportion[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="major"],
  aa_proportions$proportion[aa_proportions$gwas_class=="background" & aa_proportions$allele_class=="major"],
  aa_proportions$AA[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="major"]
)
abline(0, 1)
plot(
  aa_proportions$proportion[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$gwas_class=="background" & aa_proportions$allele_class=="minor"],
  main="Minor allele AA proportions in adaptive sites vs background", 
  xlab="AA proportion across adaptive sites",
  ylab="AA proportion across background sites",
  col="white"
)
text(
  aa_proportions$proportion[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="minor"],
  aa_proportions$proportion[aa_proportions$gwas_class=="background" & aa_proportions$allele_class=="minor"],
  aa_proportions$AA[aa_proportions$gwas_class=="adaptive" & aa_proportions$allele_class=="minor"]
)
abline(0, 1)



# ---- enrichment ----
# Identify hits and non-hits
hits <- final_res[final_res$P_res < quantile(final_res$P_res, 0.1), ]
non_hits <- final_res[final_res$P_res > quantile(final_res$P_res, 0.1), ]

stopifnot(nrow(hits) > 0, nrow(non_hits) > 0)
stopifnot(nrow(hits) + nrow(non_hits) == nrow(final_res))

# Prepare data: filter for count > 10, exclude "-"
hits_long <- hits |>
  select(Gene, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to = "AA", values_to = "count") |>
  filter(count > 10, AA != "-")

non_hits_long <- non_hits |>
  select(Gene, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to = "AA", values_to = "count") |>
  filter(count > 10, AA != "-")

stopifnot(nrow(hits_long) > 0, nrow(non_hits_long) > 0)

# Define major (FROM) and minor (TO) alleles
FROM <- hits_long[hits_long$count > 4000, ]
FROM_bg <- non_hits_long[non_hits_long$count > 4000, ]

TO <- hits_long[hits_long$count < 4000, ]
TO_bg <- non_hits_long[non_hits_long$count < 4000, ]

stopifnot(nrow(FROM) > 0, nrow(TO) > 0)
stopifnot(nrow(FROM_bg) > 0, nrow(TO_bg) > 0)

# Calculate enrichment for minor alleles (TO)
to_counts <- table(TO$AA)
to_bg_counts <- table(TO_bg$AA)

# Align amino acids
shared_aa <- intersect(names(to_counts), names(to_bg_counts))
stopifnot(length(shared_aa) > 0)

to_counts <- to_counts[shared_aa]
to_bg_counts <- to_bg_counts[shared_aa]

# Test and calculate enrichment
to_prop <- to_counts / sum(to_counts)
to_bg_prop <- to_bg_counts / sum(to_bg_counts)
to_enrich <- to_prop / to_bg_prop

pvals_to <- sapply(shared_aa, function(a) {
  binom.test(to_counts[a], sum(to_counts), 
             p = to_bg_counts[a] / sum(to_bg_counts))$p.value
})

stars_to <- ifelse(pvals_to < 0.001, "***",
                   ifelse(pvals_to < 0.01, "**",
                          ifelse(pvals_to < 0.05, "*", "")))

bp <- barplot(to_enrich, las = 2, ylim = c(0, max(to_enrich) * 1.2),
              main = "Minor allele enrichment in top 10%",
              ylab = "Enrichment")
abline(h = 1, lty = 2, col = "gray")
text(bp, to_enrich, labels = stars_to, pos = 3, cex = 0.9)

# Calculate enrichment for major alleles (FROM)
from_counts <- table(FROM$AA)
from_bg_counts <- table(FROM_bg$AA)

shared_aa_from <- intersect(names(from_counts), names(from_bg_counts))
stopifnot(length(shared_aa_from) > 0)

from_counts <- from_counts[shared_aa_from]
from_bg_counts <- from_bg_counts[shared_aa_from]

from_prop <- from_counts / sum(from_counts)
from_bg_prop <- from_bg_counts / sum(from_bg_counts)
from_enrich <- from_prop / from_bg_prop

pvals_from <- sapply(shared_aa_from, function(a) {
  binom.test(from_counts[a], sum(from_counts),
             p = from_bg_counts[a] / sum(from_bg_counts))$p.value
})

stars_from <- ifelse(pvals_from < 0.001, "***",
                     ifelse(pvals_from < 0.01, "**",
                            ifelse(pvals_from < 0.05, "*", "")))

bp <- barplot(from_enrich, las = 2, ylim = c(0, max(from_enrich) * 1.2),
              main = "Major allele enrichment in top 10%",
              ylab = "Enrichment")
abline(h = 1, lty = 2, col = "gray")
text(bp, from_enrich, labels = stars_from, pos = 3, cex = 0.9)


# ---- protein complex enrichment ---- 

hits$complex <- substr(hits$Gene, 1, 3)
non_hits$complex <- substr(non_hits$Gene, 1, 3)

# Prepare long format data
hits_long <- hits |>
  select(Gene, complex, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to = "AA", values_to = "count") |>
  filter(count > 10, AA != "-")

non_hits_long <- non_hits |>
  select(Gene, complex, all_of(aa_cols)) |>
  pivot_longer(all_of(aa_cols), names_to = "AA", values_to = "count") |>
  filter(count > 10, AA != "-")

stopifnot(nrow(hits_long) > 0, nrow(non_hits_long) > 0)

# Define major (FROM) and minor (TO) alleles
FROM <- hits_long[hits_long$count > 4000, ]
FROM_bg <- non_hits_long[non_hits_long$count > 4000, ]

TO <- hits_long[hits_long$count < 4000, ]
TO_bg <- non_hits_long[non_hits_long$count < 4000, ]

stopifnot(nrow(FROM) > 0, nrow(TO) > 0)
stopifnot(nrow(FROM_bg) > 0, nrow(TO_bg) > 0)


# Calculate global background proportions
to_bg_counts <- table(TO_bg$AA)
from_bg_counts <- table(FROM_bg$AA)

to_bg_prop <- to_bg_counts / sum(to_bg_counts)
from_bg_prop <- from_bg_counts / sum(from_bg_counts)

# Function to calculate enrichment for one complex
calc_enrichment <- function(complex_data, bg_prop, label) {
  counts <- table(complex_data$AA)
  shared_aa <- intersect(names(counts), names(bg_prop))
  
  if (length(shared_aa) == 0) return(NULL)
  
  counts <- counts[shared_aa]
  bg_prop_sub <- bg_prop[shared_aa]
  
  prop <- counts / sum(counts)
  enrich <- prop / bg_prop_sub
  
  pvals <- sapply(shared_aa, function(a) {
    binom.test(counts[a], sum(counts), p = bg_prop_sub[a])$p.value
  })
  
  stars <- ifelse(pvals < 0.001, "***",
                  ifelse(pvals < 0.01, "**",
                         ifelse(pvals < 0.05, "*", "")))
  
  list(enrich = enrich, pvals = pvals, stars = stars, n = sum(counts))
}

# Analyze by complex with global background
complexes <- c("psa", "psb", "atp", "rpo", "rps", "ndh")
complexes <- complexes[order(complexes)]

par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))

for (cmp in complexes) {
  # Minor alleles (TO)
  to_cmp <- TO[TO$complex == cmp, ]
  
  if (nrow(to_cmp) > 5) {
    res <- calc_enrichment(to_cmp, to_bg_prop, "TO")
    if (!is.null(res)) {
      bp <- barplot(res$enrich, las = 2, ylim = c(0, max(res$enrich) * 1.2),
                    main = paste0(cmp, " - Minor allele (n=", res$n, ")"),
                    ylab = "Enrichment vs global")
      abline(h = 1, lty = 2, col = "gray")
      text(bp, res$enrich, labels = res$stars, pos = 3, cex = 0.8)
    }
  }
  
  # Major alleles (FROM)
  from_cmp <- FROM[FROM$complex == cmp, ]
  
  if (nrow(from_cmp) > 5) {
    res <- calc_enrichment(from_cmp, from_bg_prop, "FROM")
    if (!is.null(res)) {
      bp <- barplot(res$enrich, las = 2, ylim = c(0, max(res$enrich) * 1.2),
                    main = paste0(cmp, " - Major allele (n=", res$n, ")"),
                    ylab = "Enrichment vs global")
      abline(h = 1, lty = 2, col = "gray")
      text(bp, res$enrich, labels = res$stars, pos = 3, cex = 0.8)
    }
  }
}

par(mfrow = c(1, 1))

# ============================================================================
# STRATIFIED ENRICHMENT (complex-specific background)
# ============================================================================

par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))

for (cmp in complexes) {
  # Minor alleles with complex-specific background
  to_cmp <- TO[TO$complex == cmp, ]
  to_bg_cmp <- TO_bg[TO_bg$complex == cmp, ]
  
  if (nrow(to_cmp) > 5 && nrow(to_bg_cmp) > 5) {
    bg_counts <- table(to_bg_cmp$AA)
    bg_prop <- bg_counts / sum(bg_counts)
    
    res <- calc_enrichment(to_cmp, bg_prop, "TO")
    if (!is.null(res)) {
      bp <- barplot(res$enrich, las = 2, ylim = c(0, max(res$enrich) * 1.2),
                    main = paste0(cmp, " - Minor (n=", res$n, ")"),
                    ylab = "Enrichment vs complex bg")
      abline(h = 1, lty = 2, col = "gray")
      text(bp, res$enrich, labels = res$stars, pos = 3, cex = 0.8)
    }
  }
  
  # Major alleles with complex-specific background
  from_cmp <- FROM[FROM$complex == cmp, ]
  from_bg_cmp <- FROM_bg[FROM_bg$complex == cmp, ]
  
  if (nrow(from_cmp) > 5 && nrow(from_bg_cmp) > 5) {
    bg_counts <- table(from_bg_cmp$AA)
    bg_prop <- bg_counts / sum(bg_counts)
    
    res <- calc_enrichment(from_cmp, bg_prop, "FROM")
    if (!is.null(res)) {
      bp <- barplot(res$enrich, las = 2, ylim = c(0, max(res$enrich) * 1.2),
                    main = paste0(cmp, " - Major (n=", res$n, ")"),
                    ylab = "Enrichment vs complex bg")
      abline(h = 1, lty = 2, col = "gray")
      text(bp, res$enrich, labels = res$stars, pos = 3, cex = 0.8)
    }
  }
}

par(mfrow = c(1, 1))
# ---- transition matrix ----

all_trans <- list()

for(f in model_files){
  models <- readRDS(f)
  n <- length(models)
  
  trans_df <- data.frame(
    site = seq_len(n),
    Gene = character(n),
    P_res = numeric(n),
    major = character(n),
    minor_list = character(n),
    minor_counts = character(n),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_len(n)){
    m <- models[[i]]
    trans_df$Gene[i] <- m$Gene
    trans_df$P_res[i] <- m$P_res
    
    counts <- unlist(m$residue_counts)
    counts <- counts[names(counts) != "-"]
    counts <- counts[names(counts) != "X"]
    if(length(counts)==0) next
    
    major <- names(which.max(counts))
    minors <- counts[counts > 10 & names(counts) != major]
    if(length(minors)==0) next
    
    trans_df$major[i] <- major
    trans_df$minor_list[i] <- paste(names(minors), collapse=",")
    trans_df$minor_counts[i] <- paste(minors, collapse=",")
  }
  
  all_trans[[f]] <- trans_df
}

trans <- do.call(rbind, all_trans)

threshold <- quantile(trans$P_res, 0.05)
trans$gwas_class <- ifelse(trans$P_res < threshold, "adaptive", "background")
table(trans$gwas_class)

long <- trans |>
  separate_rows(minor_list, minor_counts, sep=",") |>
  filter(minor_list != "", minor_counts != "", major != "") |>
  rename(minor = minor_list, minor_count = minor_counts) |>
  mutate(minor_count = as.numeric(minor_count))

transitions <- long |>
  group_by(gwas_class, major, minor) |>
  summarise(n = n(), .groups="drop")

hist(transitions$n)

mat_adap <- transitions |>
  filter(gwas_class=="adaptive") |>
  pivot_wider(names_from=minor, values_from=n, values_fill=0)

mat_bg <- transitions |>
  filter(gwas_class=="background") |>
  pivot_wider(names_from=minor, values_from=n, values_fill=0)

comp <- transitions |>
  group_by(major, minor) |>
  summarise(
    adaptive = sum(n[gwas_class=="adaptive"]),
    background = sum(n[gwas_class=="background"]),
    .groups="drop"
  ) |>
  mutate(
    total = adaptive + background,
    prop_adap = adaptive / total,
    prop_bg = background / total
  )

plot(comp$prop_adap, comp$prop_bg)

ggplot(transitions, aes(major, n, fill=gwas_class)) +
  geom_col(position="dodge") +
  facet_wrap(~minor) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(comp, aes(prop_bg, prop_adap, label=paste(major,minor,sep="→"))) +
  geom_point() +
  geom_text(nudge_y=0.01) +
  abline(0,1) +
  theme_bw()

transitions2 <- transitions |>
  filter(major!="X", minor!="X")

mat_adap <- transitions2 |>
  filter(gwas_class=="adaptive") |>
  complete(major, minor, fill=list(n=0))

mat_bg <- transitions2 |>
  filter(gwas_class=="background") |>
  complete(major, minor, fill=list(n=0))

mat_adap$gwas_class <- "adaptive"
mat_bg$gwas_class <- "background"

hm <- bind_rows(mat_adap, mat_bg)

ggplot(hm, aes(minor, major, fill=n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~gwas_class, nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1))

hm2 <- hm |>
  group_by(gwas_class) |>
  mutate(n_scaled = n / max(n)) |>
  ungroup()

ggplot(hm2, aes(minor, major, fill=n_scaled)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~gwas_class, nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1))

transitions <- long |>
  group_by(gwas_class, major, minor) |>
  summarise(
    total_minor = sum(minor_count),
    .groups="drop"
  )

comp <- transitions |>
  group_by(major, minor) |>
  summarise(
    adaptive = sum(total_minor[gwas_class=="adaptive"], na.rm=TRUE),
    background = sum(total_minor[gwas_class=="background"], na.rm=TRUE),
    .groups="drop"
  ) |>
  mutate(
    prop_adap = adaptive / sum(adaptive),
    prop_bg = background / sum(background)
  )

hist(comp$prop_adap)
hist(comp$prop_bg)


comp_total_bg <- sum(comp$background)
comp_total_ad <- sum(comp$adaptive)
comp$prop_bg <- comp$background / comp_total_bg
comp$prop_adap <- comp$adaptive / comp_total_ad
comp$transition <- paste0(comp$major, "->",comp$minor)
par(mfrow=c(1,1))
plot(comp$prop_adap, comp$prop_bg, main="Transitions, GWAS hits vs background")
text(comp$prop_adap, comp$prop_bg, labels=comp$transition, pos=3)
abline(0,1)

cor(comp$prop_adap, comp$prop_bg)
# ---- transition hydrophobicity ----

library(arrow)
library(data.table)
library(Peptides)

# Load AA properties
data(AAdata)
aa_hydrophobicity <- AAdata$Hydrophobicity$Tanford
aa_mw <- setNames(sapply(names(aa_hydrophobicity), mw), names(aa_hydrophobicity))
plot(aa_hydrophobicity, aa_mw)
cor(aa_hydrophobicity, aa_mw)
# Read all model files
model_files <- list.files("results/residue_models_clean/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

# Extract transitions from all sites
all_transitions <- list()

for (f in model_files) {
  models <- readRDS(f)
  
  for (i in seq_along(models)) {
    m <- models[[i]]
    counts <- unlist(m$residue_counts)
    counts <- counts[!names(counts) %in% c("-", "X")]
    if (length(counts) == 0) next
    
    major <- names(which.max(counts))
    minors <- counts[counts > 10 & names(counts) != major]
    if (length(minors) == 0) next
    
    for (minor in names(minors)) {
      all_transitions[[length(all_transitions) + 1]] <- data.frame(
        P_res = m$P_res,
        major = major,
        minor = minor,
        count = minors[minor],
        stringsAsFactors = FALSE
      )
    }
  }
}

trans_df <- do.call(rbind, all_transitions)
stopifnot(nrow(trans_df) > 0)

# Classify sites
threshold <- quantile(final_res$P_res, 0.05)
trans_df$gwas_class <- ifelse(trans_df$P_res < threshold, "adaptive", "background")

# Calculate property changes
trans_df$major_hydro <- aa_hydrophobicity[trans_df$major]
trans_df$minor_hydro <- aa_hydrophobicity[trans_df$minor]
trans_df$delta_hydro <- trans_df$minor_hydro - trans_df$major_hydro

trans_df$major_mw <- aa_mw[trans_df$major]
trans_df$minor_mw <- aa_mw[trans_df$minor]
trans_df$delta_mw <- trans_df$minor_mw - trans_df$major_mw

stopifnot(!any(is.na(trans_df$delta_hydro)))
stopifnot(!any(is.na(trans_df$delta_mw)))

# Classify direction of changes
trans_df$hydro_dir <- ifelse(trans_df$delta_hydro > 0, "increase", 
                             ifelse(trans_df$delta_hydro < 0, "decrease", "no_change"))
trans_df$mw_dir <- ifelse(trans_df$delta_mw > 0, "increase",
                          ifelse(trans_df$delta_mw < 0, "decrease", "no_change"))
hist(trans_df$delta_hydro[trans_df$gwas_class=="background"])
hist(trans_df$delta_hydro[trans_df$gwas_class=="adaptive"])
boxplot(delta_hydro ~ gwas_class, trans_df)
ggplot(trans_df, aes(gwas_class, delta_hydro)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=.15, outlier.shape=NA) +
  geom_jitter(width=.1, alpha=.3)

# Count transitions by direction
hydro_counts <- aggregate(count ~ gwas_class + hydro_dir, trans_df, sum)
mw_counts <- aggregate(count ~ gwas_class + mw_dir, trans_df, sum)

# Calculate proportions
hydro_counts$prop <- ave(hydro_counts$count, hydro_counts$gwas_class, 
                         FUN = function(x) x / sum(x))
mw_counts$prop <- ave(mw_counts$count, mw_counts$gwas_class,
                      FUN = function(x) x / sum(x))

ggplot(trans_df, aes(x = delta_hydro, weight = count, fill = gwas_class)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  ggtitle("Difference in hydrophobicity in transitions")

ggplot(trans_df, aes(x = gwas_class, y = delta_hydro, weight = count)) +
  geom_boxplot() +
  theme_minimal()
# Plot
par(mfrow = c(1, 2))

barplot(prop ~ gwas_class + hydro_dir, 
        data = hydro_counts[hydro_counts$hydro_dir != "no_change", ],
        beside = TRUE,
        col = c("#EF4444", "#9CA3AF"),
        las = 2,
        main = "Hydrophobicity changes",
        ylab = "Proportion of transitions",
        legend.text = c("Adaptive", "Background"),
        args.legend = list(x = "bottom"))
abline(h=0.5, col="black", lty=2)

barplot(prop ~ gwas_class + mw_dir,
        data = mw_counts,
        beside = TRUE,
        col = c("#EF4444", "#9CA3AF"),
        las = 2,
        ylim=c(0,0.6),
        main = "Molecular weight changes",
        ylab = "Proportion of transitions")
abline(h=0.5, col="black", lty=2)
# Statistical tests
cat("\n=== Hydrophobicity test ===\n")
hydro_mat <- matrix(c(
  sum(trans_df$count[trans_df$gwas_class == "adaptive" & trans_df$delta_hydro > 0]),
  sum(trans_df$count[trans_df$gwas_class == "adaptive" & trans_df$delta_hydro < 0]),
  sum(trans_df$count[trans_df$gwas_class == "background" & trans_df$delta_hydro > 0]),
  sum(trans_df$count[trans_df$gwas_class == "background" & trans_df$delta_hydro < 0])
), nrow = 2, byrow = TRUE)
rownames(hydro_mat) <- c("adaptive", "background")
colnames(hydro_mat) <- c("increase", "decrease")
print(hydro_mat)
print(chisq.test(hydro_mat))

cat("\n=== Molecular weight test ===\n")
mw_mat <- matrix(c(
  sum(trans_df$count[trans_df$gwas_class == "adaptive" & trans_df$delta_mw > 0]),
  sum(trans_df$count[trans_df$gwas_class == "adaptive" & trans_df$delta_mw < 0]),
  sum(trans_df$count[trans_df$gwas_class == "background" & trans_df$delta_mw > 0]),
  sum(trans_df$count[trans_df$gwas_class == "background" & trans_df$delta_mw < 0])
), nrow = 2, byrow = TRUE)
rownames(mw_mat) <- c("adaptive", "background")
colnames(mw_mat) <- c("increase", "decrease")
print(mw_mat)
print(chisq.test(mw_mat))


# ---- AA prop characteristic analysis. ----

library(Peptides)
data(AAdata)
aa_hydrophobicity <- AAdata$Hydrophobicity$Tanford
aa_lipophilicity <- AAdata$zScales$Z1

comp <- comp[comp$minor!="X",]
# Add amino acid properties to transitions data
comp <- comp |>
  mutate(
    major_hydro = aa_hydrophobicity[major],
    minor_hydro = aa_hydrophobicity[minor],
    delta_hydro = minor_hydro - major_hydro,
    
    major_lipo = aa_lipophilicity[major],
    minor_lipo = aa_lipophilicity[minor],
    delta_lipo = minor_lipo - major_lipo,
    
    major_mw = mw(major),
    minor_mw = mw(minor),
    delta_mw = minor_mw - major_mw
  )

stopifnot(!any(is.na(comp$delta_hydro)))
stopifnot(!any(is.na(comp$delta_lipo)))
stopifnot(!any(is.na(comp$delta_mw)))

transitions <- transitions[transitions$minor!="X",]

# Test for differences between adaptive vs background
delta_df <- transitions |>
  mutate(
    major_hydro = aa_hydrophobicity[major],
    minor_hydro = aa_hydrophobicity[minor],
    delta_hydro = minor_hydro - major_hydro,
    
    major_lipo = aa_lipophilicity[major],
    minor_lipo = aa_lipophilicity[minor],
    delta_lipo = minor_lipo - major_lipo,
    
    major_mw = mw(major),
    minor_mw = mw(minor),
    delta_mw = minor_mw - major_mw
  )

# Weighted t-tests using total_minor as weights
t.test(delta_hydro ~ gwas_class, data = delta_df, var.equal = FALSE)
t.test(delta_lipo ~ gwas_class, data = delta_df, var.equal = FALSE)
t.test(delta_mw ~ gwas_class, data = delta_df, var.equal = FALSE)

#but should 
# Classify transitions by direction of change
comp <- comp |>
  mutate(
    hydro_direction = case_when(
      delta_hydro > 0 ~ "increase",
      delta_hydro < 0 ~ "decrease",
      TRUE ~ "no_change"
    ),
    lipo_direction = case_when(
      delta_lipo > 0 ~ "increase",
      delta_lipo < 0 ~ "decrease",
      TRUE ~ "no_change"
    ),
    mw_direction = case_when(
      delta_mw > 0 ~ "increase",
      delta_mw < 0 ~ "decrease",
      TRUE ~ "no_change"
    )
  )

# Calculate proportions for hydrophobicity
hydro_summary <- comp |>
  group_by(hydro_direction) |>
  summarise(
    adaptive = sum(adaptive),
    background = sum(background),
    .groups = "drop"
  ) |>
  mutate(
    prop_adaptive = adaptive / sum(adaptive),
    prop_background = background / sum(background)
  )

print("Hydrophobicity transitions:")
print(hydro_summary)

# Test if proportion of increasing transitions differs
cont_table_hydro <- matrix(c(
  sum(comp$adaptive[comp$delta_hydro > 0]),
  sum(comp$adaptive[comp$delta_hydro < 0]),
  sum(comp$background[comp$delta_hydro > 0]),
  sum(comp$background[comp$delta_hydro < 0])
), nrow = 2, byrow = TRUE)
rownames(cont_table_hydro) <- c("adaptive", "background")
colnames(cont_table_hydro) <- c("increase", "decrease")

print(cont_table_hydro)
chisq.test(cont_table_hydro)

# Repeat for lipophilicity
lipo_summary <- comp |>
  group_by(lipo_direction) |>
  summarise(
    adaptive = sum(adaptive),
    background = sum(background),
    .groups = "drop"
  ) |>
  mutate(
    prop_adaptive = adaptive / sum(adaptive),
    prop_background = background / sum(background)
  )

print("Lipophilicity transitions:")
print(lipo_summary)

# Repeat for molecular weight
mw_summary <- comp |>
  group_by(mw_direction) |>
  summarise(
    adaptive = sum(adaptive),
    background = sum(background),
    .groups = "drop"
  ) |>
  mutate(
    prop_adaptive = adaptive / sum(adaptive),
    prop_background = background / sum(background)
  )

print("Molecular weight transitions:")
print(mw_summary)

# Combine summaries for plotting
prop_summary <- bind_rows(
  hydro_summary |> 
    select(direction = hydro_direction, prop_adaptive, prop_background) |>
    mutate(property = "Hydrophobicity"),
  lipo_summary |> 
    select(direction = lipo_direction, prop_adaptive, prop_background) |>
    mutate(property = "Lipophilicity"),
  mw_summary |> 
    select(direction = mw_direction, prop_adaptive, prop_background) |>
    mutate(property = "Molecular Weight")
) |>
  pivot_longer(
    cols = c(prop_adaptive, prop_background),
    names_to = "class",
    values_to = "proportion",
    names_prefix = "prop_"
  ) |>
  filter(direction != "no_change")

# Bar plot with error bars
ggplot(prop_summary_ci, aes(x = direction, y = proportion, fill = class)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  facet_wrap(~property) +
  scale_fill_manual(
    values = c("adaptive" = "#EF4444", "background" = "#9CA3AF"),
    labels = c("Adaptive", "Background")
  ) +
  theme_bw() +
  labs(
    title = "Direction of amino acid property changes",
    x = "Direction of change",
    y = "Proportion of transitions",
    fill = "Site class"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Visualize distributions
par(mfrow = c(3, 2))
hist(delta_df$delta_hydro[delta_df$gwas_class == "adaptive"], 
     main = "Δ Hydrophobicity (adaptive)", xlab = "")
hist(delta_df$delta_hydro[delta_df$gwas_class == "background"], 
     main = "Δ Hydrophobicity (background)", xlab = "")
hist(delta_df$delta_lipo[delta_df$gwas_class == "adaptive"], 
     main = "Δ Lipophilicity (adaptive)", xlab = "")
hist(delta_df$delta_lipo[delta_df$gwas_class == "background"], 
     main = "Δ Lipophilicity (background)", xlab = "")
hist(delta_df$delta_mw[delta_df$gwas_class == "adaptive"], 
     main = "Δ Molecular weight (adaptive)", xlab = "")
hist(delta_df$delta_mw[delta_df$gwas_class == "background"], 
     main = "Δ Molecular weight (background)", xlab = "")

# Summary statistics
aggregate(cbind(delta_hydro, delta_lipo, delta_mw) ~ gwas_class, 
          data = delta_df, FUN = mean)

# Add amino acid properties to df_filtered
df_filtered <- df_filtered |>
  mutate(
    hydro = aa_hydrophobicity[AA],
    lipo = aa_lipophilicity[AA],
    mw = mw(AA)
  )

stopifnot(!any(is.na(df_filtered$hydro)))
stopifnot(!any(is.na(df_filtered$lipo)))
stopifnot(!any(is.na(df_filtered$mw)))

# Calculate weighted means by effect class and allele class
aa_props <- df_filtered |>
  group_by(effect_class, allele_class) |>
  summarise(
    mean_hydro = weighted.mean(hydro, count),
    mean_lipo = weighted.mean(lipo, count),
    mean_mw = weighted.mean(mw, count),
    .groups = "drop"
  )

print(aa_props)

# Test differences in minor alleles across effect classes
minor_only <- df_filtered |>
  filter(allele_class == "minor")

# Pairwise comparisons for minor alleles
t.test(hydro ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "cold"), ])
t.test(lipo ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "cold"), ])
t.test(mw ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "cold"), ])

t.test(hydro ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "no_effect"), ])
t.test(lipo ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "no_effect"), ])
t.test(mw ~ effect_class, data = minor_only[minor_only$effect_class %in% c("hot", "no_effect"), ])

t.test(hydro ~ effect_class, data = minor_only[minor_only$effect_class %in% c("cold", "no_effect"), ])
t.test(lipo ~ effect_class, data = minor_only[minor_only$effect_class %in% c("cold", "no_effect"), ])
t.test(mw ~ effect_class, data = minor_only[minor_only$effect_class %in% c("cold", "no_effect"), ])


plot(aa_hydrophobicity, aa_lipophilicity)
# Visualize distributions of properties by effect class for minor alleles
par(mfrow = c(3, 1))
boxplot(hydro ~ effect_class, data = minor_only, 
        main = "Hydrophobicity (minor alleles)", ylab = "Hydrophobicity")
boxplot(lipo ~ effect_class, data = minor_only, 
        main = "Lipophilicity (minor alleles)", ylab = "Lipophilicity")
boxplot(mw ~ effect_class, data = minor_only, 
        main = "Molecular weight (minor alleles)", ylab = "MW")

# ---- positional enrichment ---- 

# Extract position data and test for positional enrichment
all_res <- list()
for(f in model_files){
  models <- readRDS(f)
  res <- data.frame(
    Aligned_Position = integer(length(models)),
    P_res = numeric(length(models)),
    Effect = numeric(length(models)),
    Gene = character(length(models)),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_along(models)){
    mod <- models[[i]]
    res$Aligned_Position[i] <- mod$Aligned_Position
    res$P_res[i] <- mod$P_res
    effs <- mod$effects
    res$Effect[i] <- mean(effs$Effect)
    res$Gene[i] <- mod$Gene
  }
  
  all_res[[f]] <- res
}

# Bind and classify
final_res <- do.call(rbind, all_res)
stopifnot(nrow(final_res) > 0)
stopifnot(all(c("Aligned_Position", "P_res", "Gene") %in% names(final_res)))

threshold <- quantile(final_res$P_res, 0.05)
final_res$gwas_class <- ifelse(final_res$P_res < threshold, "adaptive", "background")

# Calculate normalized position (0-1) within each gene
final_res <- final_res |>
  group_by(Gene) |>
  mutate(
    min_pos = min(Aligned_Position),
    max_pos = max(Aligned_Position),
    gene_length = max_pos - min_pos,
    norm_position = (Aligned_Position - min_pos) / gene_length
  ) |>
  ungroup()

stopifnot(all(final_res$norm_position >= 0 & final_res$norm_position <= 1))
stopifnot(all(final_res$gene_length >= 0))

# Distribution comparison
par(mfrow=c(1,2))
hist(final_res$norm_position[final_res$gwas_class=="adaptive"], 
     breaks=20, main="Adaptive sites", xlab="Normalized position", col="#EF4444")
hist(final_res$norm_position[final_res$gwas_class=="background"], 
     breaks=20, main="Background sites", xlab="Normalized position", col="#9CA3AF")

# Statistical test
ks_test <- ks.test(
  final_res$norm_position[final_res$gwas_class=="adaptive"],
  final_res$norm_position[final_res$gwas_class=="background"]
)
print(ks_test)

contingency_table <- table(final_res$pos_bin, final_res$gwas_class)
chi_test <- chisq.test(contingency_table)
print(chi_test)

# Summary stats by class
summary_stats <- final_res |>
  group_by(gwas_class) |>
  summarise(
    n = n(),
    mean_pos = mean(norm_position),
    median_pos = median(norm_position),
    sd_pos = sd(norm_position)
  )
print(summary_stats)

# Calculate binned position enrichment with error bars
n_bins <- 20
final_res$pos_bin <- cut(final_res$norm_position, 
                         breaks = seq(0, 1, length.out = n_bins + 1),
                         include.lowest = TRUE)

bin_summary <- final_res |>
  group_by(pos_bin, gwas_class) |>
  summarise(
    n = n(),
    mean_pos = mean(norm_position),
    .groups = "drop"
  ) |>
  group_by(gwas_class) |>
  mutate(
    total_n = sum(n),
    proportion = n / total_n,
    se = sqrt(proportion * (1 - proportion) / total_n)
  ) |>
  ungroup()

stopifnot(all(bin_summary$proportion >= 0 & bin_summary$proportion <= 1))

# Line plot with error bars
ggplot(bin_summary, aes(x = mean_pos, y = proportion, color = gwas_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = proportion - se, ymax = proportion + se), 
                width = 0.02, alpha = 0.7) +
  scale_color_manual(values = c(
    "adaptive" = "#EF4444",
    "background" = "#9CA3AF"
  )) +
  theme_bw() +
  labs(
    x = "Normalized gene position",
    y = "Proportion of sites",
    color = "GWAS class",
    title = "Positional distribution of adaptive vs background sites"
  )

# Calculate binned position enrichment with error bars, bifurcating adaptive by effect
n_bins <- 20
final_res$pos_bin <- cut(final_res$norm_position, 
                         breaks = seq(0, 1, length.out = n_bins + 1),
                         include.lowest = TRUE)

# Add effect-based classification for adaptive sites
stopifnot("Effect" %in% names(final_res))
final_res <- final_res |>
  mutate(
    effect_class = case_when(
      gwas_class == "background" ~ "background",
      gwas_class == "adaptive" & Effect > 0 ~ "hot",
      gwas_class == "adaptive" & Effect <= 0 ~ "cold",
      TRUE ~ NA_character_
    )
  )

stopifnot(!any(is.na(final_res$effect_class)))

bin_summary <- final_res |>
  group_by(pos_bin, effect_class) |>
  summarise(
    n = n(),
    mean_pos = mean(norm_position),
    .groups = "drop"
  ) |>
  group_by(effect_class) |>
  mutate(
    total_n = sum(n),
    proportion = n / total_n,
    se = sqrt(proportion * (1 - proportion) / total_n)
  ) |>
  ungroup()

stopifnot(all(bin_summary$proportion >= 0 & bin_summary$proportion <= 1))

# Line plot with error bars
ggplot(bin_summary, aes(x = mean_pos, y = proportion, color = effect_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = proportion - se, ymax = proportion + se), 
                width = 0.02, alpha = 0.7) +
  scale_color_manual(values = c(
    "hot" = "#EF4444",
    "cold" = "#3B82F6",
    "background" = "#9CA3AF"
  )) +
  theme_bw() +
  labs(
    x = "Normalized gene position",
    y = "Proportion of sites",
    color = "Site class",
    title = "Positional distribution: hot/cold adaptive vs background sites"
  )


# ---- hit classification ----
# goal is to classify each hit into general, cold, or hot 
all_res <- list()

extract_all_predictors <- function(model_files) {
  stopifnot(all(file.exists(model_files)))
  
  results <- list()
  
  for (f in model_files) {
    models <- readRDS(f)
    
    for (i in seq_along(models)) {
      mod <- models[[i]]
      
      # Skip if no effects
      if (is.null(mod$effects) || nrow(mod$effects) == 0) next
      
      effs <- mod$effects
      
      # Clean residue names (remove X_aa prefix)
      residues <- sub("^X_aa", "", effs$Residue)
      
      # Get counts for each residue
      counts <- sapply(residues, function(r) {
        if (r %in% names(mod$residue_counts)) {
          mod$residue_counts[[r]]
        } else {
          NA_integer_
        }
      })
      
      # Create one row per predictor
      for (j in seq_len(nrow(effs))) {
        results[[length(results) + 1]] <- data.frame(
          Gene = mod$Gene,
          Position = mod$Aligned_Position,
          AA = residues[j],
          Count = counts[j],
          Effect = effs$Effect[j],
          SE = effs$SE[j],
          P_value = effs$P_value[j],
          P_res = mod$P_res,
          N = mod$N,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  stopifnot(length(results) > 0)
  do.call(rbind, results)
}

# Usage:
all_predictors <- extract_all_predictors(model_files)
hit_predictors <- all_predictors[all_predictors$P_res < quantile(final_res$P_res,0.05),]
hit_predictors <- hit_predictors[hit_predictors$Count > 10,]
hit_predictors <- hit_predictors[hit_predictors$AA != "-",]
hist(hit_predictors$Effect)
hist(-log10(hit_predictors$P_value))
plot(hit_predictors$Effect,-log10(hit_predictors$P_value)
     ,main="Residue effect by significance")

major_predictors <- hit_predictors[hit_predictors$Count > 4000,]
minor_predictors <- hit_predictors[hit_predictors$Count < 4000,]
minor_predictors <- minor_predictors[minor_predictors$AA != "-",]

par(mfrow=c(1,2))
plot(major_predictors$Effect,-log10(major_predictors$P_value)
     ,main="Major Residue effect by significance")
plot(minor_predictors$Effect,-log10(minor_predictors$P_value)
     ,main="Minor Residue effect by significance")
minor_predictors$AA

aa_colors <- c(
  # Hydrophobic (browns/oranges)
  "A" = "#D2691E", "V" = "#CD853F", "I" = "#A0522D", "L" = "#8B4513", 
  "M" = "#DEB887", "F" = "#F4A460", "W" = "#BC8F8F", "P" = "#DAA520",
  # Polar uncharged (greens)
  "S" = "#90EE90", "T" = "#3CB371", "C" = "#2E8B57", "Y" = "#66CDAA",
  "N" = "#20B2AA", "Q" = "#48D1CC",
  # Positively charged (blues)
  "K" = "#4169E1", "R" = "#1E90FF", "H" = "#6495ED",
  # Negatively charged (reds)
  "D" = "#DC143C", "E" = "#B22222",
  # Special
  "G" = "#FFD700", "-" = "#808080"
)

hit_predictors$color <- aa_colors[hit_predictors$AA]
major_predictors$color <- aa_colors[major_predictors$AA]
minor_predictors$color <- aa_colors[minor_predictors$AA]

# Side-by-side comparison
par(mfrow = c(1, 2))

plot(major_predictors$Effect, -log10(major_predictors$P_value),
     col = major_predictors$color, pch = 19, cex = 1,
     main = "Major Residue Effect",
     xlab = "Effect", ylab = "-log10(P-value)")
unique_aa_major <- sort(unique(major_predictors$AA))
legend("topright", legend = unique_aa_major,
       col = aa_colors[unique_aa_major], pch = 19,
       cex = 0.6, bg = "white")

plot(minor_predictors$Effect, -log10(minor_predictors$P_value),
     col = minor_predictors$color, pch = 19, cex = 0.8,
     main = "Minor Residue Effect",
     xlab = "Effect", ylab = "-log10(P-value)")
unique_aa_minor <- sort(unique(minor_predictors$AA))
legend("topright", legend = unique_aa_minor,
       col = aa_colors[unique_aa_minor], pch = 19,
       ncol = 2, cex = 0.6, bg = "white")

par(mfrow = c(1, 1))

# Example usage:
classifications <- classify_all_sites(model_files)
# table(classifications$Classification)
# mixed_sites <- classifications[classifications$Classification == "mixed", ]
table(major_predictors$AA)




plot(to_enrich[aa_cols], from_enrich[aa_cols])
x <- to_enrich[aa_cols]
y <- from_enrich[aa_cols]

plot(as.numeric(x), as.numeric(y))

text(as.numeric(x), as.numeric(y), labels = names(x), pos = 3)
abline(0,1)


from_prop <- from_counts / sum(from_counts)
from_enrich <- from_prop / bg_prop[names(from_prop)]

barplot(from_enrich, las=2, main="FROM enrichment")
sort(from_enrich, decreasing=TRUE)

to_prop <- to_counts / sum(to_counts)
to_enrich <- to_prop / bg_prop[names(to_prop)]

barplot(to_enrich, las=2, main="Minor allele enrichment")
sort(to_enrich, decreasing=TRUE)


split_FROM <- split(FROM, FROM$complex)
split_TO   <- split(TO,   TO$complex)

for(cmp in names(split_FROM)){
  x <- table(split_FROM[[cmp]]$AA)
  barplot(x, las=2, main=paste(cmp, "FROM"))
}

for(cmp in names(split_TO)){
  x <- table(split_TO[[cmp]]$AA)
  barplot(x, las=2, main=paste(cmp, "TO"))
}


#model SUM SUMMARY OF EFFECT BY major, minor, and transition
#for a given site, is a biallelic model appropriate? Is it always the case that the 
#Minor allele affect is the same for all minor alleles
#also record site level effects, just so I can plot them 





example_model <- readRDS(model_files[1])
plot(example_model)

# Extract summary stats from each gene
all_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      N = m$N,
      R2_reduced = m$R2_reduced,
      R2_full = m$R2_full,
      Delta_R2 = m$R2_full - m$R2_reduced,
      P = m$P_res
    )
  }))
}))

stopifnot(nrow(all_results) > 0)

# Calculate -log10(P)
all_results[, neg_log10_P := -log10(P)]
all_results[is.infinite(neg_log10_P), neg_log10_P := NA]

# Create output directory
dir.create("results/manhattan_plots", showWarnings = FALSE, recursive = TRUE)

# Plot for each gene
genes <- unique(all_results$Gene)
bf <- 0.05/nrow(all_results) 

for (gene in genes) {
  gene_data <- all_results[Gene == gene]
  
  p <- ggplot(gene_data, aes(x = Position, y = neg_log10_P)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05/12000), linetype = "dashed", color = "red") +
    labs(
      title = paste0("Manhattan Plot: ", gene),
      x = "Aligned Position",
      y = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  print(p)
  
  #ggsave(
  #  filename = paste0("results/manhattan_plots/", gene, "_manhattan.png"),
  #  plot = p, width = 10, height = 6, dpi = 300)
}



# ---- correlate cleaned results with uncleaned results ----


res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)
f <- res_files[1]

res_list <- readRDS(f)
gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
res_df <- do.call(rbind, res_list)
unclean_res_df <- na.omit(res_df)

clean_res_df <- all_results[grep("atpA", all_results$Gene),]
par(mfrow=c(1,2))
plot(-log10(unclean_res_df$P_res))
plot(-log10(clean_res_df$P))
dev.off()
res <- merge(clean_res_df[,c("Position", "P")],
             unclean_res_df[,c("Aligned_Position", "P_res")],
            by.x="Position", by.y="Aligned_Position")
plot(-log10(res$P_res), -log10(res$P),
     main="Effect of removing outlying genes by MDS on GWAS results",
     xlab="Length Filter -log10(p)",
     ylab="Manual embedding MDS + Length Filter -log10(p)")
abline(a=0,b=1,col="tomato", lty=2)
text(5,1,paste("cor=",
               round(cor(-log10(res$P_res), -log10(res$P)),3)))


res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

pdf("plots/mds_control_gwas_effect.pdf", width=10, height=10)
par(mfrow=c(5,4))
for (f in res_files) {
  res_list <- readRDS(f)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  res_df <- do.call(rbind, res_list)
  unclean_res_df <- na.omit(res_df)
  
  clean_res_df <- all_results[grep(gene_name, all_results$Gene),]
  if(nrow(clean_res_df) > 1) {
  
  res <- merge(clean_res_df[,c("Position", "P")],
               unclean_res_df[,c("Aligned_Position", "P_res")],
               by.x="Position", by.y="Aligned_Position")
  plot(-log10(res$P_res), -log10(res$P),
       main=
         paste(gene_name,
                "Effect of removing outlying genes by MDS on GWAS results"),
       xlab="Length Filter -log10(p)",
       ylab="Manual embedding MDS + Length Filter -log10(p)")
  abline(a=0,b=1,col="tomato", lty=2)
  text(quantile(-log10(res$P_res), 0.99),
       quantile(-log10(res$P), 0.5),paste("pearson cor=",
                 round(cor(-log10(res$P_res), -log10(res$P)),3)))
  text(quantile(-log10(res$P_res), 0.99),
       quantile(-log10(res$P), 0.3),
       paste("spearman cor=",
             round(cor(-log10(res$P_res), -log10(res$P), method="spearman"),3)))
  }
}
dev.off()

# ---- manhattan ----

# Sort and assign cumulative positions
all_results <- all_results[order(Gene, Position)]

gene_lengths <- all_results[, .(len = max(Position, na.rm = TRUE)), by = Gene]
gene_offsets <- c(0, cumsum(head(gene_lengths$len, -1)))
names(gene_offsets) <- gene_lengths$Gene

all_results[, pos_global := Position + gene_offsets[Gene]]

# Bonferroni threshold
bf <- 0.05 / nrow(all_results)
ylim_max <- 1 + max(all_results$neg_log10_P, na.rm = TRUE)

# Alternate colors for genes
gene_colors <- setNames(rep(c("steelblue3", "grey60"),
                            length.out = length(unique(all_results$Gene))),
                        unique(all_results$Gene))

# Plot
plot(all_results$pos_global, all_results$neg_log10_P,
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "-log10(p)",
     ylim = c(0, ylim_max),
     main = "Chloroplast Protein Association Study",
     xaxt = "n")

# Draw points per gene
for (g in unique(all_results$Gene)) {
  sub <- all_results[Gene == g]
  points(sub$pos_global, sub$neg_log10_P, col = gene_colors[g], lwd = 1.5)
  # Highlight significant points
  sig <- sub[neg_log10_P > -log10(bf)]
  if (nrow(sig) > 0) {
    points(sig$pos_global, sig$neg_log10_P, col = gene_colors[g], pch = 16)
  }
}

# Bonferroni line
abline(h = -log10(bf), col = "black", lty = 2)

# Gene labels
midpoints <- gene_lengths[, .(mid = gene_offsets[Gene] + len / 2), by = Gene]
axis(1, at = midpoints$mid, labels = midpoints$Gene, las = 2, cex.axis = 0.7)




# ---- clean/dirty comparison ----


res_files <- list.files("results/", pattern = "tmp_results_nested_.*\\.rds$", full.names = TRUE)

f <- res_files[1]
res <- readRDS(f)

#list of residue level result objects
str(res)

res[1]


all_comparisons <- rbindlist(lapply(res_files, function(f) {
  res_list <- readRDS(f)
  gene_name <- sub("tmp_results_nested_(.*)\\.rds$", "\\1", basename(f))
  unclean_res_df <- na.omit(do.call(rbind, res_list))
  clean_res_df <- all_results[grep(gene_name, Gene)]
  
  if (nrow(clean_res_df) > 1) {
    merge(clean_res_df[, .(Position, P_clean = P)],
          unclean_res_df[, .(Aligned_Position, P_unclean = P_res)],
          by.x = "Position", by.y = "Aligned_Position")
  }
}))

pearson_cor <- cor(-log10(all_comparisons$P_unclean), -log10(all_comparisons$P_clean))
spearman_cor <- cor(-log10(all_comparisons$P_unclean), -log10(all_comparisons$P_clean), method = "spearman")

plot(-log10(all_comparisons$P_unclean), -log10(all_comparisons$P_clean),
     xlab = "Length filter -log10(P)", ylab = "Length + Manual ESMC MDS Z-score filter -log10(P)",
     main = "Effect of MDS Filtering on GWAS (All Genes)",
     pch = 16, col = rgb(0, 0, 0, 0.3))
abline(a = 0, b = 1, col = "blue", lty = 3)

segments(-log10(bf),0,-log10(bf),-log10(bf), col="coral", lty=2)
segments(0,-log10(bf),-log10(bf),-log10(bf), col="coral", lty=2)
legend("topleft", 
       legend = c(paste("Pearson:", round(pearson_cor, 3)),
                  paste("Spearman:", round(spearman_cor, 3)),
                  paste("Bonferroni:", round(-log10(bf), 2))),
       lty = c(NA, NA, 2),
       col = c(NA, NA, "coral"),
       pch = c(NA, NA, NA),
       bty = "n")


LSD::comparisonplot(
  -log10(all_comparisons$P_unclean),
  -log10(all_comparisons$P_clean),
  pimp = T)

# ---- manhattan-hudson plot ----  

# Extract effects with their signs from model files
all_results_signed <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (!is.null(m$effects) && nrow(m$effects) > 0) {
      # Get the dominant effect sign (could use mean, median, or max absolute)
      mean_effect <- mean(m$effects$Effect, na.rm = TRUE)
      
      data.table(
        Gene = m$Gene,
        Position = m$Aligned_Position,
        N = m$N,
        R2_reduced = m$R2_reduced,
        R2_full = m$R2_full,
        Delta_R2 = m$R2_full - m$R2_reduced,
        P = m$P_res,
        Effect_sign = sign(mean_effect)
      )
    }
  }))
}))

table(all_results_signed$Effect_sign)

# Calculate signed -log10(P)
all_results_signed[, signed_log10_P := -log10(P) * Effect_sign]
all_results_signed[is.infinite(signed_log10_P), signed_log10_P := NA]

# Sort and assign cumulative positions
all_results_signed <- all_results_signed[order(Gene, Position)]

gene_lengths <- all_results_signed[, .(len = max(Position, na.rm = TRUE)), by = Gene]
gene_offsets <- c(0, cumsum(head(gene_lengths$len, -1)))
names(gene_offsets) <- gene_lengths$Gene

all_results_signed[, pos_global := Position + gene_offsets[Gene]]

# Bonferroni threshold
bf <- 0.05 / nrow(all_results_signed)
ylim_max <- 1 + max(abs(all_results_signed$signed_log10_P), na.rm = TRUE)
ylim_min <- 1 - min(all_results_signed$signed_log10_P, na.rm = TRUE)

# devtools::install_github("AndiKur4/MaizePal")
library("MaizePal")
maize_pal("Painted")
# Alternate colors for genes
gene_colors <- setNames(rep(maize_pal("OaxacaGreen"),
                            length.out = length(unique(all_results_signed$Gene))),
                        unique(all_results_signed$Gene))

# Plot
plot(all_results_signed$pos_global, all_results_signed$signed_log10_P,
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "Signed -log10(p) [+ heat adaptive, - cold adaptive]",
     ylim = c(-ylim_min, ylim_max),
     main = "Manhattan-Hudson Plot: Chloroplast Protein Temperature GWAS",
     xaxt = "n")

# Add horizontal line at 0
abline(h = 0, col = "grey50", lty = 1)

# Draw points per gene
for (g in unique(all_results_signed$Gene)) {
  sub <- all_results_signed[Gene == g]
  points(sub$pos_global, sub$signed_log10_P, col = gene_colors[g], lwd = 1.5)
  # Highlight significant points
  sig <- sub[abs(signed_log10_P) > -log10(bf)]
  if (nrow(sig) > 0) {
    points(sig$pos_global, sig$signed_log10_P, col = gene_colors[g], pch = 16)
  }
}

# Bonferroni lines (positive and negative)
abline(h = -log10(bf), col = "coral", lty = 2)
abline(h = -(-log10(bf)), col = "coral", lty = 2)

# Gene labels
midpoints <- gene_lengths[, .(mid = gene_offsets[Gene] + len / 2), by = Gene]
axis(1, at = midpoints$mid, labels = midpoints$Gene, las = 2, cex.axis = 0.7)


top <- all_results[all_results$P < quantile(all_results$P,0.05),]






library(Peptides)
library(data.table)

aa_props <- as.data.table(AAdata)
setnames(aa_props, "AA", "Residue")

sig_residues <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (m$P_res < bf) {
      counts <- as.data.table(m$residue_counts)
      setnames(counts, c("Residue", "Count"))
      counts[, `:=`(Gene = m$Gene, Position = m$Aligned_Position, 
                    P = m$P_res, N = m$N)]
      counts
    }
  }))
}))

stopifnot(nrow(sig_residues) > 0)

# Calculate proportions
sig_residues[, Prop := Count / N]

# Merge with AA properties
sig_residues <- merge(sig_residues, aa_props, by = "Residue", all.x = TRUE)

# Summary by position
pos_summary <- sig_residues[Residue != "-", .(
  Mean_Mass = weighted.mean(MW, Prop, na.rm = TRUE),
  Mean_pI = weighted.mean(pI, Prop, na.rm = TRUE),
  Mean_Hydrophobicity = weighted.mean(Hydrophobicity, Prop, na.rm = TRUE),
  Prop_Charged = sum(Prop[Charge != 0], na.rm = TRUE),
  Prop_Polar = sum(Prop[Polarity == "polar"], na.rm = TRUE),
  Prop_Aromatic = sum(Prop[Residue %in% c("F", "Y", "W")], na.rm = TRUE)
), by = .(Gene, Position, P)]

# Summary across all significant positions
overall_summary <- sig_residues[Residue != "-", .(
  Count = sum(Count),
  Total_N = sum(Count)
), by = Residue][, Freq := Count / Total_N]

overall_summary <- merge(overall_summary, aa_props, by = "Residue", all.x = TRUE)
setorder(overall_summary, -Freq)

print(overall_summary)
print(pos_summary)
