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
  
  out_path <- paste0("results/residue_models_triple/", gene, "_effects.rds")
  
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
    
    res_factor <- factor(residues)
    X_aa <- model.matrix(~ res_factor - 1)
    
    # Remove zero variance columns
    var0 <- which(apply(X_aa, 2, var) == 0)
    if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
    if (ncol(X_aa) == 0) next
    
    # -----------------------------------------------------------------
    # FIT MODELS
    # -----------------------------------------------------------------
    
    # Model 0: intercept only (null model)
    fit_null <- tryCatch(lm(y ~ 1), error = function(e) NULL)
    if (is.null(fit_null)) next
    
    # Model 1: residue variation only (no pop structure control)
    fit_aa <- tryCatch(lm(y ~ X_aa), error = function(e) NULL)
    if (is.null(fit_aa)) next
    
    # Model 2: population structure control only
    fit_pcs <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
    if (is.null(fit_pcs)) next
    
    # Model 3: residue variation with population structure control
    fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
    if (is.null(fit_full)) next
    
    r2_aa <- summary(fit_aa)$r.squared
    r2_pcs <- summary(fit_pcs)$r.squared
    r2_full <- summary(fit_full)$r.squared
    
    # P-values for each test
    p_aa_only <- tryCatch(anova(fit_null, fit_aa)[2, "Pr(>F)"], error = function(e) NA_real_)
    p_pcs_only <- tryCatch(anova(fit_null, fit_pcs)[2, "Pr(>F)"], error = function(e) NA_real_)
    p_aa_with_pcs <- tryCatch(anova(fit_pcs, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
    
    # Extract only residue effect sizes and SEs from full model
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
      R2_aa = r2_aa,
      R2_pcs = r2_pcs,
      R2_full = r2_full,
      P_aa_only = p_aa_only,
      P_pcs_only = p_pcs_only,
      P_aa_with_pcs = p_aa_with_pcs,
      residue_counts = as.list(residue_table),
      effects = effects_dt,
      IDs = common_ids
    )
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) > 0) {
    dir.create("results/residue_models_triple", showWarnings = FALSE, recursive = TRUE)
    saveRDS(results_list, out_path)
    message("\nCompleted ", gene, ": ", length(results_list), " positions with effects")
  } else {
    message("No positions retained for gene ", gene)
  }
}

# ---- manhattans & scattering ----
library(data.table)

model_files <- list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

# Extract summary stats from each gene
all_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_only = m$P_aa_only,
      P_pcs_only = m$P_pcs_only,
      P_aa_with_pcs = m$P_aa_with_pcs
    )
  }))
}))
stopifnot(nrow(all_results) > 0)

# Calculate -log10(P) for both analyses
all_results[, neg_log10_P_aa_only := -log10(P_aa_only)]
all_results[, neg_log10_P_aa_with_pcs := -log10(P_aa_with_pcs)]
all_results[is.infinite(neg_log10_P_aa_only), neg_log10_P_aa_only := NA]
all_results[is.infinite(neg_log10_P_aa_with_pcs), neg_log10_P_aa_with_pcs := NA]
sum(is.na(all_results))
plot(all_results$neg_log10_P_aa_with_pcs,all_results$neg_log10_P_aa_only)

# Sort and assign cumulative positions
all_results <- all_results[order(Gene, Position)]
gene_lengths <- all_results[, .(len = max(Position, na.rm = TRUE)), by = Gene]
gene_offsets <- c(0, cumsum(head(gene_lengths$len, -1)))
names(gene_offsets) <- gene_lengths$Gene
all_results[, pos_global := Position + gene_offsets[Gene]]

# Bonferroni threshold
bf <- 0.05 / nrow(all_results)
thresh <- quantile(all_results$P_aa_with_pcs,0.05)
all_results$gwas_hit <- all_results$P_aa_with_pcs < thresh
# Alternate colors for genes
gene_colors <- setNames(rep(c("steelblue3", "grey60"),
                            length.out = length(unique(all_results$Gene))),
                        unique(all_results$Gene))

# Gene labels
midpoints <- gene_lengths[, .(mid = gene_offsets[Gene] + len / 2), by = Gene]

# ----- Manhattan: AA only (no pop structure control) -----
ylim_max <- 1 + max(all_results$neg_log10_P_aa_only, na.rm = TRUE)

plot(all_results$pos_global, all_results$neg_log10_P_aa_only,
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "-log10(p)",
     ylim = c(0, ylim_max),
     main = "Residue Variation (No Pop Structure Control)",
     xaxt = "n")

for (g in unique(all_results$Gene)) {
  sub <- all_results[Gene == g]
  points(sub$pos_global, sub$neg_log10_P_aa_only, col = gene_colors[g], lwd = 1.5)
  sig <- sub[neg_log10_P_aa_only > -log10(bf)]
  if (nrow(sig) > 0) {
    points(sig$pos_global, sig$neg_log10_P_aa_only, col = gene_colors[g], pch = 16)
  }
}

abline(h = -log10(bf), col = "black", lty = 2)
axis(1, at = midpoints$mid, labels = midpoints$Gene, las = 2, cex.axis = 0.7)

# ----- Manhattan: AA with PCs (pop structure controlled) -----
ylim_max <- 1 + max(all_results$neg_log10_P_aa_with_pcs, na.rm = TRUE)

plot(all_results$pos_global, all_results$neg_log10_P_aa_with_pcs,
     type = "n",
     xlab = "Aligned position across genes",
     ylab = "-log10(p)",
     ylim = c(0, ylim_max),
     main = "Residue Variation (With Pop Structure Control)",
     xaxt = "n")

for (g in unique(all_results$Gene)) {
  sub <- all_results[Gene == g]
  points(sub$pos_global, sub$neg_log10_P_aa_with_pcs, col = gene_colors[g], lwd = 1.5)
  sig <- sub[neg_log10_P_aa_with_pcs > -log10(bf)]
  if (nrow(sig) > 0) {
    points(sig$pos_global, sig$neg_log10_P_aa_with_pcs, col = gene_colors[g], pch = 16)
  }
}

abline(h = -log10(bf), col = "black", lty = 2)
axis(1, at = midpoints$mid, labels = midpoints$Gene, las = 2, cex.axis = 0.7)

# ----- Scatter plot: AA only vs AA with PCs -----
plot(all_results$neg_log10_P_aa_only, all_results$neg_log10_P_aa_with_pcs,
     xlab = "-log10(p) AA only",
     ylab = "-log10(p) AA with PCs",
     main = "No Control vs Pop Structure Control",
     pch = 16, col = rgb(0, 0, 0, 0.3))
cor(all_results$P_aa_only, all_results$P_aa_with_pcs)
     
#idx<- all_results$P_aa_with_pcs < quantile(all_results$P_aa_with_pcs,0.05)
#points(all_results$neg_log10_P_aa_only[idx], all_results$neg_log10_P_aa_with_pcs[idx] ,pch=1,col="red")

abline(0, 1, col = "red", lty = 2)
abline(h = -log10(bf), col = "grey50", lty = 3)
abline(v = -log10(bf), col = "grey50", lty = 3)

# ---- allele freq control ----
library(data.table)
library(dplyr)
library(tidyr)

model_files <- list.files("results/residue_models_triple/", pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

# Define amino acid columns
aa_cols <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
cols <- c("Gene", "P_aa_only", "P_aa_with_pcs", aa_cols)

# Extract data with amino acid counts
all_res <- list()
for(f in model_files){
  models <- readRDS(f)
  res <- data.frame(matrix(ncol=length(cols), nrow=length(models)))
  names(res) <- cols
  print(f)
  
  for(i in seq_along(models)){
    mod <- models[[i]]
    
    counts_vec <- unlist(mod$residue_counts)
    counts_vec <- counts_vec[names(counts_vec) %in% aa_cols]
    
    res[i, aa_cols] <- 0
    res[i, names(counts_vec)] <- counts_vec
    res$P_aa_only[i] <- mod$P_aa_only
    res$P_aa_with_pcs[i] <- mod$P_aa_with_pcs
    res$Gene[i] <- mod$Gene
  }
  
  all_res[[f]] <- res
}

# Combine all results
final_res <- do.call(rbind, all_res)
stopifnot(nrow(final_res) > 0)
stopifnot("P_aa_only" %in% names(final_res))
stopifnot("P_aa_with_pcs" %in% names(final_res))

# Calculate allele frequency for each site
final_res$site <- seq_len(nrow(final_res))

# Get total count per site (excluding gaps)
aa_cols_no_gap <- setdiff(aa_cols, "-")
final_res$total_alleles <- rowSums(final_res[, aa_cols_no_gap], na.rm = TRUE)

     # Calculate minor allele frequency (MAF) - frequency of second most common allele
final_res$MAF <- apply(final_res[, aa_cols_no_gap], 1, function(x) {
  counts <- sort(x[x > 0], decreasing = TRUE)
  if(length(counts) < 2) return(0)
  if(counts[2] < 10) return(0)
  counts[2] / sum(counts)
})
hist(final_res$MAF)
# Calculate number of alleles (polymorphism level)
final_res$n_alleles <- apply(final_res[, aa_cols_no_gap], 1, function(x) {
  sum(x > 0)
})
hist(final_res$n_alleles)
LSD::heatscatter(final_res$n_alleles, final_res$MAF)
# Bin by MAF
final_res$MAF_bin <- cut(final_res$MAF, 
                         breaks = c(0, 0.01, 0.05, 0.1, 0.2, 1),
                         labels = c("0.001-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.1"),
                         include.lowest = TRUE)
table(final_res$MAF_bin)
# Calculate -log10(p) difference (effect of pop structure control)
final_res$log10p_aa_only <- -log10(final_res$P_aa_only)
final_res$log10p_aa_with_pcs <- -log10(final_res$P_aa_with_pcs)
final_res$pop_control_effect <- final_res$log10p_aa_only - final_res$log10p_aa_with_pcs

# Remove infinite values
final_res <- final_res[is.finite(final_res$pop_control_effect), ]
stopifnot(nrow(final_res) > 0)

# ----- Analysis by MAF bin -----
library(ggplot2)

# Boxplot: pop structure control effect by MAF
ggplot(final_res, aes(x = MAF_bin, y = pop_control_effect)) +
  geom_boxplot(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    title = "Population Structure Control Effect by Allele Frequency",
    x = "Minor Allele Frequency",
    y = "Δ-log10(p) (AA only - AA with PCs)",
    subtitle = "Positive values = signal lost after controlling for pop structure"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter: MAF vs pop control effect
ggplot(final_res, aes(x = entropy, y = pop_control_effect)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Population Structure Effect vs Minor Allele Frequency",
    x = "Entropy",
    y = "Δ-log10(p) (AA only - AA with PCs)"
  )

# ----- Statistical test: does MAF predict pop control effect? -----
lm_maf <- lm(pop_control_effect ~ MAF + n_alleles, data = final_res)
print(summary(lm_maf))

# ----- Classify sites by significance -----
bf_threshold <- 0.05 / nrow(final_res)
thresh <- quantile(final_res$P_aa_with_pcs, 0.05)
final_res$sig_class <- case_when(
  final_res$P_aa_only < thresh & final_res$P_aa_with_pcs < thresh ~ "sig_both",
  final_res$P_aa_only < thresh & final_res$P_aa_with_pcs >= thresh ~ "sig_lost",
  final_res$P_aa_only >= thresh & final_res$P_aa_with_pcs < thresh ~ "sig_gained",
  TRUE ~ "not_sig"
)

table(final_res$sig_class)

# Compare MAF distribution across significance classes
ggplot(final_res, aes(x = sig_class, y = MAF, fill = sig_class)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Allele Frequency by Significance Class",
    x = "Significance Class",
    y = "Minor Allele Frequency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ----- Summary statistics -----
summary_stats <- final_res %>%
  group_by(sig_class) %>%
  summarise(
    n = n(),
    mean_MAF = mean(MAF, na.rm = TRUE),
    median_MAF = median(MAF, na.rm = TRUE),
    mean_pop_effect = mean(pop_control_effect, na.rm = TRUE),
    median_pop_effect = median(pop_control_effect, na.rm = TRUE)
  )

print(summary_stats)

# ---- control with allele_variance
final_res$allele_variance <- apply(final_res[, aa_cols_no_gap], 1, function(x) {
  p <- x / sum(x)
  var(p[p > 0])
})

sum(is.na(final_res$allele_variance))

hist(final_res$allele_variance)
plot(final_res$MAF,final_res$allele_variance)
plot(final_res$n_alleles,final_res$allele_variance)
cor(final_res$MAF,final_res$allele_variance)
cor(na.omit(final_res$n_alleles),na.omit(final_res$allele_variance))

final_res$entropy <- apply(final_res[, aa_cols_no_gap], 1, function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log(p))
})

hist(final_res$entropy)
plot(final_res$MAF,final_res$entropy)
plot(final_res$n_alleles,final_res$entropy)
plot(final_res$entropy,final_res$pop_control_effect)
cor(final_res$entropy,final_res$pop_control_effect, method="spearman")

lm_entropy <- lm(pop_control_effect ~ entropy + n_alleles, data = final_res)
print(summary(lm_entropy))


# ---- hydrophobicity again ---- 
library(arrow)
library(data.table)
library(Peptides)
library(ggplot2)

# Load AA properties
data(AAdata)
aa_hydrophobicity <- AAdata$Hydrophobicity$Tanford
aa_hydrophobicity <- c(A=3.14, R=3.14, D=0, N=0, C=4.19, E=0, 
                       Q=0, G=0, H=0, I=12.35, L=10.05, K=6.28, 
                       M=5.44, F=11.1, P=10.89, S=0, T=1.88, W=12.56, 
                       Y=11.93, V=7.12)
# Read all model files
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

# Define amino acid columns (exclude gaps and unknowns)
aa_cols_no_gap <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

aa_mass <- c(A=71.08, R=156.2, D=115.09, N=114.11, C=103.14, E=129.12, 
             Q=128.14, G=57.06, H=137.15, I=113.17, L=113.17, K=128.18, 
             M=131.21, F=147.18, P=97.12, S=87.08, T=101.11, W=186.21, 
             Y=163.18, V=99.14)

aa_volume <- c(A=88.6, R=173.4, D=111.1, N=117.7, C=108.5, E=138.4, 
               Q=143.9, G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, 
               M=162.9, F=189.9, P=122.7, S=89, T=116.1, W=227.8, 
               Y=193.6, V=140)

aa_surface_area <- c(A=115, R=225, D=150, N=160, C=135, E=190, 
                     Q=180, G=75, H=195, I=175, L=170, K=200, 
                     M=185, F=210, P=145, S=115, T=140, W=255, 
                     Y=230, V=155)
aa_df <- data.frame(
  aa = aa_cols_no_gap,
  mass = aa_mass[aa_cols_no_gap],
  volume = aa_volume[aa_cols_no_gap],
  surface_area = aa_surface_area[aa_cols_no_gap],
  hydrophobicity = aa_hydrophobicity[aa_cols_no_gap],
  stringsAsFactors = FALSE
)
rownames(aa_df) <- NULL
aa_pca <- prcomp(aa_df[,-1], scale. = T)
aa_pca
aa_pca$sdev^2/sum(aa_pca$sdev^2)
barplot(aa_pca$sdev^2/sum(aa_pca$sdev^2))
plot(aa_pca$x[,1],aa_pca$x[,2],
     main="AA PCA on mass, vol, surface area, and hydrophobicity",
     col="white",
     xlab="PC1, 81%, SizeVolMass, little hydro",
     ylab="PC2, 18%, Hydrophobicity, mass vol")
text(aa_pca$x[,1],aa_pca$x[,2],aa_df$aa)

# Extract site-level data with amino acid counts
all_sites <- list()
for (f in model_files) {
  models <- readRDS(f)
  
  for (i in seq_along(models)) {
    m <- models[[i]]
    counts <- unlist(m$residue_counts)
    counts <- counts[names(counts) %in% aa_cols_no_gap]
    
    if (length(counts) == 0) next
    
    # Filter to valid AAs with property data and count > 10
    valid_aa <- names(counts)[names(counts) %in% names(aa_hydrophobicity) & 
                                names(counts) %in% names(aa_mass) &
                                names(counts) %in% names(aa_volume) &
                                names(counts) %in% names(aa_surface_area)]
    if (length(valid_aa) == 0) next
    
    counts_valid <- counts[valid_aa]
    counts_valid <- counts_valid[counts_valid > 10]
    
    if (length(counts_valid) == 0) next
    
    # Get property values for valid AAs
    hydro_vals <- aa_hydrophobicity[names(counts_valid)]
    mass_vals <- aa_mass[names(counts_valid)]
    volume_vals <- aa_volume[names(counts_valid)]
    surface_vals <- aa_surface_area[names(counts_valid)]
    
    # Mean properties (all alleles)
    mean_hydro <- sum(counts_valid * hydro_vals) / sum(counts_valid)
    mean_mass <- sum(counts_valid * mass_vals) / sum(counts_valid)
    mean_volume <- sum(counts_valid * volume_vals) / sum(counts_valid)
    mean_surface <- sum(counts_valid * surface_vals) / sum(counts_valid)
    
    # Mean properties (minor alleles only, count < 4000)
    counts_minor <- counts_valid[counts_valid < 4000]
    
    if (length(counts_minor) > 0) {
      mean_minor_hydro <- sum(counts_minor * aa_hydrophobicity[names(counts_minor)]) / sum(counts_minor)
      mean_minor_mass <- sum(counts_minor * aa_mass[names(counts_minor)]) / sum(counts_minor)
      mean_minor_volume <- sum(counts_minor * aa_volume[names(counts_minor)]) / sum(counts_minor)
      mean_minor_surface <- sum(counts_minor * aa_surface_area[names(counts_minor)]) / sum(counts_minor)
    } else {
      mean_minor_hydro <- mean_minor_mass <- mean_minor_volume <- mean_minor_surface <- NA
    }
    
    # Mean properties (major alleles only, count >= 4000)
    counts_major <- counts_valid[counts_valid >= 4000]
    
    if (length(counts_major) > 0) {
      mean_major_hydro <- sum(counts_major * aa_hydrophobicity[names(counts_major)]) / sum(counts_major)
      mean_major_mass <- sum(counts_major * aa_mass[names(counts_major)]) / sum(counts_major)
      mean_major_volume <- sum(counts_major * aa_volume[names(counts_major)]) / sum(counts_major)
      mean_major_surface <- sum(counts_major * aa_surface_area[names(counts_major)]) / sum(counts_major)
    } else {
      mean_major_hydro <- mean_major_mass <- mean_major_volume <- mean_major_surface <- NA
    }
    
    all_sites[[length(all_sites) + 1]] <- data.frame(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_only = m$P_aa_only,
      P_aa_with_pcs = m$P_aa_with_pcs,
      mean_hydrophobicity = mean_hydro,
      mean_mass = mean_mass,
      mean_volume = mean_volume,
      mean_surface_area = mean_surface,
      mean_minor_hydro = mean_minor_hydro,
      mean_minor_mass = mean_minor_mass,
      mean_minor_volume = mean_minor_volume,
      mean_minor_surface = mean_minor_surface,
      mean_major_hydro = mean_major_hydro,
      mean_major_mass = mean_major_mass,
      mean_major_volume = mean_major_volume,
      mean_major_surface = mean_major_surface,
      n_alleles = length(counts_valid),
      n_minor_alleles = length(counts_minor),
      total_count = sum(counts_valid),
      stringsAsFactors = FALSE
    )
  }
}

sites_df <- do.call(rbind, all_sites)
stopifnot(nrow(sites_df) > 0)

# Calculate -log10(p) values
sites_df$log10p_aa_only <- -log10(sites_df$P_aa_only)
sites_df$log10p_aa_with_pcs <- -log10(sites_df$P_aa_with_pcs)

# Remove infinite values
sites_df <- sites_df[is.finite(sites_df$log10p_aa_only) & 
                       is.finite(sites_df$log10p_aa_with_pcs), ]
stopifnot(nrow(sites_df) > 0)

# Classify significance
bf_threshold <- 0.05 / nrow(sites_df)
threshold <- quantile(sites_df$P_aa_with_pcs,0.05)
sites_df$sig_class <- case_when(
  sites_df$P_aa_only < threshold & sites_df$P_aa_with_pcs < threshold ~ "sig_both",
  sites_df$P_aa_only < threshold & sites_df$P_aa_with_pcs >= threshold ~ "sig_lost",
  sites_df$P_aa_only >= threshold & sites_df$P_aa_with_pcs < threshold ~ "sig_gained",
  TRUE ~ "not_sig"
)
hist(sites_df$mean_hydrophobicity[sites_df$P_aa_with_pcs<quantile(sites_df$P_aa_with_pcs,0.05)])
hist(sites_df$mean_hydrophobicity[sites_df$P_aa_with_pcs>quantile(sites_df$P_aa_with_pcs,0.05)])

gwas_thresh <- quantile(sites_df$P_aa_with_pcs,0.05)
sites_df$gwas_class <- case_when(
  sites_df$P_aa_with_pcs < gwas_thresh ~ "sig_aa_with_ctrl",
  sites_df$P_aa_with_pcs >= gwas_thresh ~ "insig_aa_with_ctrl",
)
table(sites_df$gwas_class)
# ----- Hydrophobicity analysis -----
boxplot(mean_hydrophobicity ~ gwas_class, data=sites_df,
        main="Site mean hydrophobicity by gwas significance")
boxplot(mean_minor_hydro ~ gwas_class, data=sites_df)
boxplot(mean_major_hydro ~ gwas_class, data=sites_df)

par(mfrow=c(2,2))
boxplot(mean_hydrophobicity ~ gwas_class, data=sites_df,
        main="Site mean hydrophobicity by gwas significance",
        ylab="mean hydrophobicity (tanford 1962, kjMol-1)")
boxplot(mean_mass ~ gwas_class, data=sites_df,
        main="Site mean mass",
        ylab="mean mass (Da)")
boxplot(mean_volume ~ gwas_class, data=sites_df,
        main="Site mean Volume",
        ylab="mean Volume (Å^3)")
boxplot(mean_surface_area ~ gwas_class, data=sites_df,
        main="Site mean accessible surface area",
        ylab="mean surface area (Å^2)")

# ---- AA properties ananlysis ----
gwas_thresh <- quantile(sites_df$P_aa_with_pcs, 0.05)
sites_df$gwas_class <- case_when(
  sites_df$P_aa_with_pcs < gwas_thresh ~ "sig_aa_with_ctrl",
  sites_df$P_aa_with_pcs >= gwas_thresh ~ "insig_aa_with_ctrl"
)
table(sites_df$gwas_class)

# Statistical tests
test_hydro <- t.test(mean_hydro ~ gwas_class, data = sites_df)
test_mass <- t.test(mean_mass ~ gwas_class, data = sites_df)
test_volume <- t.test(mean_volume ~ gwas_class, data = sites_df)
test_surface <- t.test(mean_surface_area ~ gwas_class, data = sites_df)

# Format p-values
format_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

# Plotting
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

boxplot(mean_hydrophobicity ~ gwas_class, data = sites_df,
        col = c("lightblue", "lightgray"),
        main = "Mean Hydrophobicity",
        ylab = "Hydrophobicity (Tanford)",
        xlab = "",
        names = c("Significant", "Non-significant"),
        las = 1)
text(x = 1.5, y = max(sites_df$mean_hydrophobicity, na.rm = TRUE) * 0.95,
     labels = format_p(test_hydro$p.value), cex = 0.9)

boxplot(mean_mass ~ gwas_class, data = sites_df,
        col = c("lightblue", "lightgray"),
        main = "Mean Mass",
        ylab = "Mass (Da)",
        xlab = "",
        names = c("Significant", "Non-significant"),
        las = 1)
text(x = 1.5, y = max(sites_df$mean_mass, na.rm = TRUE) * 0.95,
     labels = format_p(test_mass$p.value), cex = 0.9)

boxplot(mean_volume ~ gwas_class, data = sites_df,
        col = c("lightblue", "lightgray"),
        main = "Mean Volume",
        ylab = expression(paste("Volume (", ring(A)^3, ")")),
        xlab = "",
        names = c("Significant", "Non-significant"),
        las = 1)
text(x = 1.5, y = max(sites_df$mean_volume, na.rm = TRUE) * 0.95,
     labels = format_p(test_volume$p.value), cex = 0.9)

boxplot(mean_surface_area ~ gwas_class, data = sites_df,
        col = c("lightblue", "lightgray"),
        main = "Mean Accessible Surface Area",
        ylab = expression(paste("Surface Area (", ring(A)^2, ")")),
        xlab = "",
        names = c("Significant", "Non-significant"),
        las = 1)
text(x = 1.5, y = max(sites_df$mean_surface_area, na.rm = TRUE) * 0.95,
     labels = format_p(test_surface$p.value), cex = 0.9)


sites_df$complex <- substr(sites_df$Gene, )

# Compare mean hydrophobicity across significance classes
ggplot(sites_df, aes(x = gwas_class, y = mean_hydrophobicity, fill = gwas_class)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw() +
  labs(
    title = "Mean Site Hydrophobicity by GWAS Class",
    x = "GWAS Class",
    y = "Mean Hydrophobicity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(sites_df, aes(x = sig_class, y = mean_major_hydro, fill = sig_class)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw() +
  labs(
    title = "Mean Site Minor Hydrophobicity by Sig Class",
    x = "GWAS Class",
    y = "Mean Hydrophobicity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter: hydrophobicity vs significance
ggplot(sites_df, aes(x = mean_hydrophobicity, y = log10p_aa_with_pcs)) +
  geom_point(alpha = 0.8, aes(color = sig_class)) +
  geom_hline(yintercept = -log10(bf_threshold), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    title = "Site Hydrophobicity vs Significance (Pop Structure Controlled)",
    x = "Mean Hydrophobicity",
    y = "-log10(p) with pop structure control"
  )

# Scatter: hydrophobicity vs pop structure effect
sites_df$pop_effect <- sites_df$log10p_aa_only - sites_df$log10p_aa_with_pcs

ggplot(sites_df, aes(x = mean_hydrophobicity, y = pop_effect)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Population Structure Effect vs Site Hydrophobicity",
    x = "Mean Hydrophobicity",
    y = "Δ-log10(p) (AA only - AA with PCs)"
  )

# Summary statistics
summary_hydro <- sites_df %>%
  group_by(sig_class) %>%
  summarise(
    n = n(),
    mean_hydro = mean(mean_hydrophobicity),
    sd_hydro = sd(mean_hydrophobicity),
    median_hydro = median(mean_hydrophobicity)
  )

print(summary_hydro)

# Statistical test: does hydrophobicity predict significance?
lm_hydro <- lm(log10p_aa_with_pcs ~ mean_hydrophobicity + n_alleles, data = sites_df)
print(summary(lm_hydro))

# Test if hydrophobicity differs between sig classes
kruskal.test(mean_hydrophobicity ~ sig_class, data = sites_df)

# ---- IQTREE prep hit vs not ----
library(data.table)
library(Biostrings)

# 1. Load GWAS site annotations (which sites are in gwas class?)

model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)

all_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene     = m$Gene,
      Position = m$Aligned_Position,
      p = m$P_aa_with_pcs
    )
  }))
}))
all_results$gwas_hit <- all_results$p < quantile(all_results$p, 0.05)

# Just keep Gene, Position, gwas_hit
setkey(all_results, Gene, Position)

# 2. Build the superalignment for all species + filtered genes

aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", 
                        full.names = TRUE)

get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes <- get_gene(aln_files)

# List of all species across all genes (intersection)
species_sets <- list()

aligned_gene_list <- list()

for (file in aln_files) {
  gene <- get_gene(file)
  aln <- readAAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))   # strip extra ID metadata
  
  aligned_gene_list[[gene]] <- aln
  species_sets[[gene]] <- names(aln)
}

# Universal species set: must appear in *all* alignments
all_species <- Reduce(union, species_sets)

super_list <- list()
partition_info <- list()

global_offset <- 0

MIN_NON_GAP <- 8500   # threshold for phylogenetic informativeness

for (gene in genes) {
  aln <- aligned_gene_list[[gene]]
  aln_mat <- as.matrix(aln)
  gene_len <- ncol(aln_mat)
  
  # Place into full species matrix (fill missing species with all gaps)
  mat <- matrix("-", nrow = length(all_species), ncol = gene_len,
                dimnames = list(all_species, NULL))
  present <- match(names(aln), all_species)
  mat[present[!is.na(present)], ] <- aln_mat
  
  # QC per column
  keep <- sapply(seq_len(gene_len), function(j) {
    col <- mat[, j]
    non_gap <- col[col != "-"]
    if (length(non_gap) < MIN_NON_GAP) return(FALSE)
    if (length(unique(non_gap)) < 2) return(FALSE)
    TRUE
  })
  
  kept_mat <- mat[, keep, drop = FALSE]
  n_kept <- ncol(kept_mat)
  
  if (n_kept > 0) {
    # Store surviving sites
    super_list[[gene]] <- kept_mat
    
    # Map kept positions to global coordinates
    kept_positions <- which(keep)
    partition_info[[gene]] <- data.table(
      Gene = gene,
      GenePos = kept_positions,
      GlobalPos = global_offset + seq_len(n_kept)
    )
    
    global_offset <- global_offset + n_kept
  } else {
    cat("Gene", gene, "dropped - no informative columns remaining.\n")
  }
}

cat("Total superalignment length:", global_offset, "sites\n")

# Concatenate horizontally
supermat <- do.call(cbind, super_list)


# 4. Write SUPERGENE alignment in FASTA format
out_fasta <- "iqtree_input/supergene.faa"
dir.create("iqtree_input", showWarnings = FALSE, recursive = TRUE)

aa_strings <- apply(supermat, 1, paste0, collapse = "")

write.fasta <- function(seq, ids, file) {
  con <- file(file, "w")
  for (i in seq_along(seq)) {
    writeLines(paste0(">", ids[i]), con)
    writeLines(seq[[i]], con)
  }
  close(con)
}

write.fasta(as.list(aa_strings), all_species, out_fasta)

cat("Wrote supergene alignment:", out_fasta, "\n")


# 5. Build IQ-TREE partition file labeling GWAS vs Other
part_DT <- rbindlist(partition_info)

# merge with GWAS flags
part_DT <- merge(
  part_DT, 
  all_results, 
  by.x = c("Gene", "GenePos"), 
  by.y = c("Gene", "Position"),
  all.x = TRUE
)

part_DT[is.na(gwas_hit), gwas_hit := FALSE]

# contiguous blocks
partition_file <- "iqtree_input/partition_filtered.nex"

gwas_sites  <- part_DT[gwas_hit == TRUE,  GlobalPos]
other_sites <- part_DT[gwas_hit == FALSE, GlobalPos]

cat("#NEXUS\nBEGIN SETS;\n", file = partition_file)

# GWAS partition
if (length(gwas_sites) > 0) {
  cat("CHARSET GWAS = ", file = partition_file, append = TRUE)
  cat(paste(gwas_sites, collapse = " "), file = partition_file, append = TRUE)
  cat(";\n", file = partition_file, append = TRUE)
}

# Other partition
if (length(other_sites) > 0) {
  cat("CHARSET Other = ", file = partition_file, append = TRUE)
  cat(paste(other_sites, collapse = " "), file = partition_file, append = TRUE)
  cat(";\n", file = partition_file, append = TRUE)
}

cat("END;\n", file = partition_file, append = TRUE)

cat("Wrote:", partition_file, "\n")

cmd1 <- "
/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 \
  -s supergene.faa \
  -q partition_filtered.nex \
  -m LG \
  -nt AUTO \
  -pre null_one_matrix
"

cmd2 <- "
iqtree2 \
  -s supergene.faa \
  -q partition_filtered.nex \
  -m Q.p1+Q.p2 \
  -mrate EQUAL \
  -nt AUTO \
  -pre alt_two_matrices
"

# this is slwo, lets try
# /programs/FastTree-2.1.11/FastTree -lg -gamma supergene.faa > fast.tree
# /programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -s supergene.faa -q partition_filtered.nex -m LG+G -te fast.tree -nt 88 --safe -asr -pre iq_asr

str(sites_df)

# ---- raxml prep ----
library(data.table)
library(Biostrings)
library(ape)
library(arrow)

# LOAD COMMON DATA
data <- as.data.table(read_parquet("data/processed_data.parquet"))

# Load embedding QC data
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))

# Create clean ID sets per gene (exclude ManualOutlier == TRUE)
setDT(embeds_with_mds)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

# Get union of all clean IDs across genes
all_clean_ids <- unique(clean_ids_by_gene$ID)
cat("Total clean species after outlier removal:", length(all_clean_ids), "\n")

# PART 1: BUILD CDS SUPERALIGNMENT

aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_CDS_aligned\\.fasta$", 
                        full.names = TRUE)
get_gene <- function(path) sub("_CDS_aligned\\.fasta", "", basename(path))
genes <- get_gene(aln_files)

species_sets <- list()
aligned_gene_list <- list()

for (file in aln_files) {
  gene <- get_gene(file)
  aln <- readDNAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  
  # Filter to clean IDs for this gene
  gene_clean_ids <- clean_ids_by_gene[Gene == gene, ID]
  aln <- aln[names(aln) %in% gene_clean_ids]
  
  if (length(aln) > 0) {
    aligned_gene_list[[gene]] <- aln
    species_sets[[gene]] <- names(aln)
  } else {
    cat("Gene", gene, "dropped - no clean sequences.\n")
  }
}

all_species <- Reduce(union, species_sets)
cat("Total species after gene union:", length(all_species), "\n")

super_list <- list()
global_offset <- 0
MIN_NON_GAP <- 8500

for (gene in names(aligned_gene_list)) {
  aln <- aligned_gene_list[[gene]]
  aln_mat <- as.matrix(aln)
  gene_len <- ncol(aln_mat)
  
  # Build full species matrix
  mat <- matrix("-", nrow = length(all_species), ncol = gene_len,
                dimnames = list(all_species, NULL))
  present <- match(names(aln), all_species)
  mat[present[!is.na(present)], ] <- aln_mat
  
  # QC per column
  keep <- sapply(seq_len(gene_len), function(j) {
    col <- mat[, j]
    non_gap <- col[col != "-"]
    if (length(non_gap) < MIN_NON_GAP) return(FALSE)
    if (length(unique(non_gap)) < 2) return(FALSE)
    TRUE
  })
  
  kept_mat <- mat[, keep, drop = FALSE]
  n_kept <- ncol(kept_mat)
  
  if (n_kept > 0) {
    super_list[[gene]] <- kept_mat
    global_offset <- global_offset + n_kept
  } else {
    cat("Gene", gene, "dropped - no informative sites.\n")
  }
}
cat("Total CDS superalignment length:", global_offset, "bp\n")

supermat <- do.call(cbind, super_list)

# Write FASTA
out_fasta <- "raxml_input/supercds.fasta"
dir.create("raxml_input", showWarnings = FALSE, recursive = TRUE)

dna_strings <- apply(supermat, 1, paste0, collapse = "")
write.fasta <- function(seq, ids, file) {
  con <- file(file, "w")
  for (i in seq_along(seq)) {
    writeLines(paste0(">", ids[i]), con)
    writeLines(seq[[i]], con)
  }
  close(con)
}
write.fasta(as.list(dna_strings), all_species, out_fasta)
cat("Wrote CDS superalignment:", out_fasta, "\n")


# PART 2: DISTANCE-INFORMED GENUS-LEVEL COLLAPSING

cat("\nCollapsing closely related congeneric species...\n")

# Map IDs to organism names
id_to_organism <- data[, .(ID, Organism)]
setkey(id_to_organism, ID)

species_info <- data.table(
  ID = all_species,
  Organism = id_to_organism[all_species, Organism]
)
species_info[, genus := sub("^([^ ]+) .*$", "\\1", Organism)]

cat("Total species:", length(all_species), "\n")
cat("Total genera:", uniqueN(species_info$genus), "\n")

# Sample subset for distance estimation (max 2000 species for speed)
set.seed(123)
SAMPLE_SIZE <- min(2000, length(all_species))
sampled_ids <- sample(all_species, SAMPLE_SIZE)

cat("Sampling", SAMPLE_SIZE, "species for distance estimation...\n")
aln_dna <- readDNAStringSet(out_fasta)
sampled_aln <- aln_dna[sampled_ids]
dist_mat_sample <- dist.dna(as.DNAbin(sampled_aln), model = "raw", pairwise.deletion = TRUE)

# Get distribution of distances
dist_vals <- as.vector(dist_mat_sample)
dist_vals <- dist_vals[!is.na(dist_vals)]
DISTANCE_THRESHOLD <- quantile(dist_vals, 0.001)  # 5th percentile = very close

cat("Distance threshold (5th percentile):", round(DISTANCE_THRESHOLD, 6), "\n")
hist(dist_vals)
abline(v=DISTANCE_THRESHOLD)
as.vector(DISTANCE_THRESHOLD) * length(aln_dna$PX136589.1)
# For each genus with multiple species, collapse based on distance
keep_species <- character()

for (g in unique(species_info$genus)) {
  genus_species <- species_info[genus == g, ID]
  
  if (length(genus_species) == 1) {
    # Only one species in genus, keep it
    keep_species <- c(keep_species, genus_species)
  } else {
    # Multiple species - compute pairwise distances within genus
    genus_aln <- aln_dna[genus_species]
    
    tryCatch({
      genus_dist <- dist.dna(as.DNAbin(genus_aln), model = "raw", pairwise.deletion = TRUE)
      genus_dist_mat <- as.matrix(genus_dist)
      
      # Greedy clustering: iteratively pick representatives
      remaining <- genus_species
      representatives <- character()
      
      while (length(remaining) > 0) {
        # Pick first remaining species as representative
        rep <- remaining[1]
        representatives <- c(representatives, rep)
        
        # Remove species within distance threshold of this representative
        if (length(remaining) > 1) {
          rep_idx <- which(genus_species == rep)
          rem_idx <- match(remaining, genus_species)
          distances_to_rep <- genus_dist_mat[rep_idx, rem_idx]
          
          # Keep species that are far enough away
          far_enough <- remaining[distances_to_rep > DISTANCE_THRESHOLD | is.na(distances_to_rep)]
          remaining <- far_enough[far_enough != rep]
        } else {
          remaining <- character()
        }
      }
      
      keep_species <- c(keep_species, representatives)
      
      if (length(representatives) < length(genus_species)) {
        cat("Genus", g, ": collapsed", length(genus_species), "to", length(representatives), "species\n")
      }
      
    }, error = function(e) {
      # If distance calculation fails, just keep first species
      keep_species <<- c(keep_species, genus_species[1])
      cat("Genus", g, ": distance calculation failed, keeping 1 representative\n")
    })
  }
}

cat("\nCollapsed from", length(all_species), "to", length(keep_species), "species\n")

# Filter alignment
filtered_aln <- aln_dna[keep_species]
out_fasta_filtered <- "raxml_input/supercds_collapsed.fasta"
writeXStringSet(filtered_aln, out_fasta_filtered)
cat("Wrote collapsed alignment:", out_fasta_filtered, "\n")

# PART 3: RAxML COMMAND GENERATION

raxml_cmd <- sprintf(
  "/programs/raxml-ng_v1.2.0/raxml-ng  -T 20 -s %s -n supercds_tree -m GTRGAMMA -p 12345 -x 12345 -# 100 -f a",
  out_fasta_filtered
)

cat("\n========================================\n")
cat("RAxML command:\n")
cat(raxml_cmd, "\n")
cat("========================================\n")

writeLines(raxml_cmd, "raxml_input/run_raxml.sh")
cat("\nSaved command to: raxml_input/run_raxml.sh\n")
cat("Run with: bash raxml_input/run_raxml.sh\n")

# SUMMARY STATS

cat("\n=== PIPELINE SUMMARY ===\n")
cat("Original species (after outlier filter):", length(all_clean_ids), "\n")
cat("Species in superalignment:", length(all_species), "\n")
cat("Final species (after tip collapse):", length(keep_species), "\n")
cat("Alignment length:", global_offset, "bp\n")
cat("Number of genes:", length(names(aligned_gene_list)), "\n")
cat("Ready for RAxML!\n")


# ---- STRUCTURAL ASSESSMENT ----
library(data.table)
library(stringr)
library(Biostrings)

netsurf_df_1 <- read.csv("data/693ABEC1003AD9BE13062023.csv")
netsurf_df_2 <- read.csv("data/6939DA900036F1B58359D2BE.csv")
netsurf_df_3 <- read.csv("data/693C4D9200013958CF3F917B.csv")
netsurf_df_4 <- read.csv("data/693C6D5E00019BDB56427C57.csv")
netsurf_df <- rbind(netsurf_df_1, netsurf_df_2,netsurf_df_3,netsurf_df_4)

# Parse netsurf IDs
netsurf_df$ID <- sub("^>", "", str_split_i(netsurf_df$id, "_Gene", 1))
gene_part <- str_split_i(netsurf_df$id, "_Gene", 2)
netsurf_df$Gene <- str_split_i(gene_part, "_", 2)
table(netsurf_df$Gene)
# Basic checks
stopifnot(all(!is.na(netsurf_df$ID)))
stopifnot(all(!is.na(netsurf_df$Gene)))
stopifnot(nrow(netsurf_df) > 0)

# Load alignments to map netsurf position (n) to alignment position
aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", 
                        full.names = TRUE)
get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))

# Create position mapping for each gene
netsurf_mapped_list <- list()

for (aln_file in aln_files) {
  gene <- get_gene(aln_file)
  
  # Get netsurf data for this gene
  netsurf_gene <- netsurf_df[netsurf_df$Gene == gene, ]
  if (nrow(netsurf_gene) == 0) next
  
  # Read alignment
  aln <- readAAStringSet(aln_file)
  aln_ids <- sub("\\|.*", "", names(aln))
  names(aln) <- aln_ids
  
  # Get netsurf IDs present in alignment
  netsurf_ids <- unique(netsurf_gene$ID)
  common_ids <- intersect(netsurf_ids, aln_ids)
  
  if (length(common_ids) == 0) {
    message("Gene ", gene, ": no overlap between netsurf and alignment")
    next
  }
  
  message("Gene ", gene, ": ", length(common_ids), " sequences with netsurf data")
  
  # For each ID, map ungapped position to alignment position
  for (id in common_ids) {
    seq_str <- as.character(aln[[id]])
    chars <- strsplit(seq_str, "")[[1]]
    
    # Build mapping: ungapped position -> alignment position
    ungapped_pos <- 0
    pos_map <- integer()
    
    for (aln_pos in seq_along(chars)) {
      if (chars[aln_pos] != "-") {
        ungapped_pos <- ungapped_pos + 1
        pos_map[ungapped_pos] <- aln_pos
      }
    }
    
    # Map netsurf positions to alignment positions
    netsurf_id <- netsurf_gene[netsurf_gene$ID == id, ]
    netsurf_id$Position <- pos_map[netsurf_id$n]
    
    stopifnot(all(!is.na(netsurf_id$Position)))
    
    netsurf_mapped_list[[paste(gene, id, sep = "_")]] <- netsurf_id
  }
}

netsurf_mapped <- rbindlist(netsurf_mapped_list, fill = TRUE)

atpA <- netsurf_mapped[netsurf_mapped$Gene=="atpA",]
hist(atpA$p.q3_H.)
hist(atpA$p.q3_E.)
hist(atpA$p.q3_C.)
table(atpA$q3)
boxplot(p.q3_H. ~ Position, atpA)
boxplot(p.q3_E. ~ Position, atpA)
boxplot(p.q3_C. ~ Position, atpA)
rbcL <- netsurf_mapped[netsurf_mapped$Gene=="rbcL",]
hist(rbcL$p.q3_H.)
hist(rbcL$p.q3_E.)
hist(rbcL$p.q3_C.)
table(rbcL$q3)


# Calculate structural variation by alignment position
# For each Gene + Position, compute variation in q3 states
struct_var <- netsurf_mapped[, .(
  n_seqs = .N,
  n_H = sum(q3 == "H"),
  n_E = sum(q3 == "E"),
  n_C = sum(q3 == "C"),
  prop_H = mean(q3 == "H"),
  prop_E = mean(q3 == "E"),
  prop_C = mean(q3 == "C"),
  entropy = {
    props <- c(mean(q3 == "H"), mean(q3 == "E"), mean(q3 == "C"))
    props <- props[props > 0]
    -sum(props * log2(props))
  },
  consensus = names(which.max(table(q3))),
  consensus_freq = max(table(q3)) / .N
), by = .(Gene, Position)]

# Classify positions by structural conservation
# High conservation: consensus_freq >= 0.8
# Moderate: 0.5 <= consensus_freq < 0.8
# Variable: consensus_freq < 0.5
struct_var$conservation <- cut(struct_var$consensus_freq,
                               breaks = c(0, 0.5, 0.95, 1),
                               labels = c("variable", "moderate", "high"),
                               include.lowest = TRUE)
rbcL_struct_var <- struct_var[struct_var$Gene=="rbcL",]

#struct_var <- struct_var[struct_var$Gene!="petG",]
#struct_var <- struct_var[struct_var$Gene!="rbcL",]
hist(rbcL_struct_var$consensus_freq, main="rbcL Q3 consensus freq, n=1000")

length(unique(struct_var$Gene))
hist(struct_var$consensus_freq, breaks = 50,
     main="Frequency of consensus structural assignment in subsamples, \n n=60 Genes, 100 samples~gene")
summary(struct_var$consensus_freq)
table(struct_var$conservation)
# Merge with GWAS sites
merged <- merge(sites_df, struct_var, 
                by = c("Gene", "Position"), 
                all.x = TRUE)

# Summary statistics
message("\n=== Structural Conservation Summary ===")
print(table(struct_var$conservation))
message("\nMean entropy: ", round(mean(struct_var$entropy), 3))
message("Mean consensus frequency: ", round(mean(struct_var$consensus_freq), 3))

# Test: are GWAS sites enriched for certain structural states?
# Only use positions with high conservation (consensus_freq >= 0.8)
high_cons <- merged[!is.na(merged$consensus_freq) & 
                      merged$consensus_freq >= 0.8, ]
plot(-log10(merged$P_aa_with_pcs), merged$consensus_freq, main="GWAS significance against consensus structure frequency")
message("\n=== Testing ", nrow(high_cons), " high-confidence structural sites ===")

results <- data.frame()
gwas_classes <- unique(sites_df$gwas_class)
q3_states <- c("H", "E", "C")
boxplot(-log10(P_aa_with_pcs) ~ consensus, high_cons,
        subset = gwas_class == "sig_aa_with_ctrl")
boxplot(P_aa_with_pcs ~ consensus, high_cons,
        subset = gwas_class == "sig_aa_with_ctrl",
        main="p of significant sites with pop structure control by Q3 structure",
        xlab="Structure, (Coil, Helix, E--> Beta Strand)",
        ylab = "p anova(temp ~ pcs, temp ~ aa + pcs)"
        )


hist(high_cons$P_aa_only,breaks=100,main="P value distribution of temp ~ aa")
abline(v=
quantile(high_cons$P_aa_only,0.25),col="red"
)
hist(high_cons$P_aa_with_pcs)
quantile(high_cons$P_aa_only,0.05)
hist(-log10(high_cons$P_aa_only),breaks=100,main="P value distribution of temp ~ aa")
abline(v=
         -log10(quantile(high_cons$P_aa_only,0.25)),col="red"
)

table(high_cons$P_aa_only < quantile(high_cons$P_aa_only,0.25))

high_cons$site_class <- NA
high_cons$site_class[high_cons$P_aa_only < quantile(high_cons$P_aa_only,0.25)] <- "sig_no_control"
high_cons$site_class[high_cons$P_aa_with_pcs < quantile(high_cons$P_aa_with_pcs,0.05)] <- "sig_with_control"
high_cons$site_class[is.na(high_cons$site_class)] <- "not_sig"

table(high_cons$site_class)
sum(high_cons$P_aa_only < quantile(high_cons$P_aa_only,0.25) &
  high_cons$P_aa_with_pcs < quantile(high_cons$P_aa_with_pcs,0.05))


table(high_cons$site_class)
plot(
  -log10(high_cons$P_aa_only),
  -log10(high_cons$P_aa_with_pcs),
  ylim=c(0,33)
)
abline(5,-1/20,col="red")


library(ggplot2)
library(ggExtra)

df <- transform(
  high_cons,
  x = -log10(P_aa_only),
  y = -log10(P_aa_with_pcs),
  site_class = site_class
)

p <- ggplot(df, aes(x, y, color=site_class)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(x = "-log10(P aa only)", y = "-log10(P aa with PCs)") +
  theme_classic()
p
table(df$site_class)
ggMarginal(p, type = "histogram", bins = 40)

par(mfrow=c(1,1))
hist(high_cons$P_aa_only[high_cons$P_aa_only<0.00000000001],breaks=60, main="P-val, no control")
abline(v=quantile(high_cons$P_aa_only,0.25), col="red")

par(mfrow=c(2,2))
hist(high_cons$P_aa_only,breaks=60, main="P-val, no control")
abline(v=quantile(high_cons$P_aa_only,0.25), col="red")
hist(-log10(high_cons$P_aa_only),breaks=30, main="-log10(P), no control")
abline(v=-log10(quantile(high_cons$P_aa_only,0.25)), col="red")
LSD::heatscatter(1:nrow(high_cons),-log10(high_cons$P_aa_only), main="log10p", xlab="index")
abline(h=-log10(quantile(high_cons$P_aa_only,0.25)),col="red")


x <- -log10(sort(high_cons$P_aa_only))
exp <- -log10(ppoints(length(x)))
-log10(0.01/10900)
plot(exp, x,
     xlab = "Expected -log10(P)",
     ylab = "Observed -log10(P)",
     main="QQplot for p values, no control")
abline(0, 1, col = "blue")
abline(h=-log10(quantile(high_cons$P_aa_only,0.25)),col="red")

boxplot(P_aa_only ~ consensus + site_class, 
        data = high_cons,
        col = rep(c("lightgray", "gold", "darkred"), each = 3),
        las = 2,
        ylab = "-log10(P-value)",
        xlab = "",
        main = "GWAS significance by structure",
        names = rep(c("C", "E", "H"), 3),
        at = c(1:3, 5:7, 9:11))

legend("topright", 
       legend = c("not_sig", "sig_no_control", "sig_with_control"),
       fill = c("lightgray", "gold", "darkred"),
       bty = "n")

axis(1, at = c(2, 6, 10), 
     labels = c("not_sig", "sig_no_ctrl", "sig_with_ctrl"),
     line = 1, tick = FALSE)


boxplot(-log10(P_aa_with_pcs) ~ consensus + site_class, 
        data = high_cons,
        col = rep(c("lightgray", "salmon", "lightblue", "lightgreen"), each = 3),
        las = 2,
        ylab = "-log10(P-value)",
        xlab = "",
        main = "GWAS significance by structure and significance class",
        names = rep(c("C", "E", "H"), 4),
        at = c(1:3, 5:7, 9:11, 13:15))
legend("topright", 
       legend = c("not_sig", "sig_both", "sig_gained", "sig_lost"),
       fill = c("lightgray", "salmon", "lightblue", "lightgreen"),
       bty = "n")

structure_table <- table(high_cons$site_class, high_cons$consensus)
barplot(structure_table, 
        beside = TRUE,
        col = c("lightgray", "gold", "darkred"),
        legend.text = TRUE,
        args.legend = list(x = "top", bty = "n"),
        ylab = "Count",
        xlab = "Structure",
        main = "Structure distribution by site class")




structure_table
structure_prop <- prop.table(structure_table, margin = 1)
round(structure_prop,3)

# Statistical tests for enrichment
results <- data.frame()
site_classes <- c("sig_no_control", "sig_with_control")
q3_states <- c("H", "E", "C")

plot_data <- high_cons
for (site_class in site_classes) {
  test_sites <- plot_data[plot_data$site_class == site_class, ]
  background <- plot_data[plot_data$site_class == "not_sig", ]
  
  stopifnot(nrow(test_sites) > 0)
  stopifnot(nrow(background) > 0)
  
  for (state in q3_states) {
    obs_count <- sum(test_sites$consensus == state)
    obs_prop <- obs_count / nrow(test_sites)
    
    bg_count <- sum(background$consensus == state)
    bg_prop <- bg_count / nrow(background)
    
    cont_table <- matrix(c(obs_count, nrow(test_sites) - obs_count,
                           bg_count, nrow(background) - bg_count),
                         nrow = 2, byrow = TRUE)
    
    fisher_test <- fisher.test(cont_table)
    
    results <- rbind(results, data.frame(
      site_class = site_class,
      q3_state = state,
      n_sites = nrow(test_sites),
      obs_count = obs_count,
      obs_prop = obs_prop,
      bg_prop = bg_prop,
      fold_enrichment = obs_prop / bg_prop,
      p_value = fisher_test$p.value,
      odds_ratio = as.numeric(fisher_test$estimate)
    ))
  }
}

results$p_adj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$p_adj), ]
print(results)

# Create visualization
par(mfrow = c(1, 2))

# Plot 1: Fold enrichment
results_wide <- reshape(results[, c("site_class", "q3_state", "fold_enrichment")],
                        direction = "wide",
                        idvar = "site_class",
                        timevar = "q3_state")
rownames(results_wide) <- results_wide$site_class
results_wide$site_class <- NULL
names(results_wide)
str_split_i(colnames(results_wide),"\\.",2)
colnames(results_wide) <- str_split_i(colnames(results_wide),"\\.",2)
names(results_wide)

barplot(as.matrix(t(results_wide)),
        beside = TRUE,
        col = c("steelblue", "orange", "purple"),
        ylab = "Fold enrichment vs not_sig",
        main = "Structure enrichment in significant sites",
        legend.text = F,
        args.legend = list(x = "topright", bty = "n", title = "Structure"))
abline(h = 1, lty = 2, col = "red")

# Plot 2: -log10(p-value) with significance threshold
results$neg_log_p <- -log10(results$p_adj)

results_wide_p <- reshape(results[, c("site_class", "q3_state", "neg_log_p")],
                          direction = "wide",
                          idvar = "site_class",
                          timevar = "q3_state")
rownames(results_wide_p) <- results_wide_p$site_class
results_wide_p$site_class <- NULL
names(results_wide_p)
colnames(results_wide_p) <- str_split_i(colnames(results_wide_p),"\\.",2)


barplot(as.matrix(t(results_wide_p)),
        beside = TRUE,
        col = c("steelblue", "orange", "purple"),
        ylab = "-log10(adjusted p-value)",
        main = "Significance of enrichment",
        legend.text = T,
        args.legend = list(x = "topright", bty = "n", title = "Structure"))
abline(h = -log10(0.05), lty = 2, col = "red", lwd = 2)
text(x = 1, y = -log10(0.05) + 0.5, labels = "p = 0.05", pos = 3, col = "red")

par(mfrow = c(1, 1))


# Print results more clearly
cat("\n=== Enrichment Analysis vs not_sig background ===\n")
print(results[, c("site_class", "q3_state", "fold_enrichment", "p_adj")])

# Summary interpretation
cat("\n=== Summary ===\n")
for (site_class in c("sig_no_control", "sig_with_control")) {
  cat("\n", site_class, ":\n", sep = "")
  subset_res <- results[results$site_class == site_class, ]
  for (i in 1:nrow(subset_res)) {
    with(subset_res[i, ], {
      direction <- ifelse(fold_enrichment > 1, "ENRICHED", "DEPLETED")
      sig <- ifelse(p_adj < 0.05, "***", 
                    ifelse(p_adj < 0.10, "*", "ns"))
      cat(sprintf("  %s: %.2fx %s (p_adj=%.4f) %s\n", 
                  q3_state, fold_enrichment, direction, p_adj, sig))
    })
  }
}

# Quick visual check of proportions
cat("\n=== Proportion comparison ===\n")
prop_table <- prop.table(structure_table, margin = 1)
print(round(prop_table, 3))

cor(-log10(high_cons$P_aa_only), -log10(high_cons$P_aa_with_pcs))

boxplot(P_aa_with_pcs ~ consensus + sig_class, 
        data = high_cons,
        col = rep(c("lightgray", "salmon", "lightblue", "lightgreen"), each = 3),
        las = 2,
        ylab = "(P-value",
        xlab = "",
        main = "GWAS significance by structure and significance class",
        names = rep(c("C", "E", "H"), 4),
        at = c(1:3, 5:7, 9:11, 13:15))
legend("topright", 
       legend = c("not_sig", "sig_both", "sig_gained", "sig_lost"),
       fill = c("lightgray", "salmon", "lightblue", "lightgreen"),
       bty = "n")

boxplot(-log10(P_aa_only) ~ consensus + sig_class, 
        data = high_cons,
        col = rep(c("lightgray", "salmon", "lightblue", "lightgreen"), each = 3),
        las = 2,
        ylab = "(P-value",
        xlab = "",
        main = "GWAS significance by structure and significance class",
        names = rep(c("C", "E", "H"), 4),
        at = c(1:3, 5:7, 9:11, 13:15))
legend("topright", 
       legend = c("not_sig", "sig_both", "sig_gained", "sig_lost"),
       fill = c("lightgray", "salmon", "lightblue", "lightgreen"),
       bty = "n")



table(high_cons$consensus)
results <- data.frame()
gwas_classes <- unique(sites_df$gwas_class)
q3_states <- c("H", "E", "C")

for (gwas_class in gwas_classes) {
  gwas_sites <- high_cons[high_cons$gwas_class == gwas_class, ]
  background <- high_cons[high_cons$gwas_class != gwas_class, ]
  
  stopifnot(nrow(gwas_sites) > 0)
  stopifnot(nrow(background) > 0)
  
  for (state in q3_states) {
    obs_count <- sum(gwas_sites$consensus == state)
    obs_prop <- obs_count / nrow(gwas_sites)
    
    bg_count <- sum(background$consensus == state)
    bg_prop <- bg_count / nrow(background)
    
    cont_table <- matrix(c(obs_count, nrow(gwas_sites) - obs_count,
                           bg_count, nrow(background) - bg_count),
                         nrow = 2, byrow = TRUE)
    
    fisher_test <- fisher.test(cont_table)
    
    results <- rbind(results, data.frame(
      gwas_class = gwas_class,
      q3_state = state,
      n_gwas_sites = nrow(gwas_sites),
      obs_count = obs_count,
      obs_prop = obs_prop,
      bg_prop = bg_prop,
      fold_enrichment = obs_prop / bg_prop,
      p_value = fisher_test$p.value,
      odds_ratio = fisher_test$estimate
    ))
  }
}

results$p_adj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$p_adj), ]
print(results)

# Show that structural properties correlate even when p-values don't

# 1. Correlation of p-values by structure
cat("\n=== P-value correlations by structure ===\n")
for (state in c("H", "E", "C")) {
  subset_data <- high_cons[high_cons$consensus == state, ]
  cor_val <- cor(-log10(subset_data$P_aa_only), -log10(subset_data$P_aa_with_pcs))
  cat(sprintf("%s: r = %.3f\n", state, cor_val))
}

# 2. Visual comparison
par(mfrow = c(2, 2))

# Overall correlation
plot(-log10(high_cons$P_aa_only), -log10(high_cons$P_aa_with_pcs),
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.2),
     xlab = "-log10(P without control)",
     ylab = "-log10(P with control)",
     main = sprintf("Overall correlation\nr = %.3f", 
                    cor(-log10(high_cons$P_aa_only), -log10(high_cons$P_aa_with_pcs))))
abline(0, 1, col = "red", lty = 2)

# By structure
colors <- c("H" = "steelblue", "E" = "orange", "C" = "purple")
plot(-log10(high_cons$P_aa_only), -log10(high_cons$P_aa_with_pcs),
     col = colors[high_cons$consensus],
     pch = 16, cex = 0.5,
     xlab = "-log10(P without control)",
     ylab = "-log10(P with control)",
     main = "Colored by structure")
legend("topleft", legend = c("H", "E", "C"), 
       col = colors, pch = 16, bty = "n")
abline(0, 1, col = "red", lty = 2)

# 3. Show that structure enrichment is consistent across both tests
cat("\n=== Structure enrichment in top sites (both tests) ===\n")

# Top sites from each test
top_no_ctrl <- high_cons[high_cons$P_aa_only < 0.01, ]
top_with_ctrl <- high_cons[high_cons$P_aa_with_pcs < 0.01, ]

cat("\nTop sites (P < 0.01) without control:\n")
print(round(prop.table(table(top_no_ctrl$consensus)), 3))

cat("\nTop sites (P < 0.01) with control:\n")
print(round(prop.table(table(top_with_ctrl$consensus)), 3))

cat("\nBackground (all sites):\n")
print(round(prop.table(table(high_cons$consensus)), 3))

# 4. Barplot comparison
struct_comparison <- rbind(
  all = prop.table(table(high_cons$consensus)),
  top_no_ctrl = prop.table(table(top_no_ctrl$consensus)),
  top_with_ctrl = prop.table(table(top_with_ctrl$consensus))
)

barplot(struct_comparison,
        beside = TRUE,
        col = c("lightgray", "gold", "darkred"),
        legend.text = c("All sites", "Top (no ctrl)", "Top (with ctrl)"),
        args.legend = list(x = "topright", bty = "n"),
        ylab = "Proportion",
        xlab = "Structure",
        main = "Structure enrichment is consistent\nacross both significance tests")

par(mfrow = c(1, 1))


# ---- plastome mapping /circos plot ---- 
library(data.table)
library(Biostrings)
library(genbankr)  # for parsing GenBank files

 ---- CONFIG ----
ARABIDOPSIS_ID <- "AP000423.1"
GBF_FILE <- "data/gbfs/AP0004231fa.gbf"  # adjust if different naming

# ---- LOAD GWAS RESULTS ----
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Aligned_Position = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs,
      P_aa_only = m$P_aa_only
    )
  }))
}))
stopifnot(nrow(gwas_results) > 0)

# ---- PARSE GENBANK FILE ----
# Extract gene coordinates from Arabidopsis chloroplast
gb <- readGenBank(GBF_FILE)
features <- otherFeatures(gb)  # or genes(gb), cds(gb) depending on structure

# Extract CDS features with gene names and coordinates
# Adjust parsing based on actual GBF structure
gene_coords <- data.table(
  gene = mcols(features)$gene,
  start = start(features),
  end = end(features),
  strand = as.character(strand(features))
)
gene_coords <- gene_coords[!is.na(gene)]
setkey(gene_coords, gene)

# ---- MAP ALIGNMENT POSITIONS TO ARABIDOPSIS SEQUENCE ----
aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", 
                        full.names = TRUE)

mapping_list <- list()

for (aln_file in aln_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(aln_file))
  
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(aln)) next
  
  names(aln) <- sub("\\|.*", "", names(aln))
  
  # Find Arabidopsis sequence
  at_idx <- which(names(aln) == ARABIDOPSIS_ID)
  if (length(at_idx) == 0) {
    message("Gene ", gene, ": Arabidopsis not found in alignment")
    next
  }
  
  at_seq <- as.character(aln[[at_idx]])
  chars <- strsplit(at_seq, "")[[1]]
  
  # Build mapping: alignment position -> Arabidopsis ungapped position
  ungapped_pos <- 0
  aln_to_at <- integer(length(chars))
  
  for (i in seq_along(chars)) {
    if (chars[i] != "-") {
      ungapped_pos <- ungapped_pos + 1
      aln_to_at[i] <- ungapped_pos
    } else {
      aln_to_at[i] <- NA  # gap in Arabidopsis
    }
  }
  
  mapping_list[[gene]] <- data.table(
    Gene = gene,
    Aligned_Position = seq_along(chars),
    At_AA_Position = aln_to_at
  )
}

pos_mapping <- rbindlist(mapping_list)
stopifnot(nrow(pos_mapping) > 0)

# ---- MERGE GWAS WITH MAPPING ----
gwas_mapped <- merge(gwas_results, pos_mapping, 
                     by = c("Gene", "Aligned_Position"),
                     all.x = TRUE)

# ---- CONVERT AA POSITION TO GENOMIC COORDINATE ----
# For each gene, AA position * 3 gives codon start within the gene
# Then add gene start coordinate

gwas_mapped <- merge(gwas_mapped, gene_coords, 
                     by.x = "Gene", by.y = "gene",
                     all.x = TRUE)

# Calculate genomic position (codon start)
gwas_mapped[, Genomic_Position := ifelse(
  strand == "+",
  start + (At_AA_Position - 1) * 3,
  end - (At_AA_Position - 1) * 3 - 2  # for minus strand, count from end
)]

# Add significance classification
gwas_thresh <- quantile(gwas_mapped$P_aa_with_pcs, 0.05, na.rm = TRUE)
gwas_mapped[, gwas_hit := P_aa_with_pcs < gwas_thresh]

# ---- OUTPUT ----
gwas_final <- gwas_mapped[!is.na(Genomic_Position), .(
  Gene, Aligned_Position, At_AA_Position, 
  Genomic_Position, strand,
  P_aa_with_pcs, P_aa_only, gwas_hit
)]

setorder(gwas_final, Genomic_Position)

message("Mapped ", nrow(gwas_final), " positions to Arabidopsis genome")
message("GWAS hits: ", sum(gwas_final$gwas_hit))

saveRDS(gwas_final, "results/gwas_arabidopsis_mapped.rds")
fwrite(gwas_final, "results/gwas_arabidopsis_mapped.csv")

# ---- QUICK VISUALIZATION ----
# Manhattan-style plot on genomic coordinates
plot(gwas_final$Genomic_Position, -log10(gwas_final$P_aa_with_pcs),
     pch = 16, cex = 0.5,
     col = ifelse(gwas_final$gwas_hit, "red", "grey50"),
     xlab = "Arabidopsis chloroplast position (bp)",
     ylab = "-log10(P)",
     main = "GWAS mapped to Arabidopsis chloroplast")

# Add gene labels
gene_mids <- gwas_final[, .(mid = median(Genomic_Position)), by = Gene]
text(gene_mids$mid, par("usr")[4] * 0.95, gene_mids$Gene, 
     cex = 0.6, srt = 45, adj = 1)


