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

final_res$sig_class <- case_when(
  final_res$P_aa_only < bf_threshold & final_res$P_aa_with_pcs < bf_threshold ~ "sig_both",
  final_res$P_aa_only < bf_threshold & final_res$P_aa_with_pcs >= bf_threshold ~ "sig_lost",
  final_res$P_aa_only >= bf_threshold & final_res$P_aa_with_pcs < bf_threshold ~ "sig_gained",
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

sites_df$sig_class <- case_when(
  sites_df$P_aa_only < bf_threshold & sites_df$P_aa_with_pcs < bf_threshold ~ "sig_both",
  sites_df$P_aa_only < bf_threshold & sites_df$P_aa_with_pcs >= bf_threshold ~ "sig_lost",
  sites_df$P_aa_only >= bf_threshold & sites_df$P_aa_with_pcs < bf_threshold ~ "sig_gained",
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





