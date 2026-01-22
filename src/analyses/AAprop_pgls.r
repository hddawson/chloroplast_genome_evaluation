#!/usr/bin/env Rscript
# ==============================================================================
# PGLS: Temperature ~ Amino Acid Proportions by Site Significance Class
# ==============================================================================

library(ape)
library(data.table)
library(caper)
library(arrow)

# ==== CONFIGURATION ====
aln_file <- "raxml_input/superaa_collapsed.fasta"
data_file <- "data/processed_data.parquet"
gwas_dir <- "results/residue_models_triple/"
partition_map_file <- "raxml_input/partitionMap.rds"
tree_file <- "raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"

stopifnot(file.exists(aln_file), file.exists(data_file), file.exists(tree_file))

# ==== LOAD DATA ====

cat("Loading alignment...\n")
aln <- read.FASTA(aln_file, type = "AA")
aln_mat <- do.call(rbind, lapply(aln, function(x) rawToChar(x, multiple = TRUE)))
rownames(aln_mat) <- names(aln)
stopifnot(is.matrix(aln_mat))
cat("Alignment:", nrow(aln_mat), "taxa,", ncol(aln_mat), "positions\n")

cat("Loading phenotype data...\n")
data <- as.data.table(read_parquet(data_file))
pheno <- data$pheno_wc2.1_2.5m_bio_8_p50
names(pheno) <- data$ID

cat("Loading tree...\n")
tree <- read.tree(tree_file)
# Drop internal node labels if present
tree$node.label <- NULL

# ==== LOAD GWAS RESULTS & CLASSIFY SITES ====

cat("Loading GWAS results...\n")
partition_map <- readRDS(partition_map_file)
setDT(partition_map)

model_files <- list.files(gwas_dir, pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_effects <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    if (is.null(m$effects) || nrow(m$effects) == 0) return(NULL)
    data.table(
      Gene = m$Gene,
      GenePos = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs,
      P_aa_only = m$P_aa_only
    )
  }), fill = TRUE)
}), fill = TRUE)

gwas_effects <- unique(gwas_effects)
gwas_effects <- merge(gwas_effects, partition_map[, .(Gene, GenePos, GlobalPos)],
                      by = c("Gene", "GenePos"), all.x = TRUE)
gwas_effects <- gwas_effects[!is.na(GlobalPos)]

# Define thresholds
thresh_control <- quantile(gwas_effects$P_aa_with_pcs, 0.05, na.rm = TRUE)
thresh_nocontrol <- quantile(gwas_effects$P_aa_only, 0.20, na.rm = TRUE)

gwas_effects[, sig_class := fcase(
  P_aa_with_pcs < thresh_control & P_aa_only < thresh_nocontrol, "sig_both",
  P_aa_with_pcs < thresh_control & P_aa_only >= thresh_nocontrol, "sig_control",
  P_aa_with_pcs >= thresh_control & P_aa_only < thresh_nocontrol, "sig_nocontrol",
  default = "not_sig"
)]

cat("\nSite classification:\n")
print(gwas_effects[, .N, by = sig_class][order(-N)])

# ==== CALCULATE AA PROPORTIONS BY SIG CLASS ====

aa_chars <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

calc_aa_props <- function(aln_mat, positions) {
  sub_mat <- aln_mat[, positions, drop = FALSE]
  
  props <- t(apply(sub_mat, 1, function(row) {
    counts <- table(factor(row, levels = aa_chars))
    counts / sum(counts)
  }))
  
  colnames(props) <- paste0("p_", aa_chars)
  dt <- as.data.table(props)
  dt[, taxon := rownames(aln_mat)[1:nrow(dt)]]  # or simpler:
  # dt <- data.table(taxon = rownames(sub_mat), props)
  dt
}

cat("\nCalculating AA proportions by sig_class...\n")

sig_classes <- c("sig_both", "sig_control", "sig_nocontrol", "not_sig")
aa_props_list <- list()

for (sc in sig_classes) {
  pos <- gwas_effects[sig_class == sc, GlobalPos]
  pos <- pos[pos <= ncol(aln_mat)]  # ensure valid positions
  
  if (length(pos) < 10) {
    cat("Skipping", sc, "- too few positions\n")
    next
  }
  
  cat(sc, ":", length(pos), "positions\n")
  props <- calc_aa_props(aln_mat, pos)
  props[, sig_class := sc]
  aa_props_list[[sc]] <- props
}

# ==== PREPARE PGLS DATA ====

cat("\nPreparing PGLS data...\n")

# Find common taxa
common_taxa <- Reduce(intersect, list(
  rownames(aln_mat),
  names(pheno)[!is.na(pheno)],
  tree$tip.label
))
cat("Common taxa:", length(common_taxa), "\n")
stopifnot(length(common_taxa) > 100)

# Prune tree
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))
stopifnot(length(tree_pruned$tip.label) == length(common_taxa))

pinales <- data$ID[grep("Pinales", data$Order)]
pinales_in_tree <- intersect(pinales, tree_pruned$tip.label)
tree_pruned <- root(tree_pruned, outgroup = pinales_in_tree, resolve.root = TRUE)

summary(tree_pruned)
is.rooted(tree_pruned)
# ==== RUN PGLS BY SIG CLASS ====

cat("\n=== PGLS Results ===\n")

pgls_results <- list()

for (sc in names(aa_props_list)) {
  cat("\n---", sc, "---\n")
  
  props <- aa_props_list[[sc]]
  props <- props[taxon %in% common_taxa]
  
  # Add phenotype
  props[, temp := pheno[taxon]]
  props <- props[!is.na(temp)]
  
  # Remove constant columns (zero variance)
  aa_cols <- paste0("p_", aa_chars)
  var_check <- props[, lapply(.SD, var, na.rm = TRUE), .SDcols = aa_cols]
  keep_cols <- aa_cols[unlist(var_check) > 1e-10]
  
  if (length(keep_cols) < 2) {
    cat("Too few variable AA columns, skipping\n")
    next
  }
  
  # Build formula - drop one AA to avoid collinearity (proportions sum to 1)
  keep_cols <- keep_cols[-1]  # drop first
  form <- as.formula(paste("temp ~", paste(keep_cols, collapse = " + ")))
  
  # Prepare comparative data
  df <- as.data.frame(props)
  rownames(df) <- df$taxon
  
  comp_data <- comparative.data(
    phy = tree_pruned,
    data = df,
    names.col = taxon,
    vcv = TRUE,
    warn.dropped = FALSE
  )
  
  # Fit PGLS
  pgls_fit <- tryCatch({
    pgls(form, data = comp_data, lambda = "ML")
  }, error = function(e) {
    cat("PGLS failed:", e$message, "\n")
    NULL
  })
  
  if (is.null(pgls_fit)) next
  
  # Summarize
  summ <- summary(pgls_fit)
  cat("Lambda:", pgls_fit$param["lambda"], "\n")
  cat("R-squared:", summ$r.squared, "\n")
  cat("F-statistic:", summ$fstatistic[1], "on", summ$fstatistic[2], "and", summ$fstatistic[3], "DF\n")
  cat("Model p-value:", pf(summ$fstatistic[1], summ$fstatistic[2], summ$fstatistic[3], lower.tail = FALSE), "\n")
  
  # Store significant coefficients
  coefs <- as.data.table(summ$coefficients, keep.rownames = "term")
  setnames(coefs, c("term", "estimate", "se", "t", "p"))
  sig_coefs <- coefs[p < 0.05 & term != "(Intercept)"]
  
  if (nrow(sig_coefs) > 0) {
    cat("\nSignificant AA effects (p < 0.05):\n")
    print(sig_coefs[order(p)])
  }
  
  pgls_results[[sc]] <- list(
    model = pgls_fit,
    summary = summ,
    coefs = coefs,
    lambda = pgls_fit$param["lambda"],
    r2 = summ$r.squared,
    n_sites = gwas_effects[sig_class == sc, .N]
  )
}
saveRDS(pgls_results, "data/tmp/pgls_aa_props_results.rds")
# ==== SUMMARY TABLE ====

cat("\n=== Summary Across Significance Classes ===\n")

summary_dt <- rbindlist(lapply(names(pgls_results), function(sc) {
  res <- pgls_results[[sc]]
  fstat <- summary(res$model)$fstatistic
  data.table(
    sig_class = sc,
    n_sites = res$n_sites,
    lambda = round(res$lambda, 3),
    r_squared = round(res$r2, 4),
    f_stat = round(fstat[1], 2),
    model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE),
    n_sig_aa = res$coefs[p < 0.05 & term != "(Intercept)", .N]
  )
}))

print(summary_dt[order(-r_squared)])

# ==== SAVE ====

saveRDS(pgls_results, "data/tmp/pgls_aa_props_results.rds")
cat("\nResults saved to data/tmp/pgls_aa_props_results.rds\n")


#read results

pgls_results <- readRDS("data/tmp/pgls_aa_props_results.rds")

#extract coefficients for each amino acid, make a heatmap of them 
library(ggplot2)
library(reshape2)

coef_dt <- rbindlist(lapply(names(pgls_results), function(sc) {
  res <- pgls_results[[sc]]
  coefs <- res$coefs
  coefs[, sig_class := sc]
  coefs
}))

coef_dt <- coef_dt[term != "(Intercept)"]
head(coef_dt)
summary(coef_dt)

ggplot(coef_dt, aes(x = term, y = estimate, fill = sig_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "PGLS Coefficients for Amino Acid Proportions by Significance Class",
       x = "Amino Acid",
       y = "Coefficient Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#heatmap of relative coefficient strengths
coef_mat <- dcast(coef_dt, term ~ sig_class, value.var = "estimate", fill = 0)
rownames(coef_mat) <- coef_mat$term
names_vec <- coef_mat$term
coef_mat$term <- NULL
coef_mat_scaled <- as.data.frame(scale(coef_mat))
rownames(coef_mat_scaled) <- names_vec

pheatmap(coef_mat, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Heatmap of  PGLS Coefficients for Amino Acid Proportions")


library(Peptides)
data(AAdata)

hydrophobicity_scale <- AAdata$Hydrophobicity$Tanford
#    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
# 0.62 -2.53 -0.78 -0.09  0.29 -0.85 -0.74  0.48 -0.40  1.38  1.53 -1.50  0.64  1.19  0.12 -0.18 -0.05  0.81  0.26  1.80 

#does hydrophobicity correlate with coefficient strength?
coef_dt[, aa := sub("p_", "", term)]
coef_dt[, hydrophobicity := hydrophobicity_scale[aa]]
ggplot(coef_dt, aes(x = hydrophobicity, y = estimate, color = sig_class)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Correlation of Amino Acid Hydrophobicity with PGLS Coefficients",
       x = "Hydrophobicity (Tanford Scale)",
       y = "PGLS Coefficient Estimate")

#correlations please 
cor_dt <- coef_dt[, .(correlation = cor(hydrophobicity, estimate, use = "complete.obs")), by = sig_class]
print(cor_dt)

#multivariate analysis 

kideras <- AAdata$kideraFactors
#$KF1
#    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
#-1.56  0.22  1.14  0.58  0.12 -0.47 -1.45  1.46 -0.41 -0.73 -1.04 -0.34 -1.40 -0.21  2.06  0.81  0.26  0.30  1.38 -0.74 

#$KF2
#    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
#-1.67  1.27 -0.07 -0.22 -0.89  0.24  0.19 -1.96  0.52 -0.16  0.00  0.82  0.18  0.98 -0.33 -1.08 -0.70  2.10  1.48 -0.71 

#$KF3
#    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
#-0.97  1.37 -0.12 -1.58  0.45  0.07 -1.61 -0.23 -0.28  1.79 -0.24 -0.23 -0.42 -0.36 -1.15  0.16  1.21 -0.72  0.80  2.04 
#etc. until 10 - do any correlate?

for (i in 1:10) {
  factor_scale <- kideras[[paste0("KF", i)]]
  coef_dt[, kidera := factor_scale[aa]]
  
  p <- ggplot(coef_dt, aes_string(x = "kidera", y = "estimate", color = "sig_class")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    labs(title = paste("Correlation of Amino Acid Kidera Factor", i, "with PGLS Coefficients"),
         x = paste("Kidera Factor", i),
         y = "PGLS Coefficient Estimate")
  print(p)
  
  cor_dt <- coef_dt[, .(correlation = cor(kidera, estimate, use = "complete.obs")), by = sig_class]
  cat("\nKidera Factor", i, "correlations:\n")
  print(cor_dt)
}


coef_dt <- rbindlist(lapply(names(pgls_results), function(sc) {
  res <- pgls_results[[sc]]
  coefs <- res$coefs
  coefs[, sig_class := sc]
  coefs
}))
coef_dt <- coef_dt[term != "(Intercept)"]

# Create estimate matrix
coef_mat <- dcast(coef_dt, term ~ sig_class, value.var = "estimate", fill = 0)
rn <- coef_mat$term
coef_mat[, term := NULL]
coef_mat <- as.matrix(coef_mat)
rownames(coef_mat) <- rn

# Create p-value matrix for annotations
p_mat <- dcast(coef_dt, term ~ sig_class, value.var = "p", fill = 1)
p_mat[, term := NULL]
p_mat <- as.matrix(p_mat)
rownames(p_mat) <- rn

# Convert p-values to significance stars
sig_stars <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_stars[p_mat < 0.05] <- "*"
sig_stars[p_mat < 0.01] <- "**"
sig_stars[p_mat < 0.001] <- "***"
rownames(sig_stars) <- rownames(p_mat)
colnames(sig_stars) <- colnames(p_mat)

max_abs <- max(abs(coef_mat))
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(coef_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         display_numbers = sig_stars,
         fontsize_number = 12,
         breaks = breaks,
         color = colors,
         main = "PGLS Coefficients for AA Proportions\n(* p<0.05, ** p<0.01, *** p<0.001)")

coef_mat
sapply(pgls_results, function(res) res$coefs[term == "(Intercept)", estimate])

# ---- corr b w kfs and effects by sig class ----
coef_dt <- rbindlist(lapply(names(pgls_results), function(sc) {
  res <- pgls_results[[sc]]
  coefs <- res$coefs
  coefs[, sig_class := sc]
  coefs
}))
coef_dt <- coef_dt[term != "(Intercept)"]
coef_dt[, aa := sub("p_", "", term)]

# KF annotations
kf_labels <- c(
  KF1 = "KF1: Helix/bend",
  KF2 = "KF2: Side-chain size",
  KF3 = "KF3: Extended structure",
  KF4 = "KF4: Hydrophobicity",
  KF5 = "KF5: Double-bend",
  KF6 = "KF6: Partial specific vol",
  KF7 = "KF7: Flat extended",
  KF8 = "KF8: Alpha region",
  KF9 = "KF9: pK-C",
  KF10 = "KF10: Surrounding hydro"
)

# Calculate correlations and p-values for each KF x sig_class
cor_results <- rbindlist(lapply(1:10, function(i) {
  factor_scale <- kideras[[paste0("KF", i)]]
  coef_dt[, kidera := factor_scale[aa]]
  
  coef_dt[, {
    ct <- cor.test(kidera, estimate, use = "complete.obs")
    .(KF = paste0("KF", i), correlation = ct$estimate, p = ct$p.value)
  }, by = sig_class]
}))

# Create correlation matrix
cor_mat <- dcast(cor_results, KF ~ sig_class, value.var = "correlation")
rn <- cor_mat$KF
cor_mat[, KF := NULL]
cor_mat <- as.matrix(cor_mat)
rownames(cor_mat) <- kf_labels[rn]

# Create p-value matrix
p_mat <- dcast(cor_results, KF ~ sig_class, value.var = "p")
p_mat[, KF := NULL]
p_mat <- as.matrix(p_mat)
rownames(p_mat) <- rn

# Significance stars
sig_stars <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_stars[p_mat < 0.05] <- "*"
sig_stars[p_mat < 0.01] <- "**"
sig_stars[p_mat < 0.001] <- "***"

# Symmetric color scale centered on 0
max_abs <- max(abs(cor_mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("orange", "white", "green"))(100)

pheatmap(cor_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = sig_stars,
         fontsize_number = 12,
         breaks = breaks,
         color = colors,
         main = "Correlation of Kidera Factors with AA PGLS Coefficients\n(* p<0.05, ** p<0.01, *** p<0.001)")


# Build long-form data with all KF values
plot_dt <- rbindlist(lapply(1:10, function(i) {
  factor_scale <- kideras[[paste0("KF", i)]]
  dt <- copy(coef_dt)
  dt[, kidera := factor_scale[aa]]
  dt[, KF := factor(kf_labels[paste0("KF", i)], levels = kf_labels)]
  dt
}))

# Calculate correlations and p-values for annotations
cor_results <- plot_dt[, {
  ct <- cor.test(kidera, estimate, use = "complete.obs")
  .(r = round(ct$estimate, 2), p = ct$p.value)
}, by = .(sig_class, KF)]

cor_results[, label := paste0("r=", r, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))]

# Merge for label positioning
label_pos <- plot_dt[, .(x = max(kidera), y = max(estimate)), by = .(sig_class, KF)]
cor_results <- merge(cor_results, label_pos, by = c("sig_class", "KF"))

ggplot(plot_dt, aes(x = kidera, y = estimate)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  geom_text(data = cor_results, aes(x = x, y = y, label = label), 
            hjust = 1, vjust = 1, size = 2.5) +
  facet_grid(KF ~ sig_class, scales = "free") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "grey80")) +
  labs(x = "Kidera Factor Value", 
       y = "PGLS Coefficient Estimate",
       title = "AA PGLS Coefficients vs Kidera Factors by Significance Class")


cor_results <- rbindlist(lapply(1:10, function(i) {
  factor_scale <- kideras[[paste0("KF", i)]]
  coef_dt[, kidera := factor_scale[aa]]
  
  coef_dt[, {
    ct <- cor.test(kidera, estimate, use = "complete.obs")
    .(KF = paste0("KF", i), r = ct$estimate, p = ct$p.value)
  }, by = sig_class]
}))

# Get top 3 most significant
top3 <- cor_results[order(p)][1:3]
print(top3)

# Generate plots
plots <- lapply(1:3, function(j) {
  row <- top3[j]
  kf_name <- row$KF
  sc <- row$sig_class
  
  factor_scale <- kideras[[kf_name]]
  dt <- coef_dt[sig_class == sc]
  dt[, kidera := factor_scale[aa]]
  
  ct <- cor.test(dt$kidera, dt$estimate)
  
  ggplot(dt, aes(x = kidera, y = estimate)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "lightblue") +
    geom_text(aes(label = aa), hjust = -0.3, vjust = -0.3, size = 3) +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("r = %.3f\np = %.4f", ct$estimate, ct$p.value),
             hjust = 1.1, vjust = 1.5, size = 4) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = paste0(kf_name, ": ", kf_labels[kf_name]),
         y = "PGLS Coefficient Estimate (°C)",
         title = paste0("Sig class: ", sc))
})

# Print or arrange
library(patchwork)
plots[[1]] / plots[[2]] / plots[[3]]
