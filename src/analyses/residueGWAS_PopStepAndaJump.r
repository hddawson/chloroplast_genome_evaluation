library(arrow)
library(data.table)
library(Biostrings)

# LOAD COMMON DATA
data <- as.data.table(read_parquet("data/processed_data.parquet"))
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"

ev_pcs <- readRDS("data/tmp/majMinor_aln_pca.rds")
pvar <- ev_pcs$sdev^2 / sum(ev_pcs$sdev^2)
pvar[1:10]
aln_index <- read_parquet("data/tmp/majMinor_aln.pq")

embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))

setDT(embeds_with_mds)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

# Prepare PC table - keep ALL PCs
pcs_IDS <- aln_index$index
scores <- as.data.table(ev_pcs$x)
max_pcs <- ncol(scores)
setnames(scores, paste0("PC", seq_len(max_pcs)))
scores[, ID := pcs_IDS]

message("Total PCs available: ", max_pcs)

# Pre-join phenotype + all PCs
all_pc_names <- paste0("PC", seq_len(max_pcs))
pheno_pcs_full <- data[, .(ID, pheno = get(pheno_col))][scores, on = "ID", nomatch = 0]
stopifnot(nrow(pheno_pcs_full) > 0)

# PREPARE FILE LISTS

aln_files <- list.files("data/tmp/alignedGenes/", pattern = "_AA_aligned\\.fasta$", full.names = TRUE)
get_gene <- function(path) sub("_AA_aligned\\.fasta", "", basename(path))
genes_to_process <- get_gene(aln_files)
stopifnot(length(genes_to_process) > 0)

message("Processing ", length(genes_to_process), " genes")


# N_PCS EXPERIMENT SETUP

# Define n_pcs values: 0, 250, 500, ..., up to max_pcs
n_pcs_values <- seq(0, max_pcs, by = 250)
if (max(n_pcs_values) < max_pcs) {
  n_pcs_values <- c(n_pcs_values, max_pcs)
}

message("Testing ", length(n_pcs_values), " n_pcs values: ", 
        paste(range(n_pcs_values), collapse = " to "))

# Create output directory structure
base_out_dir <- "results/npc_experiment"
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)


# MAIN LOOP: iterate over n_pcs values


for (n_pcs in n_pcs_values) {
  
  message("\n========== n_pcs = ", n_pcs, " ==========")
  
  out_dir <- file.path(base_out_dir, paste0("npc_", n_pcs))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Subset PCs for this run
  if (n_pcs == 0) {
    pc_names_use <- character(0)
  } else {
    pc_names_use <- paste0("PC", seq_len(n_pcs))
  }
  
  for (gene in genes_to_process) {
    
    out_path <- file.path(out_dir, paste0(gene, "_effects.rds"))
    
    if (file.exists(out_path)) {
      message("Skipping ", gene, " (output already exists)")
      next
    }
    
    # Get clean IDs for this gene
    clean_ids_gene <- clean_ids_by_gene[Gene == gene, ID]
    if (length(clean_ids_gene) == 0) {
      message("Skipping ", gene, ": no clean IDs")
      next
    }
    
    aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
    
    aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
    if (is.null(aln) || length(aln) == 0) {
      message("Skipping ", gene, ": could not read alignment")
      next
    }
    
    names(aln) <- sub("\\|.*", "", names(aln))
    aln_mat <- as.matrix(aln)
    seq_ids <- names(aln)
    
    common_ids <- intersect(intersect(seq_ids, pheno_pcs_full$ID), clean_ids_gene)
    if (length(common_ids) < 8500L) {
      message("Skipping ", gene, ": insufficient overlap (n=", length(common_ids), ")")
      next
    }
    
    aln_mat <- aln_mat[match(common_ids, seq_ids), , drop = FALSE]
    pheno_pcs_sub <- pheno_pcs_full[ID %in% common_ids]
    setkey(pheno_pcs_sub, ID)
    pheno_pcs_sub <- pheno_pcs_sub[match(common_ids, ID)]
    
    y <- pheno_pcs_sub$pheno
    
    # Prepare PC matrix (or NULL if n_pcs == 0)
    if (n_pcs > 0) {
      X_pcs <- scale(as.matrix(pheno_pcs_sub[, ..pc_names_use]))
    } else {
      X_pcs <- NULL
    }
    
    # -------------------------------------------------------------------
    # LOOP OVER POSITIONS
    # -------------------------------------------------------------------
    
    positions <- seq_len(ncol(aln_mat))
    results_list <- vector("list", length(positions))
    
    for (pos in positions) {
      
      residues <- aln_mat[, pos]
      residue_table <- table(residues)
      
      if (length(residue_table) < 2L || all(names(residue_table) == "-")) next
      
      non_gap <- residues != "-"
      if (sum(non_gap) < 8500L) next
      
      res_factor <- factor(residues)
      X_aa <- model.matrix(~ res_factor - 1)
      
      var0 <- which(apply(X_aa, 2, var) == 0)
      if (length(var0) > 0) X_aa <- X_aa[, -var0, drop = FALSE]
      if (ncol(X_aa) == 0) next
      
      # -----------------------------------------------------------------
      # FIT MODELS
      # -----------------------------------------------------------------
      
      fit_null <- tryCatch(lm(y ~ 1), error = function(e) NULL)
      if (is.null(fit_null)) next
      
      fit_aa <- tryCatch(lm(y ~ X_aa), error = function(e) NULL)
      if (is.null(fit_aa)) next
      
      if (n_pcs > 0) {
        fit_pcs <- tryCatch(lm(y ~ X_pcs), error = function(e) NULL)
        fit_full <- tryCatch(lm(y ~ X_aa + X_pcs), error = function(e) NULL)
        if (is.null(fit_pcs) || is.null(fit_full)) next
        
        r2_pcs <- summary(fit_pcs)$r.squared
        r2_full <- summary(fit_full)$r.squared
        p_pcs_only <- tryCatch(anova(fit_null, fit_pcs)[2, "Pr(>F)"], error = function(e) NA_real_)
        p_aa_with_pcs <- tryCatch(anova(fit_pcs, fit_full)[2, "Pr(>F)"], error = function(e) NA_real_)
        
        coef_summary <- summary(fit_full)$coefficients
      } else {
        r2_pcs <- NA_real_
        r2_full <- summary(fit_aa)$r.squared
        p_pcs_only <- NA_real_
        p_aa_with_pcs <- NA_real_
        
        coef_summary <- summary(fit_aa)$coefficients
      }
      
      r2_aa <- summary(fit_aa)$r.squared
      p_aa_only <- tryCatch(anova(fit_null, fit_aa)[2, "Pr(>F)"], error = function(e) NA_real_)
      
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
        n_pcs = n_pcs,
        R2_aa = r2_aa,
        R2_pcs = r2_pcs,
        R2_full = r2_full,
        P_aa_only = p_aa_only,
        P_pcs_only = p_pcs_only,
        P_aa_with_pcs = p_aa_with_pcs,
        residue_counts = as.list(residue_table),
        effects = effects_dt
      )
    }
    
    results_list <- Filter(Negate(is.null), results_list)
    
    if (length(results_list) > 0) {
      saveRDS(results_list, out_path)
      message("  ", gene, ": ", length(results_list), " positions")
    }
  }
}

message("\n========== EXPERIMENT COMPLETE ==========")
message("Results saved to: ", base_out_dir)

# ---- analysis ---- 
library(data.table)


# LOAD RESULTS


base_out_dir <- "results/npc_experiment"
npc_dirs <- list.dirs(base_out_dir, recursive = FALSE)
stopifnot(length(npc_dirs) > 0)

# Extract n_pcs values from directory names
get_npc <- function(d) as.integer(sub("npc_", "", basename(d)))
npc_values <- sort(sapply(npc_dirs, get_npc))
message("Found results for n_pcs: ", paste(npc_values, collapse = ", "))

# Load all results into one table
load_npc_results <- function(n_pcs) {
  dir_path <- file.path(base_out_dir, paste0("npc_", n_pcs))
  rds_files <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) return(NULL)
  
  all_res <- rbindlist(lapply(rds_files, function(f) {
    res_list <- readRDS(f)
    rbindlist(lapply(res_list, function(r) {
      data.table(
        Gene = r$Gene,
        Position = r$Aligned_Position,
        n_pcs = r$n_pcs,
        N = r$N,
        R2_aa = r$R2_aa,
        R2_pcs = r$R2_pcs,
        R2_full = r$R2_full,
        P_aa_only = r$P_aa_only,
        P_aa_with_pcs = r$P_aa_with_pcs
      )
    }))
  }), fill = TRUE)
  return(all_res)
}

results <- rbindlist(lapply(npc_values, load_npc_results))
stopifnot(nrow(results) > 0)

# Create site ID
results[, site := paste(Gene, Position, sep = "_")]
message("Loaded ", nrow(results), " site-npc combinations")
message("Unique sites: ", uniqueN(results$site))

results <- results[results$n_pcs!=2250]

# ---------------------------------------------------------------------
# 1. NUMBER OF SIGNIFICANT SITES VS N_PCS
# ---------------------------------------------------------------------

alpha <- 0.05

# For n_pcs=0, use P_aa_only; otherwise use P_aa_with_pcs
results[, P_relevant := ifelse(n_pcs == 0, P_aa_only, P_aa_with_pcs)]

sig_counts <- results[, .(
  n_sig_nominal = sum(P_relevant < alpha, na.rm = TRUE),
  n_sig_bonf = sum(P_relevant < alpha / .N, na.rm = TRUE),
  n_total = .N,
  median_P = median(P_relevant, na.rm = TRUE),
  median_R2 = median(R2_full, na.rm = TRUE)
), by = n_pcs][order(n_pcs)]

message("\n=== Significant sites by n_pcs ===")
print(sig_counts)

# Plot
pdf("results/npc_experiment/sig_sites_vs_npcs.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sig_counts$n_pcs, sig_counts$n_sig_nominal, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = "Significant sites (p < 0.05)",
     main = "Significance vs Number PCs in GWAS")
abline(h = sig_counts[n_pcs == max(n_pcs), n_sig_nominal], lty = 2, col = "red")


plot(sig_counts$n_pcs, sig_counts$n_sig_bonf, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = "Significant sites (Bonferroni)",
     main = "Significance vs Number PCs in GWAS")
abline(h = sig_counts[n_pcs == max(n_pcs), n_sig_bonf], lty = 2, col = "red")

dev.off()

par(mfrow=c(2,1))
plot(sig_counts$n_pcs, sig_counts$n_sig_bonf, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = "Significant sites (Bonferroni)",
     main = "Significance vs Number PCs in GWAS")
abline(h = sig_counts[n_pcs == max(n_pcs), n_sig_bonf], lty = 2, col = "red")
barplot(cumsum(pvar)[0:2000], main="Percent Variance Explained By PCs")

# Calculate cumulative variance explained at each n_pcs value
cumvar <- cumsum(pvar)
var_at_npc <- sapply(sig_counts$n_pcs, function(x) {
  if (x == 0) return(0)
  round(cumvar[x] * 100, 1)
})

# Plot with annotations
plot(sig_counts$n_pcs, sig_counts$n_sig_bonf, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = "Significant sites (Bonferroni)",
     main = "Significance vs Number PCs in GWAS")
abline(h = sig_counts[n_pcs == max(n_pcs), n_sig_bonf], lty = 2, col = "red")

# Add variance labels
text(sig_counts$n_pcs, sig_counts$n_sig_bonf, 
     labels = paste0(var_at_npc, "%"), 
     pos = 3, cex = 0.7, col = "blue")
par(mfrow=c(1,1))
# ---------------------------------------------------------------------
# 2. INFLECTION POINT DETECTION
# ---------------------------------------------------------------------

# Simple approach: find where rate of change slows
if (nrow(sig_counts) > 2) {
  delta_sig <- -diff(sig_counts$n_sig_nominal)
  delta_pcs <- diff(sig_counts$n_pcs)
  rate <- delta_sig / delta_pcs
  
  # Inflection = where rate drops substantially
  rate_dt <- data.table(
    n_pcs_from = sig_counts$n_pcs[-nrow(sig_counts)],
    n_pcs_to = sig_counts$n_pcs[-1],
    rate_of_decrease = rate
  )
  message("\n=== Rate of decrease in significant sites ===")
  print(rate_dt)
}

# ---------------------------------------------------------------------
# 3. STABILITY OF EFFECT ESTIMATES ACROSS N_PCS
# ---------------------------------------------------------------------

# Wide format for R2_full correlation
r2_wide <- dcast(results, site ~ n_pcs, value.var = "R2_full")
r2_mat <- as.matrix(r2_wide[, -1, with = FALSE])
rownames(r2_mat) <- r2_wide$site

# Pairwise correlations of R2 across n_pcs levels
r2_cor <- cor(r2_mat, use = "pairwise.complete.obs")
message("\n=== R2 correlation matrix across n_pcs ===")
print(round(r2_cor, 3))

library(pheatmap)
pheatmap(r2_cor)


# Same for -log10(P)
results[, neglogP := -log10(P_relevant + 1e-300)]
p_wide <- dcast(results, site ~ n_pcs, value.var = "neglogP")
p_mat <- as.matrix(p_wide[, -1, with = FALSE])
p_cor <- cor(p_mat, use = "pairwise.complete.obs")
message("\n=== -log10(P) correlation matrix across n_pcs ===")
print(round(p_cor, 3))
pheatmap(p_cor, main = "Correlation between -log10p across num PCs controlled for")

# ---------------------------------------------------------------------
# 4. SITE-LEVEL STABILITY: CORRELATION WITH MAX N_PCS
# ---------------------------------------------------------------------
npc_values <- npc_values[1:9]
max_npc <- max(npc_values)
ref_col <- as.character(max_npc)

stability_vs_max <- data.table(
  n_pcs = npc_values,
  r2_cor_with_max = sapply(as.character(npc_values), function(x) cor(r2_mat[, x], r2_mat[, ref_col], use = "complete.obs")),
  neglogP_cor_with_max = sapply(as.character(npc_values), function(x) cor(p_mat[, x], p_mat[, ref_col], use = "complete.obs"))
)
message("\n=== Correlation with max n_pcs (", max_npc, ") ===")
print(stability_vs_max)

pdf("results/npc_experiment/stability_vs_npcs.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(stability_vs_max$n_pcs, stability_vs_max$r2_cor_with_max, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = paste0("Correlation with n_pcs=", max_npc),
     main = "R² stability", ylim = c(0, 1))

plot(stability_vs_max$n_pcs, stability_vs_max$neglogP_cor_with_max, type = "b", pch = 19,
     xlab = "Number of PCs", ylab = paste0("Correlation with n_pcs=", max_npc),
     main = "-log10(P) stability", ylim = c(0, 1))

dev.off()

# ---------------------------------------------------------------------
# 5. IDENTIFY CONSISTENTLY SIGNIFICANT SITES
# ---------------------------------------------------------------------

# Sites significant at all n_pcs >= some threshold
min_npc_stable <- 500  # adjust as needed
high_npc <- npc_values[npc_values >= min_npc_stable]

consistent_sig <- results[n_pcs %in% high_npc, .(
  n_sig = sum(P_relevant < alpha, na.rm = TRUE),
  n_tested = .N,
  mean_R2 = mean(R2_full, na.rm = TRUE)
), by = site]

consistent_sig[, prop_sig := n_sig / n_tested]
consistent_sites <- consistent_sig[prop_sig == 1][order(-mean_R2)]

message("\n=== Sites significant across all n_pcs >= ", min_npc_stable, " ===")
message("N = ", nrow(consistent_sites))
if (nrow(consistent_sites) > 0) {
  print(head(consistent_sites, 20))
}



# Sites significant at ALL n_pcs levels (0 through 2000)
all_npc <- unique(results$n_pcs)

consistent_all <- results[, .(
  n_sig = sum(P_relevant < alpha, na.rm = TRUE),
  n_tested = .N
), by = site]

consistent_all[, prop_sig := n_sig / n_tested]
always_sig <- consistent_all[prop_sig == 1][order(-n_sig)]

message("Sites significant at ALL ", length(all_npc), " n_pcs levels: ", nrow(always_sig))
print(always_sig)

# SAVE SUMMARY






saveRDS(list(
  sig_counts = sig_counts,
  r2_correlation = r2_cor,
  p_correlation = p_cor,
  stability_vs_max = stability_vs_max,
  consistent_sites = consistent_sites
), file.path(base_out_dir, "analysis_summary.rds"))

message("\n=== Analysis complete ===")
message("PDFs saved to: ", base_out_dir)

