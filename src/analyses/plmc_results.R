library(ggplot2)


couplings <- read.table("results/plmc/rbcL_couplings.txt", 
                        col.names = c("i", "sep1", "j", "sep2", "gap", "score"))

summary(couplings)
couplings <- couplings[, c("i", "j", "score")]

# Filter: top 5% AND long-range (|i-j| >= 5)
couplings$sep <- abs(couplings$i - couplings$j)
long_range <- couplings[couplings$sep >= 5, ]
threshold <- quantile(long_range$score, 0.95)
top_couplings <- long_range[long_range$score >= threshold, ]

stopifnot(nrow(top_couplings) > 0)

n <- max(couplings$i, couplings$j)

# Plot contact map - upper left to bottom right convention
ggplot(top_couplings, aes(i, j)) +
  geom_point(size = 1, shape = 21, fill = "black") +
  geom_point(aes(x = j, y = i), size = 1, shape = 21, fill = "black") +  # mirror
  scale_y_reverse() +  # flip y so (1,1) is top-left
  coord_fixed(xlim = c(1, n), ylim = c(n, 1)) +
  theme_minimal() +
  labs(title = "rbcL top 5% long-range couplings (|i-j| ≥ 5)", 
       x = "position", y = "position")


# APC correction (Average Product Correction)
apc_correct <- function(df) {
  # Calculate mean score for each position
  mean_i <- tapply(df$score, df$i, mean)
  mean_j <- tapply(df$score, df$j, mean)
  overall_mean <- mean(df$score)
  
  df$apc <- mean_i[as.character(df$i)] * mean_j[as.character(df$j)] / overall_mean
  df$cn <- df$score - df$apc  # corrected norm (CN) score
  df
}

couplings <- apc_correct(couplings)

# Filter to |i-j| >= 5 (remove local contacts)
couplings_filt <- couplings[abs(couplings$i - couplings$j) >= 5, ]

# Plot top N couplings
nrow(couplings) * 0.05
top_n <- 11850
top_pairs <- head(couplings[order(-couplings$cn), ], top_n)

ggplot(top_pairs, aes(i, j)) +
  geom_point(aes(color = cn), size = 1) +
  geom_point(aes(x = j, y = i, color = cn), size = 1) +  # mirror
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "rbcL EC contact map", color = "CN score")

ggplot(top_pairs, aes(i, j)) +
  geom_point(aes(color = cn), size = 0.4, alpha = 0.1) +
  geom_point(aes(x = j, y = i, color = cn),
             size = 0.4, alpha = 0.1) +
  scale_color_viridis_c(option = "magma") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "rbcL EC contact map - plmc, top 5% couplings", color = "CN score")


rbcL_couplings <- read.csv("results/evcouplings/rbcL/couplings/rbcL_CouplingScores.csv")
summary(rbcL_couplings)
hist(rbcL_couplings$probability)

hist(rbcL_couplings$fn)
hist(rbcL_couplings$cn)
plot(rbcL_couplings$i, rbcL_couplings$cn)

top_couples <- rbcL_couplings[rbcL_couplings$cn > quantile(rbcL_couplings$cn, 0.99),]
plot(top_couples$i, top_couples$j)

library(data.table)

ec <- fread("results/evcouplings/rbcL/couplings/rbcL_CouplingScores.csv")

# With good data, top 1% should be >> mean
summary(ec$cn)
top_cn <- quantile(ec$cn, 0.99)
mean_cn <- mean(ec$cn)
cat("Top 1% / Mean ratio:", top_cn/mean_cn, "\n")
stopifnot(top_cn/mean_cn > 3)  # Should be much higher with real signal

# Check separation
hist(ec$cn, breaks=100, main="CN distribution (bad if unimodal)")



#!/usr/bin/env Rscript
# Diagnostic script for PLMC/EVcouplings output quality
# Checks whether coupling scores show real coevolutionary signal

library(data.table)
library(ggplot2)

# ============================================================================
# LOAD DATA
# ============================================================================

# Try EVcouplings output first (more complete), fall back to plmc raw output
ec_file <- "results/evcouplings/rbcL/couplings/rbcL_CouplingScores.csv"
plmc_file <- "results/plmc/rbcL_couplings.txt"

if (file.exists(ec_file)) {
  ec <- fread(ec_file)
  cat("Loaded EVcouplings output:", nrow(ec), "pairs\n")
  use_evcouplings <- TRUE
} else if (file.exists(plmc_file)) {
  ec <- fread(plmc_file, col.names = c("i", "sep1", "j", "sep2", "gap", "score"))
  ec <- ec[, .(i, j, score)]
  use_evcouplings <- FALSE
  cat("Loaded plmc raw output:", nrow(ec), "pairs\n")
} else {
  stop("No coupling file found")
}

# ============================================================================
# 1. BASIC STATISTICS
# ============================================================================

cat("\n=== BASIC STATISTICS ===\n")
L <- max(ec$i, ec$j)
cat("Sequence length (L):", L, "\n")
cat("Total pairs:", nrow(ec), "\n")
cat("Expected L*(L-1)/2:", L*(L-1)/2, "\n")

# Score column depends on source
score_col <- if (use_evcouplings) "cn" else "score"
scores <- ec[[score_col]]

cat("\nScore distribution (", score_col, "):\n", sep = "")
print(summary(scores))
cat("SD:", sd(scores), "\n")

# ============================================================================
# 2. SIGNAL-TO-NOISE DIAGNOSTICS
# ============================================================================

cat("\n=== SIGNAL-TO-NOISE CHECKS ===\n")

# Check 1: Top scores should be much higher than mean
top1 <- quantile(scores, 0.99)
top5 <- quantile(scores, 0.95)
mean_score <- mean(scores)

cat("Top 1% threshold:", round(top1, 4), "\n")
cat("Top 5% threshold:", round(top5, 4), "\n")
cat("Mean score:", round(mean_score, 4), "\n")
cat("Top 1% / Mean ratio:", round(top1/mean_score, 2), "\n")
cat("Top 5% / Mean ratio:", round(top5/mean_score, 2), "\n")

# Interpretation
if (top1/mean_score < 3) {
  cat("WARNING: Low top1/mean ratio (<3) suggests weak signal\n")
} else if (top1/mean_score < 5) {
  cat("MARGINAL: Ratio 3-5 suggests moderate signal\n")
} else {
  cat("GOOD: Ratio >5 suggests strong coevolutionary signal\n")
}

# Check 2: Distribution shape - should be heavy-tailed, not normal
skewness <- mean((scores - mean_score)^3) / sd(scores)^3
kurtosis <- mean((scores - mean_score)^4) / sd(scores)^4 - 3

cat("\nSkewness:", round(skewness, 2), "(want >1 for right tail)\n")
cat("Excess kurtosis:", round(kurtosis, 2), "(want >3 for heavy tails)\n")

if (skewness < 1) {
  cat("WARNING: Low skewness suggests no distinct high-scoring tail\n")
}

# ============================================================================
# 3. SEQUENCE SEPARATION ANALYSIS
# ============================================================================

cat("\n=== SEPARATION ANALYSIS ===\n")

ec$sep <- abs(ec$i - ec$j)

# True contacts are typically |i-j| >= 5 (not local backbone)
local <- ec[sep < 5]
longrange <- ec[sep >= 5]

cat("Local pairs (|i-j| < 5):", nrow(local), "\n")
cat("Long-range pairs (|i-j| >= 5):", nrow(longrange), "\n")

# Compare score distributions
cat("\nLocal mean:", round(mean(local[[score_col]]), 4), "\n")
cat("Long-range mean:", round(mean(longrange[[score_col]]), 4), "\n")

# Top pairs should be enriched for medium separation (5-30 typically)
top_pairs <- ec[ec[[score_col]] >= top5]
cat("\nTop 5% pairs - separation distribution:\n")
print(summary(top_pairs$sep))

# ============================================================================
# 4. COVERAGE / EFFECTIVE SEQUENCES (if EVcouplings)
# ============================================================================

if (use_evcouplings) {
  cat("\n=== ALIGNMENT QUALITY METRICS ===\n")
  
  # Check for statistics file
  stats_file <- "results/evcouplings/rbcL/align/rbcL_alignment_statistics.csv"
  if (file.exists(stats_file)) {
    stats <- fread(stats_file)
    print(stats)
    
    # Key metrics
    if ("num_seqs" %in% names(stats)) {
      cat("\nNumber of sequences:", stats$num_seqs, "\n")
    }
    if ("num_eff" %in% names(stats)) {
      neff <- stats$num_eff
      cat("Effective sequences (Neff):", neff, "\n")
      cat("Neff/L ratio:", round(neff/L, 2), "\n")
      
      if (neff/L < 1) {
        cat("CRITICAL: Neff/L < 1 - insufficient diversity!\n")
        cat("  Need Neff/L > 1, ideally > 5 for reliable contacts\n")
      } else if (neff/L < 5) {
        cat("WARNING: Neff/L < 5 - limited diversity\n")
      } else {
        cat("GOOD: Neff/L > 5\n")
      }
    }
  } else {
    cat("Alignment statistics file not found\n")
  }
}

# ============================================================================
# 5. DIAGNOSTIC PLOTS
# ============================================================================

cat("\n=== GENERATING DIAGNOSTIC PLOTS ===\n")

# Plot 1: Score distribution with thresholds
p1 <- ggplot(ec, aes_string(x = score_col)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = top5, color = "red", linetype = "dashed") +
  geom_vline(xintercept = top1, color = "darkred", linetype = "dashed") +
  annotate("text", x = top5, y = Inf, label = "95%", vjust = 2, color = "red") +
  annotate("text", x = top1, y = Inf, label = "99%", vjust = 2, color = "darkred") +
  theme_minimal() +
  labs(title = "Score distribution (should have heavy right tail)",
       subtitle = "Unimodal/symmetric = weak signal",
       x = score_col, y = "count")
p1
# Plot 2: QQ plot - departure from normal indicates signal
theoretical_q <- qnorm(ppoints(length(scores)))
sample_q <- sort(scale(scores))
qq_df <- data.frame(theoretical = theoretical_q, sample = sample_q)

p2 <- ggplot(qq_df, aes(theoretical, sample)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "QQ plot vs normal",
       subtitle = "Upper tail above line = real signal",
       x = "theoretical quantiles", y = "sample quantiles")

# Plot 3: Score vs separation
p3 <- ggplot(ec[sample(.N, min(10000, .N))], aes_string(x = "sep", y = score_col)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "Score vs sequence separation",
       subtitle = "Vertical line = |i-j|=5 cutoff",
       x = "|i - j|", y = score_col)

# Plot 4: Contact map of top pairs
top_plot <- ec[ec[[score_col]] >= top5]
p4 <- ggplot(top_plot, aes_string(x = "i", y = "j", color = score_col)) +
  geom_point(size = 0.3, alpha = 0.5) +
  geom_point(aes_string(x = "j", y = "i", color = score_col), size = 0.3, alpha = 0.5) +
  scale_color_viridis_c(option = "magma") +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Contact map (top 5%)",
       subtitle = "Diagonal band only = no long-range signal")

# Save plots
pdf("coupling_diagnostics.pdf", width = 10, height = 8)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
cat("Saved: coupling_diagnostics.pdf\n")

# ============================================================================
# 6. SUMMARY VERDICT
# ============================================================================

cat("\n=== SUMMARY ===\n")
issues <- character(0)

if (top1/mean_score < 3) issues <- c(issues, "- Weak signal (low top/mean ratio)")
if (skewness < 1) issues <- c(issues, "- No heavy tail in distribution")
if (exists("neff") && neff/L < 1) issues <- c(issues, "- CRITICAL: Neff/L < 1")
if (exists("neff") && neff/L >= 1 && neff/L < 5) issues <- c(issues, "- Limited diversity (Neff/L < 5)")

if (length(issues) == 0) {
  cat("Data quality appears adequate for contact prediction\n")
} else {
  cat("Potential issues detected:\n")
  cat(paste(issues, collapse = "\n"), "\n")
  cat("\nTo improve: add more diverse sequences to alignment\n")
}


# Calculate Neff from alignment
# Neff = sum of sequence weights, where weight = 1/(number of sequences >= 80% identical)

library(Biostrings)

aln_file <- "results/evcouplings/rbcL/align/rbcL.a2m"  # or .sto, .fasta
aln <- readAAStringSet(aln_file)

# Convert to matrix
aln_mat <- as.matrix(aln)
n_seqs <- nrow(aln_mat)
cat("Sequences:", n_seqs, "\n")

# Pairwise identity (slow for large alignments - sample if needed)

# Calculate weights
threshold <- 0.9
weights <- numeric(n_seqs)
#initialize progress bar
pb <- txtProgressBar(min = 0, max = n_seqs, style = 3)
for (i in 1:n_seqs) {
  # Count sequences >= 80% identical to sequence i
  identities <- rowMeans(aln_mat == aln_mat[i, ], na.rm = TRUE)
  n_similar <- sum(identities >= threshold)
  weights[i] <- 1 / n_similar
  setTxtProgressBar(pb, i)
}
close(pb)


neff <- sum(weights)
L <- ncol(aln_mat)

cat("Neff:", round(neff, 1), "\n")
cat("L:", L, "\n")
cat("Neff/L:", round(neff/L, 2), "\n")

if (neff/L < 1) cat("CRITICAL: Neff/L < 1\n")
if (neff/L >= 1 & neff/L < 5) cat("WARNING: Neff/L 1-5, marginal\n")
if (neff/L >= 5) cat("GOOD: Neff/L >= 5\n")


library(data.table)

# Structure hits - predicted contacts validated against PDB
hits <- fread("results/evcouplings/rbcL/compare/rbcL_structure_hits.csv")
cat("Columns:\n")
print(names(hits))
cat("\nDimensions:", nrow(hits), "x", ncol(hits), "\n")
cat("\nFirst rows:\n")
print(head(hits))

# Key metrics to look for:
# - precision: fraction of top predicted contacts that are true contacts
# - typically reported at L, L/2, L/5 top predictions

if ("precision" %in% names(hits)) {
  cat("\nPrecision:\n")
  print(summary(hits$precision))
}

# Check the unfiltered version too
hits_unfilt <- fread("results/evcouplings/rbcL/compare/rbcL_structure_hits_unfiltered.csv")
cat("\nUnfiltered hits:", nrow(hits_unfilt), "\n")
print(head(hits_unfilt))