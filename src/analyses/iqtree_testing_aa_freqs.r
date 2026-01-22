# ---- iqtree eda ----
"""
h0

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -s supergene.faa -p partition_filtered.nex -m LG+FO+G -ft fast.tree -n 0 -nt 20 -pre H0_fixed


h1 

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -s supergene.faa -p partition_filtered.nex -m LG+FO+G -ft fast.tree -n 0 -nt 20 -pre H1_fixed
"""
# Extract log-likelihoods
lnL_H0 <- -4727634.696  # Shared frequencies

lnL_H1 <- -4727612.807  # Partition-specific frequencies

# LRT statistic
LRT <- 2 * (lnL_H1 - lnL_H0)
cat("LRT statistic:", LRT, "\n")

# Degrees of freedom
# How many partitions do you have?
n_partitions <- 2  # UPDATE THIS with your actual number

df <- 19 * (n_partitions - 1)  # 19 free AA frequency parameters per partition
cat("Degrees of freedom:", df, "\n")

# P-value
p_value <- pchisq(LRT, df, lower.tail = FALSE)
cat("P-value:", p_value, "\n")

if (p_value < 0.05) {
  cat("\n✓ SIGNIFICANT: AA frequencies differ between site classes (P <", p_value, ")\n")
} else {
  cat("\n✗ NOT SIGNIFICANT: No evidence for different AA frequencies (P =", p_value, ")\n")
}
