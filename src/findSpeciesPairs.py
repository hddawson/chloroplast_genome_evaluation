#!/usr/bin/env python3
"""
Find pairs of samples that maximize differences at significant association sites.
Identifies optimal accession pairs for reciprocal crosses based on chloroplast variants
at statistically important residues.
"""

import pandas as pd
import numpy as np
from itertools import combinations
import sys

# Configuration
VARIANT_FILE = 'data/chloroplast_mitochondria_variant_calls.txt'
RESULTS_FILE = 'results/reference_mapped_results.csv'
VARIANT_ANALYSIS_FILE = 'results/variant_analysis.csv'
OUTPUT_FILE = 'results/optimal_cross_pairs.csv'

# Load data
print("Loading data...")
variant_df = pd.read_csv(VARIANT_FILE, sep='\t')
c_variants = variant_df[variant_df['CHROM'] == 'chloroplast'].copy()
results_df = pd.read_csv(RESULTS_FILE)
variant_analysis = pd.read_csv(VARIANT_ANALYSIS_FILE)

# Get significant sites (lower quartile P_res)
p_res_threshold = results_df['P_res'].quantile(0.1)
significant_sites = variant_analysis[
    (variant_analysis['in_CDS']) & 
    (variant_analysis['P_res'] <= p_res_threshold)
]['POS'].unique()

print(f"\nFound {len(significant_sites)} significant sites in CDS")
print(f"P_res threshold (lower quartile): {p_res_threshold:.2e}")

# Filter variants to significant sites only
sig_variants = c_variants[c_variants['POS'].isin(significant_sites)].copy()
print(f"Variants at significant sites: {len(sig_variants)}")

stopifnot = lambda cond, msg: None if cond else sys.exit(f"Assertion failed: {msg}")
stopifnot(len(sig_variants) > 0, "No variants at significant sites")

# Get sample columns
sample_cols = [col for col in sig_variants.columns if col not in ['CHROM', 'POS']]
print(f"Analyzing {len(sample_cols)} samples")

# Build genotype matrix for significant sites
# Rows = sites, Columns = samples
genotype_matrix = sig_variants[sample_cols].values

# Calculate pairwise differences
print("\nCalculating pairwise differences...")
n_samples = len(sample_cols)
diff_matrix = np.zeros((n_samples, n_samples))

for i in range(n_samples):
    for j in range(i + 1, n_samples):
        # Count positions where samples differ
        diffs = np.sum(genotype_matrix[:, i] != genotype_matrix[:, j])
        diff_matrix[i, j] = diffs
        diff_matrix[j, i] = diffs

# Find top pairs
print("\nFinding optimal pairs...")
pair_results = []

for i, j in combinations(range(n_samples), 2):
    n_diffs = int(diff_matrix[i, j])

    pair_results.append({
        'sample1': sample_cols[i],
        'sample2': sample_cols[j],
        'n_differences': n_diffs,
        'n_sites_tested': len(sig_variants),
        'percent_different': 100 * n_diffs / len(sig_variants)
    })

pair_df = pd.DataFrame(pair_results)
pair_df = pair_df.sort_values('n_differences', ascending=False)

# Save results
pair_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nSaved pairwise comparison to: {OUTPUT_FILE}")

# Summary
print("\n=== TOP 20 MOST DIVERGENT PAIRS ===")
print(pair_df.head(20).to_string(index=False))

print("\n=== SUMMARY STATISTICS ===")
print(f"Mean differences per pair: {pair_df['n_differences'].mean():.1f}")
print(f"Median differences per pair: {pair_df['n_differences'].median():.1f}")
print(f"Max differences: {pair_df['n_differences'].max()}")
print(f"Min differences: {pair_df['n_differences'].min()}")

# Find samples that appear most frequently in top pairs
print("\n=== SAMPLES IN TOP 50 MOST DIVERGENT PAIRS ===")
top_50 = pair_df.head(50)
sample_counts = pd.concat([
    top_50['sample1'].value_counts(),
    top_50['sample2'].value_counts()
]).groupby(level=0).sum().sort_values(ascending=False)

print(sample_counts.head(20))

# Additional analysis: identify samples with most unique haplotypes
print("\n=== HAPLOTYPE DIVERSITY ===")
unique_haplotypes = []
for sample in sample_cols:
    haplotype = tuple(sig_variants[sample].values)
    unique_haplotypes.append(haplotype)

from collections import Counter
haplotype_counts = Counter(unique_haplotypes)
print(f"Total unique haplotypes: {len(haplotype_counts)}")
print(f"Samples sharing most common haplotype: {max(haplotype_counts.values())}")