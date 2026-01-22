"""
Create genome structure matrix: rows = genomes, cols = gene positions (relative)
"""

import pandas as pd
import numpy as np

df = pd.read_csv("data/gbf_gene_coords.csv")

# Normalize positions to relative (0-1) scale within each genome
df["rel_start"] = df["start"] / df["genome_size"]
df["rel_mid"] = (df["start"] + df["end"]) / 2 / df["genome_size"]

# Pivot: rows = genomes, cols = genes, values = relative midpoint
structure_mat = df.pivot_table(
    index="genome_id", 
    columns="gene", 
    values="rel_mid",
    aggfunc="first"  # some genomes may have duplicates (IR regions)
)

print(f"Shape: {structure_mat.shape}")
print(f"Genomes: {structure_mat.shape[0]}")
print(f"Genes: {structure_mat.shape[1]}")

# Missingness summary
missing = structure_mat.isna().sum()
print(f"\nGenes with most missing: \n{missing.sort_values(ascending=False).head(10)}")
print(f"\nGenes present in all: {(missing == 0).sum()}")

# Also get strand info (could indicate inversions)
strand_mat = df.pivot_table(
    index="genome_id",
    columns="gene",
    values="strand",
    aggfunc="first"
)

# Save both
structure_mat.to_csv("data/genome_structure_positions.csv")
strand_mat.to_csv("data/genome_structure_strands.csv")

# Quick check: how many unique gene orders?
# Use rank within genome as proxy for order
df["rank"] = df.groupby("genome_id")["start"].rank()
order_mat = df.pivot_table(
    index="genome_id",
    columns="gene",
    values="rank",
    aggfunc="first"
)
order_mat.to_csv("data/genome_structure_order.csv")

# Count unique order patterns (for genes present in >90% of genomes)
common_genes = missing[missing < 0.1 * len(structure_mat)].index.tolist()
order_common = order_mat[common_genes].dropna()
unique_orders = order_common.drop_duplicates()
print(f"\nUnique gene orders (among {len(common_genes)} common genes): {len(unique_orders)}")