import pyarrow.parquet as pq
import pandas as pd
import numpy as np
from pathlib import Path
import pickle

# ---- 1. LOAD PCA MODEL ----
print("Loading PCA model...")
pca_file = Path("data/tmp/pca_model_standard_sample100.pkl")

with open(pca_file, 'rb') as f:
    pca = pickle.load(f)

print(f"PCA model loaded: {pca.n_components_} components")

# ---- 2. LOAD TREE TAXA ----
print("\nLoading tree taxa...")
tree_taxa_file = Path("data/tree_tip_labels.txt")

with open(tree_taxa_file, 'r') as f:
    tree_taxa = set(line.strip() for line in f)
print(f"Loaded {len(tree_taxa)} taxa from tree")

# ---- 3. PROJECT ALL EMBEDDINGS ----
embed_dir = Path("data/embeddings/")
parquet_files = sorted(list(embed_dir.glob("*_residue_embeddings.parquet")))
assert len(parquet_files) > 0, "No embedding files found"

print(f"\nFound {len(parquet_files)} files to project\n")

# Collect all projected data in memory
all_projected = []

for idx, pf in enumerate(parquet_files, 1):
    print(f"[{idx}/{len(parquet_files)}] Projecting: {pf.name}")
    
    # Load file
    table = pq.read_table(pf)
    df = table.to_pandas()
    
    # Filter to tree taxa
    n_before = len(df)
    df = df[df['ID'].isin(tree_taxa)]
    n_after = len(df)
    
    if n_after == 0:
        print(f"  Skipping - no matching taxa\n")
        continue
    
    print(f"  Kept {n_after:,} / {n_before:,} residues")
    
    # Extract embeddings and metadata
    embed_cols = [c for c in df.columns if c.startswith('embedding_')]
    metadata_cols = [c for c in df.columns if not c.startswith('embedding_')]
    
    X = df[embed_cols].values
    
    # Project to PC space (this is prcomp()$x in R)
    pc_scores = pca.transform(X)
    
    print(f"  Projected to {pc_scores.shape[1]} PCs")
    print(f"  PC1 range: [{pc_scores[:, 0].min():.2f}, {pc_scores[:, 0].max():.2f}]")
    
    # Create PC dataframe
    pc_df = pd.DataFrame(
        pc_scores,
        columns=[f'PC{i+1}' for i in range(pc_scores.shape[1])]
    )
    
    # Combine with metadata
    result_df = pd.concat([df[metadata_cols].reset_index(drop=True), pc_df], axis=1)
    
    all_projected.append(result_df)
    print()

# ---- 4. COMBINE AND SAVE ----
print("Combining all projections...")
combined_df = pd.concat(all_projected, ignore_index=True)
print(f"Combined dataset: {len(combined_df):,} rows")

output_file = Path("data/tmp/all_embeddings_projected.parquet")
combined_df.to_parquet(output_file, index=False)
print(f"\nSaved to: {output_file}")

# ---- 5. SUMMARY STATS ----
print(f"\n{'='*60}")
print("PROJECTION SUMMARY:")
print(f"  Total residues: {len(combined_df):,}")
print(f"  Genes: {combined_df['Gene'].nunique()}")
print(f"  PC1 range: [{combined_df['PC1'].min():.2f}, {combined_df['PC1'].max():.2f}]")
print(f"  PC2 range: [{combined_df['PC2'].min():.2f}, {combined_df['PC2'].max():.2f}]")
print(f"{'='*60}")

print("\nColumns in output:")
print(combined_df.columns.tolist()[:10], "...")