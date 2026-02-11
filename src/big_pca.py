import pyarrow.parquet as pq
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from pathlib import Path
import pickle

# ---- 1. LOAD TREE TAXA ----
print("Loading tree taxa...")
tree_taxa_file = Path("data/tmp/tree_tip_labels.txt")

if tree_taxa_file.exists():
    with open(tree_taxa_file, 'r') as f:
        tree_taxa = set(line.strip() for line in f)
    print(f"Loaded {len(tree_taxa)} taxa from tree")
else:
    print("ERROR: Create tree_tip_labels.txt first")
    print("In R: writeLines(tree$tip.label, 'data/tmp/tree_tip_labels.txt')")
    exit(1)

# ---- 2. SETUP ----
embed_dir = Path("data/embeddings/")
parquet_files = sorted(list(embed_dir.glob("*_residue_embeddings.parquet")))
assert len(parquet_files) > 0, "No embedding files found"
print(f"\nFound {len(parquet_files)} embedding files\n")

output_dir = Path("data/tmp/")
output_dir.mkdir(parents=True, exist_ok=True)
"""
# ---- 3. RANDOM SAMPLING STRATEGY ----
# First pass: count total rows after filtering
print("First pass: counting rows...\n")
total_rows = 0
file_row_counts = []

for idx, pf in enumerate(parquet_files, 1):
    print(f"[{idx}/{len(parquet_files)}] Counting: {pf.name}")
    table = pq.read_table(pf, columns=['ID'])
    df = table.to_pandas()
    
    n_before = len(df)
    n_after = df['ID'].isin(tree_taxa).sum()
    
    print(f"  Rows: {n_after:,} (filtered from {n_before:,})")
    
    file_row_counts.append({
        'file': pf,
        'n_rows': n_after
    })
    total_rows += n_after

print(f"\n{'='*60}")
print(f"Total rows across all files: {total_rows:,}")
print(f"{'='*60}\n")

# ---- 4. DECIDE SAMPLING ----
# Option 1: Load everything (if < 50M rows)
# Option 2: Sample uniformly (if 50M-200M rows)
# Option 3: Sample heavily (if > 200M rows)

LOAD_ALL_THRESHOLD = 50_000_000
SAMPLE_THRESHOLD = 200_000_000

if total_rows < LOAD_ALL_THRESHOLD:
    print(f"Dataset is manageable ({total_rows:,} rows). Loading all data.\n")
    sample_frac = 1.0
    n_components = 100
elif total_rows < SAMPLE_THRESHOLD:
    sample_frac = LOAD_ALL_THRESHOLD / total_rows
    print(f"Dataset is large ({total_rows:,} rows). Sampling {sample_frac:.2%} of data.\n")
    n_components = 100
else:
    sample_frac = LOAD_ALL_THRESHOLD / total_rows
    print(f"Dataset is very large ({total_rows:,} rows). Sampling {sample_frac:.2%} of data.\n")
    n_components = 100

print(f"Expected rows to load: {int(total_rows * sample_frac):,}")
print(f"Expected memory (embeddings only): ~{int(total_rows * sample_frac * 960 * 4 / 1e9):.1f} GB")
print(f"Number of PCs: {n_components}\n")

# ---- 5. LOAD DATA WITH RANDOM SAMPLING ----
print("Loading and sampling data...\n")

np.random.seed(42)
all_data = []
load_stats = []

for idx, file_info in enumerate(file_row_counts, 1):
    pf = file_info['file']
    n_rows = file_info['n_rows']
    
    if n_rows == 0:
        continue
    
    print(f"[{idx}/{len(parquet_files)}] Loading: {pf.name}")
    
    table = pq.read_table(pf)
    df = table.to_pandas()
    
    # Filter to tree taxa
    df_filtered = df[df['ID'].isin(tree_taxa)].copy()
    
    # Random sample
    if sample_frac < 1.0:
        n_sample = max(1, int(len(df_filtered) * sample_frac))
        df_sampled = df_filtered.sample(n=n_sample, random_state=42 + idx)
        print(f"  Sampled {n_sample:,} / {len(df_filtered):,} rows")
    else:
        df_sampled = df_filtered
        print(f"  Loaded {len(df_sampled):,} rows")
    
    all_data.append(df_sampled)
    
    load_stats.append({
        'gene': pf.stem.replace('_residue_embeddings', ''),
        'n_total': len(df_filtered),
        'n_sampled': len(df_sampled)
    })

print(f"\nConcatenating all data...")
combined_df = pd.concat(all_data, ignore_index=True)
print(f"Combined dataset: {len(combined_df):,} rows")

# Save load stats
stats_df = pd.DataFrame(load_stats)
stats_file = output_dir / "pca_load_stats.csv"
stats_df.to_csv(stats_file, index=False)
print(f"Saved load stats to {stats_file}")

# ---- 6. EXTRACT EMBEDDINGS ----
print("\nExtracting embedding columns...")
embed_cols = [c for c in combined_df.columns if c.startswith('embedding_')]
metadata_cols = [c for c in combined_df.columns if not c.startswith('embedding_')]

X = combined_df[embed_cols].values
print(f"Embedding matrix: {X.shape}")
assert X.shape[1] == 960, f"Expected 960 dimensions, got {X.shape[1]}"

# ---- 7. RUN STANDARD PCA ----
pca_file = output_dir / f"pca_model_standard_sample{int(sample_frac*100)}.pkl"

print(f"\nRunning standard PCA with {n_components} components...")
print("This may take a few minutes...\n")

pca = PCA(n_components=n_components, random_state=42)
pca.fit(X)

print(f"PCA complete!")

# Save model
with open(pca_file, 'wb') as f:
    pickle.dump(pca, f)
print(f"Saved PCA model to {pca_file}")

# ---- 8. VARIANCE EXPLAINED ----
var_explained = pca.explained_variance_ratio_
cumvar = np.cumsum(var_explained)

print(f"\n{'='*60}")
print("VARIANCE EXPLAINED:")
print(f"  PC1:     {var_explained[0]*100:.2f}%")
print(f"  PC1-10:  {var_explained[:10].sum()*100:.2f}%")
print(f"  PC1-20:  {var_explained[:20].sum()*100:.2f}%")
print(f"  PC1-50:  {var_explained[:50].sum()*100:.2f}%")
print(f"  PC1-100: {var_explained[:100].sum()*100:.2f}%")
print(f"{'='*60}\n")

# Save variance explained
var_df = pd.DataFrame({
    'PC': range(1, n_components + 1),
    'variance_explained': var_explained,
    'cumulative_variance': cumvar
})
var_file = output_dir / "pca_variance_explained.csv"
var_df.to_csv(var_file, index=False)
print(f"Saved variance explained to {var_file}")

# ---- 9. TRANSFORM SAMPLED DATA ----
print("\nTransforming sampled data...")
#load the pca from the pickle"""
pca_file = "results/pca_model_standard_sample100.pkl"
pca = pickle.load(open(pca_file, 'rb'))
pc_scores = pca.transform(X)

pc_df = pd.DataFrame(
    pc_scores,
    columns=[f'PC{i+1}' for i in range(n_components)]
)

result_df = pd.concat([combined_df[metadata_cols].reset_index(drop=True), pc_df], axis=1)

sampled_file = output_dir / "pca_scores_sampled_combined.parquet"
result_df.to_parquet(sampled_file, index=False)
print(f"Saved combined sampled PC scores to {sampled_file}")

# ---- 10. TRANSFORM ALL FILES ----
print("\nTransforming all files using fitted PCA...\n")

transform_stats = []

for idx, pf in enumerate(parquet_files, 1):
    print(f"[{idx}/{len(parquet_files)}] Transforming: {pf.name}")
    table = pq.read_table(pf)
    df = table.to_pandas()
    
    # Filter to tree taxa
    n_before = len(df)
    df = df[df['ID'].isin(tree_taxa)]
    n_after = len(df)
    
    if n_after == 0:
        print(f"  Skipping - no matching taxa")
        print()
        continue
    
    print(f"  Kept {n_after:,} / {n_before:,} residues")
    
    embed_cols = [c for c in df.columns if c.startswith('embedding_')]
    metadata_cols = [c for c in df.columns if not c.startswith('embedding_')]
    
    X = df[embed_cols].values
    print(f"  Transforming {X.shape[0]:,} x {X.shape[1]} -> {X.shape[0]:,} x {n_components}")
    
    pc_scores = pca.transform(X)
    
    # Per-file stats
    pc_means = pc_scores.mean(axis=0)
    pc_stds = pc_scores.std(axis=0)
    print(f"  PC1 range: [{pc_scores[:, 0].min():.2f}, {pc_scores[:, 0].max():.2f}]")
    print(f"  PC1 mean ± std: {pc_means[0]:.2f} ± {pc_stds[0]:.2f}")
    
    # Create PC dataframe
    pc_df = pd.DataFrame(
        pc_scores,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    
    # Combine metadata with PCs
    result_df = pd.concat([df[metadata_cols].reset_index(drop=True), pc_df], axis=1)
    
    # Save as parquet
    gene_name = pf.stem.replace('_residue_embeddings', '')
    out_file = output_dir / f"{gene_name}_pca_scores.parquet"
    result_df.to_parquet(out_file, index=False)
    print(f"  Saved: {out_file.name}")
    
    transform_stats.append({
        'gene': gene_name,
        'n_residues': n_after,
        'pc1_mean': pc_means[0],
        'pc1_std': pc_stds[0],
        'pc1_min': pc_scores[:, 0].min(),
        'pc1_max': pc_scores[:, 0].max()
    })
    print()

# Save transform stats
trans_stats_df = pd.DataFrame(transform_stats)
trans_stats_file = output_dir / "pca_transform_stats.csv"
trans_stats_df.to_csv(trans_stats_file, index=False)
print(f"Saved transform stats to {trans_stats_file}")

print(f"\n{'='*60}")
print("DONE!")
print(f"PCA fitted on {len(combined_df):,} randomly sampled residues")
print(f"All {len(parquet_files)} genes transformed")
print(f"PC scores saved to data/tmp/")
print(f"{'='*60}")
