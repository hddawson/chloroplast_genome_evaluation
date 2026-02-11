import pyarrow.parquet as pq
import pandas as pd
import numpy as np
from sklearn.decomposition import IncrementalPCA
from pathlib import Path
import pickle

# ---- 1. LOAD TREE TAXA ----
print("Loading tree taxa...")
# Read tree tip labels from R object or file
# Assuming you have a simple text file with one ID per line
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

# ---- 3. INCREMENTAL PCA FIT ----
n_components = 100
batch_size = 1_000_000

pca_file = output_dir / "ipca_model_5m_batchsize.pkl"

if not pca_file.exists():
    print(f"Fitting IncrementalPCA with {n_components} components...")
    print(f"Batch size: {batch_size}\n")
    ipca = IncrementalPCA(n_components=n_components, batch_size=batch_size)
    
    total_rows = 0
    total_filtered = 0
    gene_stats = []
    
    for idx, pf in enumerate(parquet_files, 1):
        print(f"[{idx}/{len(parquet_files)}] Fitting: {pf.name}")
        table = pq.read_table(pf)
        df = table.to_pandas()
        
        # Filter to tree taxa
        n_before = len(df)
        df_filtered = df[df['ID'].isin(tree_taxa)]
        n_after = len(df_filtered)
        n_removed = n_before - n_after
        
        print(f"  Rows before filter: {n_before:,}")
        print(f"  Rows after filter: {n_after:,}")
        print(f"  Removed: {n_removed:,} ({n_removed/n_before*100:.1f}%)")
        
        if n_after == 0:
            print(f"  WARNING: No taxa match tree for {pf.name}")
            continue
        
        embed_cols = [c for c in df_filtered.columns if c.startswith('embedding_')]
        X = df_filtered[embed_cols].values
        
        assert X.shape[1] == 960, f"Expected 960 embeddings, got {X.shape[1]}"
        
        # Fit in batches
        batches_processed = 0
        for i in range(0, len(X), batch_size):
            batch = X[i:i+batch_size]
            ipca.partial_fit(batch)
            batches_processed += 1
        
        print(f"  Processed {batches_processed} batches")
        
        total_rows += n_after
        total_filtered += n_removed
        
        gene_stats.append({
            'gene': pf.stem.replace('_residue_embeddings', ''),
            'n_residues': n_after,
            'n_filtered': n_removed
        })
        print()
    
    print(f"=" * 60)
    print(f"Total residues used for PCA: {total_rows:,}")
    print(f"Total residues filtered out: {total_filtered:,}")
    print(f"=" * 60)
    
    # Save model
    with open(pca_file, 'wb') as f:
        pickle.dump(ipca, f)
    print(f"\nSaved PCA model to {pca_file}")
    
    # Save gene stats
    stats_df = pd.DataFrame(gene_stats)
    stats_file = output_dir / "pca_gene_stats.csv"
    stats_df.to_csv(stats_file, index=False)
    print(f"Saved gene stats to {stats_file}")
else:
    print("Loading existing PCA model...")
    with open(pca_file, 'rb') as f:
        ipca = pickle.load(f)

# Variance explained
var_explained = ipca.explained_variance_ratio_
cumvar = np.cumsum(var_explained)

print(f"\n{'='*60}")
print("VARIANCE EXPLAINED:")
print(f"  PC1-10:  {var_explained[:10].sum()*100:.2f}%")
print(f"  PC1-20:  {var_explained[:20].sum()*100:.2f}%")
print(f"  PC1-50:  {var_explained[:50].sum()*100:.2f}%")
print(f"  PC1-100: {var_explained[:100].sum()*100:.2f}%")
print(f"{'='*60}\n")

# ---- 4. TRANSFORM AND SAVE PC SCORES ----
print("Transforming data and saving PC scores...\n")

transform_stats = []

for idx, pf in enumerate(parquet_files, 1):
    print(f"[{idx}/{len(parquet_files)}] Processing: {pf.name}")
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
    
    pc_scores = ipca.transform(X)
    
    # Per-file PCA stats
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
    result_df = pd.concat([df[metadata_cols].reset_index(drop=True), 
                           pc_df], axis=1)
    
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
print("DONE! PC scores saved to data/tmp/")
print(f"Total genes processed: {len(transform_stats)}")
print(f"{'='*60}")
