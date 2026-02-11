




# %% Cell 1 - Load data
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

# Load saved results
data = np.load("data/tmp/covariance/multigene.npz", allow_pickle=True)
res = {k: data[k] for k in data.files}

# Fix types
res["gene_names"] = list(res["gene_names"])
res["gene_lengths"] = list(res["gene_lengths"])
res["ids"] = list(res["ids"])
res["n"] = int(res["n"])

print(f"Loaded: {len(res['gene_names'])} genes, {res['n']} samples")
print(f"Covariance matrix: {res['cov'].shape}")

# %% Cell 2 - Quick look at matrix
plt.figure(figsize=(14, 14))
plt.imshow(res["cov"], cmap="YlOrRd", aspect="auto")
plt.colorbar(label="RV coefficient")
plt.title(f"Chloroplast multi-gene covariance ({len(res['gene_names'])} genes)")
plt.tight_layout()
#plt.show()
#save figure
plt.savefig("data/tmp/covariance/multigene_covariance_matrix.png", dpi=300)

# %% Cell 3 - Summarize within vs between
def summarize_covariance(res):
    cov = res["cov"]
    names = res["gene_names"]
    bounds = res["gene_bounds"]
    n_genes = len(names)
    
    rows = []
    for i in range(n_genes):
        for j in range(i, n_genes):
            start_i, end_i = bounds[i], bounds[i + 1]
            start_j, end_j = bounds[j], bounds[j + 1]
            
            block = cov[start_i:end_i, start_j:end_j]
            
            if i == j:
                mask = np.triu(np.ones_like(block, dtype=bool), k=1)
                vals = block[mask]
            else:
                vals = block.flatten()
            
            vals = vals[~np.isnan(vals)]
            
            rows.append({
                "gene1": names[i],
                "gene2": names[j],
                "type": "within" if i == j else "between",
                "mean_rv": np.mean(vals) if len(vals) > 0 else np.nan,
                "median_rv": np.median(vals) if len(vals) > 0 else np.nan,
                "std_rv": np.std(vals) if len(vals) > 0 else np.nan,
                "n_pairs": len(vals)
            })
    
    return pl.DataFrame(rows)

summary = summarize_covariance(res)
print(summary)

# %% Cell 4 - Within vs between comparison
within = summary.filter(pl.col("type") == "within")["mean_rv"].to_numpy()
between = summary.filter(pl.col("type") == "between")["mean_rv"].to_numpy()

print(f"Within-gene mean RV:  {np.nanmean(within):.4f} ± {np.nanstd(within):.4f}")
print(f"Between-gene mean RV: {np.nanmean(between):.4f} ± {np.nanstd(between):.4f}")

plt.figure(figsize=(8, 5))
plt.boxplot([within[~np.isnan(within)], between[~np.isnan(between)]], 
            labels=["Within gene", "Between genes"])
plt.ylabel("Mean RV coefficient")
plt.title("Within vs Between Gene Covariance")
plt.show()

# %% Cell 5 - Gene-level heatmap (mean RV per gene pair)
names = res["gene_names"]
n_genes = len(names)
bounds = res["gene_bounds"]
cov = res["cov"]

gene_cov = np.full((n_genes, n_genes), np.nan)

for i in range(n_genes):
    for j in range(n_genes):
        start_i, end_i = bounds[i], bounds[i + 1]
        start_j, end_j = bounds[j], bounds[j + 1]
        block = cov[start_i:end_i, start_j:end_j]
        vals = block.flatten()
        vals = vals[~np.isnan(vals)]
        gene_cov[i, j] = np.mean(vals) if len(vals) > 0 else np.nan

plt.figure(figsize=(12, 10))
plt.imshow(gene_cov, cmap="YlOrRd", aspect="auto")
plt.colorbar(label="Mean RV coefficient")
plt.xticks(range(n_genes), names, rotation=90, fontsize=7)
plt.yticks(range(n_genes), names, fontsize=7)
plt.title("Gene-level mean covariance")
plt.tight_layout()
#plt.show()
#save figure
plt.savefig("data/tmp/covariance/gene_covariance_heatmap.png", dpi=300)

# %% Cell 6 - Cluster genes by covariance
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

# Convert similarity to distance
gene_dist = 1 - gene_cov
np.fill_diagonal(gene_dist, 0)
gene_dist = np.nan_to_num(gene_dist, nan=1.0)
gene_dist = (gene_dist + gene_dist.T) / 2  # symmetrize

# Cluster
Z = linkage(squareform(gene_dist), method="ward")

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Dendrogram
dendrogram(Z, labels=names, ax=axes[0], leaf_rotation=90)
axes[0].set_title("Gene clustering by covariance")

# Reordered heatmap
order = leaves_list(Z)
gene_cov_ordered = gene_cov[np.ix_(order, order)]
names_ordered = [names[i] for i in order]

im = axes[1].imshow(gene_cov_ordered, cmap="YlOrRd", aspect="auto")
plt.colorbar(im, ax=axes[1], label="Mean RV")
axes[1].set_xticks(range(n_genes))
axes[1].set_xticklabels(names_ordered, rotation=90, fontsize=7)
axes[1].set_yticks(range(n_genes))
axes[1].set_yticklabels(names_ordered, fontsize=7)
axes[1].set_title("Clustered gene covariance")

plt.tight_layout()
plt.show()
#save figure   
fig.savefig("data/tmp/covariance/gene_covariance_clustering.png", dpi=300)

# %% Cell - Sanity checks

import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.load("data/tmp/covariance/multigene.npz", allow_pickle=True)
res = {k: data[k] for k in data.files}
res["gene_names"] = list(res["gene_names"])
res["gene_lengths"] = list(res["gene_lengths"])

cov = res["cov"]
names = res["gene_names"]
bounds = res["gene_bounds"]

# %% Check 1: Basic stats
print("=== BASIC STATS ===")
print(f"Matrix shape: {cov.shape}")
print(f"NaN count: {np.isnan(cov).sum()} ({100*np.isnan(cov).sum()/cov.size:.1f}%)")
print(f"Min: {np.nanmin(cov):.4f}")
print(f"Max: {np.nanmax(cov):.4f}")
print(f"Mean: {np.nanmean(cov):.4f}")
print(f"Median: {np.nanmedian(cov):.4f}")

# %% Check 2: Distribution of values
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
vals = cov.flatten()
vals = vals[~np.isnan(vals)]
plt.hist(vals, bins=100)
plt.xlabel("RV coefficient")
plt.ylabel("Count")
plt.title(f"Distribution (n={len(vals):,})")

plt.subplot(1, 2, 2)
plt.hist(vals, bins=100, log=True)
plt.xlabel("RV coefficient")
plt.ylabel("Count (log)")
plt.title("Log scale")
plt.tight_layout()
plt.show()

# %% Check 3: How many values are exactly 1 or very close?
print("\n=== VALUES NEAR 1 ===")
for thresh in [1.0, 0.9999, 0.999, 0.99, 0.95, 0.9]:
    n = np.sum(vals >= thresh)
    print(f"RV >= {thresh}: {n:,} ({100*n/len(vals):.2f}%)")

# %% Check 4: Diagonal check (should be 1.0)
diag = np.diag(cov)
print("\n=== DIAGONAL ===")
print(f"Diagonal mean: {np.nanmean(diag):.6f}")
print(f"Diagonal min: {np.nanmin(diag):.6f}")
print(f"Diagonal max: {np.nanmax(diag):.6f}")
print(f"All diagonal == 1? {np.allclose(diag[~np.isnan(diag)], 1.0)}")

# %% Check 5: Off-diagonal distribution
off_diag = cov[~np.eye(cov.shape[0], dtype=bool)]
off_diag = off_diag[~np.isnan(off_diag)]
print("\n=== OFF-DIAGONAL ===")
print(f"Mean: {np.mean(off_diag):.4f}")
print(f"Median: {np.median(off_diag):.4f}")
print(f"Std: {np.std(off_diag):.4f}")

# %% Check 6: Within-gene vs between-gene (excluding diagonal)
print("\n=== WITHIN VS BETWEEN (excluding diagonal) ===")
n_genes = len(names)

within_vals = []
between_vals = []

for i in range(n_genes):
    for j in range(n_genes):
        start_i, end_i = bounds[i], bounds[i + 1]
        start_j, end_j = bounds[j], bounds[j + 1]
        block = cov[start_i:end_i, start_j:end_j]
        
        if i == j:
            # Within: upper triangle excluding diagonal
            mask = np.triu(np.ones_like(block, dtype=bool), k=1)
            within_vals.extend(block[mask].tolist())
        elif i < j:
            # Between: all values
            between_vals.extend(block.flatten().tolist())

within_vals = np.array(within_vals)
between_vals = np.array(between_vals)
within_vals = within_vals[~np.isnan(within_vals)]
between_vals = between_vals[~np.isnan(between_vals)]

print(f"Within-gene:  mean={np.mean(within_vals):.4f}, median={np.median(within_vals):.4f}, n={len(within_vals):,}")
print(f"Between-gene: mean={np.mean(between_vals):.4f}, median={np.median(between_vals):.4f}, n={len(between_vals):,}")

plt.figure(figsize=(8, 5))
plt.hist(within_vals, bins=50, alpha=0.5, label=f"Within (n={len(within_vals):,})", density=True)
plt.hist(between_vals, bins=50, alpha=0.5, label=f"Between (n={len(between_vals):,})", density=True)
plt.xlabel("RV coefficient")
plt.ylabel("Density")
plt.legend()
plt.title("Within vs Between gene RV distribution")
plt.show()

# %% Check 7: Look at a specific gene pair - are values sensible?
print("\n=== SPOT CHECK: rbcL vs psbA ===")
if "rbcL" in names and "psbA" in names:
    i = names.index("rbcL")
    j = names.index("psbA")
    start_i, end_i = bounds[i], bounds[i + 1]
    start_j, end_j = bounds[j], bounds[j + 1]
    block = cov[start_i:end_i, start_j:end_j]
    print(f"Block shape: {block.shape}")
    print(f"Block mean: {np.nanmean(block):.4f}")
    print(f"Block min: {np.nanmin(block):.4f}")
    print(f"Block max: {np.nanmax(block):.4f}")
    print(f"Block % NaN: {100*np.isnan(block).sum()/block.size:.1f}%")
    
    plt.figure(figsize=(8, 6))
    plt.imshow(block, cmap="YlOrRd", aspect="auto")
    plt.colorbar(label="RV")
    plt.xlabel("psbA position")
    plt.ylabel("rbcL position")
    plt.title("rbcL vs psbA covariance block")
    plt.show()

# %% Check 8: Are monomorphic positions driving high RV?
# RV of two constant vectors is undefined/1.0
print("\n=== CHECK FOR LOW VARIANCE POSITIONS ===")

# Reload kidera to check variance
from Bio import SeqIO
import polars as pl

def load_kidera(path="data/kidera_factors.csv"):
    df = pl.read_csv(path)
    aa_order = df.get_column(df.columns[0]).to_list()
    arr = df.select(df.columns[1:]).to_numpy()
    return arr, {aa: i for i, aa in enumerate(aa_order)}

KIDERA_ARR, AA_TO_IDX = load_kidera()

def encode_kidera(seqs):
    n_samples, n_pos = seqs.shape
    out = np.full((n_samples, n_pos, 10), np.nan)
    for aa, idx in AA_TO_IDX.items():
        mask = seqs == aa
        out[mask] = KIDERA_ARR[idx]
    return out

# Check one gene
gene = "rbcL"
fasta = f"data/tmp/alignedGenes/{gene}_AA_aligned.fasta"
records = list(SeqIO.parse(fasta, "fasta"))
seqs = np.array([list(str(r.seq)) for r in records])

# Filter gaps
gap_frac = np.mean(seqs == "-", axis=0)
keep = gap_frac <= 0.5
seqs = seqs[:, keep]

kidera = encode_kidera(seqs)

# Variance per position (across samples, summed over 10 KFs)
pos_var = np.nanvar(kidera, axis=0).sum(axis=1)
print(f"Position variance stats for {gene}:")
print(f"  Min: {pos_var.min():.6f}")
print(f"  Max: {pos_var.max():.6f}")
print(f"  Mean: {pos_var.mean():.6f}")
print(f"  Positions with var < 0.01: {np.sum(pos_var < 0.01)}")
print(f"  Positions with var < 0.1: {np.sum(pos_var < 0.1)}")

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(pos_var, bins=50)
plt.xlabel("Position variance (sum of 10 KFs)")
plt.ylabel("Count")
plt.title(f"{gene}: Position variance distribution")

plt.subplot(1, 2, 2)
plt.plot(pos_var)
plt.xlabel("Position")
plt.ylabel("Variance")
plt.title(f"{gene}: Variance along sequence")
plt.tight_layout()
plt.show()

# %% Check 9: RV formula sanity - test on random data
print("\n=== RV FORMULA SANITY ===")

def rv_coefficient(X, Y):
    X = X - np.nanmean(X, axis=0)
    Y = Y - np.nanmean(Y, axis=0)
    Sxy = X.T @ Y
    Sxx = X.T @ X
    Syy = Y.T @ Y
    return np.sum(Sxy**2) / np.sqrt(np.sum(Sxx**2) * np.sum(Syy**2))

# Test 1: Identical data should give RV = 1
X = np.random.randn(100, 10)
print(f"RV(X, X) = {rv_coefficient(X, X):.6f} (should be 1.0)")

# Test 2: Independent random data should give RV << 1
Y = np.random.randn(100, 10)
print(f"RV(X, random Y) = {rv_coefficient(X, Y):.6f} (should be ~0)")

# Test 3: Linearly related data
Y = X @ np.random.randn(10, 10)  # linear transform
print(f"RV(X, linear transform of X) = {rv_coefficient(X, Y):.6f} (should be ~1)")

# Test 4: Low variance data
X_low = np.random.randn(100, 10) * 0.001
Y_low = np.random.randn(100, 10) * 0.001
print(f"RV(low var X, low var Y) = {rv_coefficient(X_low, Y_low):.6f}")

# Test 5: Nearly constant data
X_const = np.ones((100, 10)) + np.random.randn(100, 10) * 1e-10
Y_const = np.ones((100, 10)) + np.random.randn(100, 10) * 1e-10
print(f"RV(near-constant X, near-constant Y) = {rv_coefficient(X_const, Y_const):.6f}")