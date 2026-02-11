import polars as pl
import numpy as np
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
from multiprocessing import shared_memory, Pool
import os

# ---------------------------------------------------------------------
# KIDERA FACTORS
# ---------------------------------------------------------------------

def load_kidera(path: str = "data/kidera_factors.csv") -> tuple[np.ndarray, dict]:
    df = pl.read_csv(path)
    aa_order = df.get_column(df.columns[0]).to_list()
    arr = df.select(df.columns[1:]).to_numpy().astype(np.float64)
    aa_to_idx = {aa: i for i, aa in enumerate(aa_order)}
    return arr, aa_to_idx

KIDERA_ARR, AA_TO_IDX = load_kidera()

# ---------------------------------------------------------------------
# READ ALIGNMENT
# ---------------------------------------------------------------------

def read_alignment(fasta_path: str, max_gap_frac: float = 0.5) -> tuple[np.ndarray, list[str], np.ndarray]:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    assert len(records) > 0
    
    ids = [r.id.split("|")[0] for r in records]
    seqs = np.array([list(str(r.seq)) for r in records])
    
    gap_frac = np.mean(seqs == "-", axis=0)
    keep = gap_frac <= max_gap_frac
    assert keep.sum() > 1
    
    return seqs[:, keep], ids, np.where(keep)[0]

# ---------------------------------------------------------------------
# ENCODE TO KIDERA
# ---------------------------------------------------------------------

def encode_kidera(seqs: np.ndarray) -> np.ndarray:
    n_samples, n_pos = seqs.shape
    out = np.full((n_samples, n_pos, 10), np.nan, dtype=np.float64)
    
    for aa, idx in AA_TO_IDX.items():
        mask = seqs == aa
        out[mask] = KIDERA_ARR[idx]
    
    return out

# ---------------------------------------------------------------------
# VARIANCE FILTER
# ---------------------------------------------------------------------

def filter_low_variance_positions(kidera: np.ndarray, min_var: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    pos_var = np.nanvar(kidera, axis=0).sum(axis=1)
    keep = pos_var >= min_var
    return kidera[:, keep, :], np.where(keep)[0]

# ---------------------------------------------------------------------
# RV COEFFICIENT
# ---------------------------------------------------------------------

def rv_coefficient(X: np.ndarray, Y: np.ndarray) -> float:
    X = X - np.nanmean(X, axis=0)
    Y = Y - np.nanmean(Y, axis=0)
    
    Sxy = X.T @ Y
    Sxx = X.T @ X
    Syy = Y.T @ Y
    
    denom = np.sqrt(np.sum(Sxx**2) * np.sum(Syy**2))
    return np.sum(Sxy**2) / denom if denom > 0 else 0.0

# ---------------------------------------------------------------------
# SHARED MEMORY PARALLEL
# ---------------------------------------------------------------------

_shm = None
_shape = None
_dtype = None

def _init_worker(shm_name, shape, dtype):
    global _shm, _shape, _dtype
    _shm = shared_memory.SharedMemory(name=shm_name)
    _shape = shape
    _dtype = dtype

def _compute_row(i, min_obs=10):
    global _shm, _shape, _dtype
    kidera = np.ndarray(_shape, dtype=_dtype, buffer=_shm.buf)
    
    n_pos = kidera.shape[1]
    Xi = kidera[:, i, :]
    valid_i = ~np.isnan(Xi).any(axis=1)
    
    results = []
    for j in range(i, n_pos):
        Xj = kidera[:, j, :]
        valid = valid_i & ~np.isnan(Xj).any(axis=1)
        
        if valid.sum() < min_obs:
            results.append((i, j, np.nan))
        else:
            rv = rv_coefficient(Xi[valid], Xj[valid])
            results.append((i, j, rv))
    
    return results

def calc_position_covariance(kidera: np.ndarray, min_obs: int = 10, 
                              n_workers: int = None) -> np.ndarray:
    n_samples, n_pos, _ = kidera.shape
    cov_mat = np.full((n_pos, n_pos), np.nan)
    
    if n_workers is None:
        n_workers = max(1, os.cpu_count() - 2)
    
    shm = shared_memory.SharedMemory(create=True, size=kidera.nbytes)
    shm_arr = np.ndarray(kidera.shape, dtype=kidera.dtype, buffer=shm.buf)
    np.copyto(shm_arr, kidera)
    
    try:
        with Pool(n_workers, initializer=_init_worker, 
                  initargs=(shm.name, kidera.shape, kidera.dtype)) as pool:
            
            results = list(tqdm(
                pool.imap(_compute_row, range(n_pos)),
                total=n_pos,
                desc=f"Covariance ({n_workers} workers)"
            ))
        
        for row_results in results:
            for i, j, rv in row_results:
                cov_mat[i, j] = rv
                cov_mat[j, i] = rv
    finally:
        shm.close()
        shm.unlink()
    
    return cov_mat

# ---------------------------------------------------------------------
# MULTI-GENE COVARIANCE
# ---------------------------------------------------------------------

def calc_multigene_covariance(fasta_paths: list[str], gene_names: list[str],
                               max_gap_frac: float = 0.5, min_samples: int = 100,
                               min_var: float = 0.1, n_workers: int = None, 
                               save_path: str = None) -> dict:
    assert len(fasta_paths) == len(gene_names)
    
    all_seqs = []
    all_ids = []
    
    for fasta in fasta_paths:
        seqs, ids, _ = read_alignment(fasta, max_gap_frac)
        all_seqs.append(seqs)
        all_ids.append(set(ids))
    
    common = set.intersection(*all_ids)
    assert len(common) >= min_samples, f"Only {len(common)} shared samples"
    print(f"Found {len(common)} samples common to all {len(gene_names)} genes")
    
    seqs_aligned = []
    for i, fasta in enumerate(fasta_paths):
        seqs, ids, _ = read_alignment(fasta, max_gap_frac)
        id_to_idx = {x: j for j, x in enumerate(ids)}
        
        if i == 0:
            id_order = [x for x in ids if x in common]
        
        idx_sorted = [id_to_idx[x] for x in id_order]
        seqs_aligned.append(seqs[idx_sorted])
    
    # Encode and filter by variance
    kidera_filtered = []
    gene_lengths_filtered = []
    
    for i, seqs in enumerate(seqs_aligned):
        kidera = encode_kidera(seqs)
        kidera_filt, kept_pos = filter_low_variance_positions(kidera, min_var)
        kidera_filtered.append(kidera_filt)
        gene_lengths_filtered.append(kidera_filt.shape[1])
        print(f"  {gene_names[i]}: {seqs.shape[1]} -> {kidera_filt.shape[1]} positions (var >= {min_var})")
    
    kidera_concat = np.concatenate(kidera_filtered, axis=1)
    print(f"Concatenated: {kidera_concat.shape[0]} samples x {kidera_concat.shape[1]} positions")
    
    cov_mat = calc_position_covariance(kidera_concat, n_workers=n_workers)
    
    gene_bounds = np.cumsum([0] + gene_lengths_filtered)
    
    res = {
        "cov": cov_mat,
        "gene_names": gene_names,
        "gene_lengths": gene_lengths_filtered,
        "gene_bounds": gene_bounds,
        "ids": id_order,
        "n": len(id_order)
    }
    
    if save_path:
        np.savez(save_path, 
                 cov=cov_mat,
                 gene_names=np.array(gene_names),
                 gene_lengths=np.array(gene_lengths_filtered),
                 gene_bounds=gene_bounds,
                 ids=np.array(id_order),
                 n=np.array(res["n"]))
        print(f"Saved to {save_path}")
    
    return res

# ---------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------

def summarize_covariance(res: dict) -> pl.DataFrame:
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

# ---------------------------------------------------------------------
# PLOTTING (all save to file)
# ---------------------------------------------------------------------

def plot_multigene_covariance(res: dict, outfile: str = "data/tmp/plots/multigene_cov.png"):
    import matplotlib.pyplot as plt
    
    names = res["gene_names"]
    bounds = res["gene_bounds"]
    n_genes = len(names)
    
    fig, ax = plt.subplots(figsize=(14, 14))
    im = ax.imshow(res["cov"], cmap="YlOrRd", aspect="auto")
    plt.colorbar(im, ax=ax, label="RV coefficient", shrink=0.8)
    
    for b in bounds[1:-1]:
        ax.axhline(b - 0.5, color="black", linewidth=0.5, alpha=0.5)
        ax.axvline(b - 0.5, color="black", linewidth=0.5, alpha=0.5)
    
    mids = [(bounds[i] + bounds[i + 1]) / 2 for i in range(n_genes)]
    ax.set_xticks(mids)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_yticks(mids)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_title(f"Multi-gene covariance (n={res['n']}, {n_genes} genes)")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved: {outfile}")

def plot_gene_level_heatmap(res: dict, outfile: str = "data/tmp/plots/gene_level_cov.png"):
    import matplotlib.pyplot as plt
    
    names = res["gene_names"]
    bounds = res["gene_bounds"]
    cov = res["cov"]
    n_genes = len(names)
    
    gene_cov = np.full((n_genes, n_genes), np.nan)
    for i in range(n_genes):
        for j in range(n_genes):
            start_i, end_i = bounds[i], bounds[i + 1]
            start_j, end_j = bounds[j], bounds[j + 1]
            block = cov[start_i:end_i, start_j:end_j]
            vals = block.flatten()
            vals = vals[~np.isnan(vals)]
            gene_cov[i, j] = np.mean(vals) if len(vals) > 0 else np.nan
    
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(gene_cov, cmap="YlOrRd", aspect="auto")
    plt.colorbar(im, ax=ax, label="Mean RV coefficient")
    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title("Gene-level mean covariance")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Saved: {outfile}")

def plot_clustered_heatmap(res: dict, outfile: str = "data/tmp/plots/gene_clustered_cov.png"):
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform
    
    names = res["gene_names"]
    bounds = res["gene_bounds"]
    cov = res["cov"]
    n_genes = len(names)
    
    # Gene-level covariance
    gene_cov = np.full((n_genes, n_genes), np.nan)
    for i in range(n_genes):
        for j in range(n_genes):
            start_i, end_i = bounds[i], bounds[i + 1]
            start_j, end_j = bounds[j], bounds[j + 1]
            block = cov[start_i:end_i, start_j:end_j]
            vals = block.flatten()
            vals = vals[~np.isnan(vals)]
            gene_cov[i, j] = np.mean(vals) if len(vals) > 0 else np.nan
    
    # Distance and cluster
    gene_dist = 1 - gene_cov
    np.fill_diagonal(gene_dist, 0)
    gene_dist = np.nan_to_num(gene_dist, nan=1.0)
    gene_dist = (gene_dist + gene_dist.T) / 2
    
    Z = linkage(squareform(gene_dist), method="ward")
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    
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
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Saved: {outfile}")

def plot_within_vs_between(res: dict, outfile: str = "data/tmp/plots/within_vs_between.png"):
    import matplotlib.pyplot as plt
    
    cov = res["cov"]
    names = res["gene_names"]
    bounds = res["gene_bounds"]
    n_genes = len(names)
    
    within_vals = []
    between_vals = []
    
    for i in range(n_genes):
        for j in range(i, n_genes):
            start_i, end_i = bounds[i], bounds[i + 1]
            start_j, end_j = bounds[j], bounds[j + 1]
            block = cov[start_i:end_i, start_j:end_j]
            
            if i == j:
                mask = np.triu(np.ones_like(block, dtype=bool), k=1)
                vals = block[mask]
                within_vals.extend(vals[~np.isnan(vals)].tolist())
            else:
                vals = block.flatten()
                between_vals.extend(vals[~np.isnan(vals)].tolist())
    
    within_vals = np.array(within_vals)
    between_vals = np.array(between_vals)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    axes[0].hist(within_vals, bins=50, alpha=0.6, label=f"Within (n={len(within_vals):,})", density=True)
    axes[0].hist(between_vals, bins=50, alpha=0.6, label=f"Between (n={len(between_vals):,})", density=True)
    axes[0].set_xlabel("RV coefficient")
    axes[0].set_ylabel("Density")
    axes[0].legend()
    axes[0].set_title("Within vs Between gene RV")
    
    # Boxplot
    axes[1].boxplot([within_vals, between_vals], labels=["Within", "Between"])
    axes[1].set_ylabel("RV coefficient")
    axes[1].set_title(f"Within: {np.mean(within_vals):.3f} vs Between: {np.mean(between_vals):.3f}")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Saved: {outfile}")

def plot_sanity_checks(res: dict, outfile: str = "data/tmp/plots/sanity_checks.png"):
    import matplotlib.pyplot as plt
    
    cov = res["cov"]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Distribution
    vals = cov.flatten()
    vals = vals[~np.isnan(vals)]
    axes[0, 0].hist(vals, bins=100)
    axes[0, 0].set_xlabel("RV coefficient")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title(f"Distribution (n={len(vals):,})")
    
    # Log scale
    axes[0, 1].hist(vals, bins=100, log=True)
    axes[0, 1].set_xlabel("RV coefficient")
    axes[0, 1].set_ylabel("Count (log)")
    axes[0, 1].set_title("Log scale")
    
    # Cumulative near 1
    thresholds = np.linspace(0.5, 1.0, 50)
    pct = [100 * np.sum(vals >= t) / len(vals) for t in thresholds]
    axes[1, 0].plot(thresholds, pct)
    axes[1, 0].set_xlabel("Threshold")
    axes[1, 0].set_ylabel("% values >= threshold")
    axes[1, 0].set_title("Cumulative from 0.5")
    
    # Stats text
    stats_text = f"""
    Matrix shape: {cov.shape}
    NaN: {np.isnan(cov).sum()} ({100*np.isnan(cov).sum()/cov.size:.1f}%)
    Min: {np.nanmin(cov):.4f}
    Max: {np.nanmax(cov):.4f}
    Mean: {np.nanmean(cov):.4f}
    Median: {np.nanmedian(cov):.4f}
    
    RV >= 0.99: {100*np.sum(vals >= 0.99)/len(vals):.2f}%
    RV >= 0.9: {100*np.sum(vals >= 0.9)/len(vals):.2f}%
    RV >= 0.5: {100*np.sum(vals >= 0.5)/len(vals):.2f}%
    """
    axes[1, 1].text(0.1, 0.5, stats_text, fontsize=10, family="monospace", va="center")
    axes[1, 1].axis("off")
    axes[1, 1].set_title("Summary stats")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Saved: {outfile}")


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--min-var", type=float, default=0.1)
    args = parser.parse_args()
    
    os.makedirs("data/tmp/plots", exist_ok=True)
    os.makedirs("data/tmp/covariance", exist_ok=True)
    
    gene_dir = Path("data/tmp/alignedGenes")
    genes = [
        "atpA", "atpB", "atpE", "atpF", "atpH", "atpI",
        "ccsA", "cemA",
        "matK",
        "ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhG", "ndhH", "ndhI", "ndhJ", "ndhK",
        "petA", "petG", "petN",
        "psaA", "psaB", "psaC", "psaJ",
        "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbM", "psbN", "psbT", "psbZ",
        "rbcL",
        "rpl14", "rpl16", "rpl20", "rpl23", "rpl2", "rpl33", "rpl36",
        "rpoA", "rpoB", "rpoC1",
        "rps11", "rps14", "rps18", "rps19", "rps3", "rps4", "rps7", "rps8",
        "ycf3", "ycf4",
    ]
    fasta_paths = [str(gene_dir / f"{g}_AA_aligned.fasta") for g in genes]
    
    # Compute
    res = calc_multigene_covariance(
        fasta_paths, genes,
        min_var=args.min_var,
        n_workers=args.workers,
        save_path="data/tmp/covariance/multigene_filtered.npz"
    )
    
    # Summary table
    summary = summarize_covariance(res)
    print(summary)
    summary.write_csv("data/tmp/covariance/summary.csv")
    
    # All plots
    plot_sanity_checks(res)
    plot_multigene_covariance(res)
    plot_gene_level_heatmap(res)
    plot_clustered_heatmap(res)
    plot_within_vs_between(res)
    
    print("\nDone!")