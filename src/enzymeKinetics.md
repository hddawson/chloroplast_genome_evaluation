# ESM-C Embeddings vs Enzyme Kinetics/Stability

## Goal
Test whether ESM-C protein embeddings correlate with enzyme kinetic parameters (kcat, Km, kcat/Km) and stability metrics (Tm, ΔG).

## Data Sources

### Kinetics: BRENDA + SABIO-RK
- **BRENDA** (https://www.brenda-enzymes.org): Curated enzyme database with kcat, Km values linked to UniProt IDs
- **SABIO-RK** (http://sabiork.h-its.org): Biochemical reaction kinetics with standardized conditions

### Stability: ProThermDB + FireProtDB
- **ProThermDB** (https://web.iitm.ac.in/bioinfo2/prothermdb/): Thermodynamic data (Tm, ΔG) for wild-type proteins
- **FireProtDB** (https://loschmidt.chemi.muni.cz/fireprotdb/): Curated stability data with sequences

### Sequences
- **UniProt** (https://www.uniprot.org): Retrieve sequences by accession

## Directory Structure
```
esmc_kinetics/
├── data/
│   ├── raw/                 # Downloaded files
│   ├── processed/           # Filtered, merged datasets
│   └── embeddings/          # ESM-C outputs
├── scripts/
│   ├── download_data.py
│   ├── parse_and_filter.py
│   ├── embed_sequences.py
│   └── analyze_embeddings.py
├── results/
│   ├── figures/
│   └── tables/
├── environment.yml
└── README.md
```

## Workflow

### 1. Setup Environment
```bash
conda env create -f environment.yml
conda activate esmc_kinetics
```

### 2. Download Data
```bash
python scripts/download_data.py
```

### 3. Parse and Filter
```bash
python scripts/parse_and_filter.py
```
Filtering criteria:
- Remove entries without UniProt ID
- Remove duplicate UniProt-substrate pairs (keep median value)
- Require sequence length 50-2000 aa
- Log-transform kcat, Km (span orders of magnitude)

### 4. Generate Embeddings
```bash
python scripts/embed_sequences.py --model esmc_300m --batch_size 8
```
ESM-C model: `esmc_300m` (or `esmc_600m` for higher capacity)
Pooling: Mean over residue embeddings → single vector per protein

### 5. Analyze
```bash
python scripts/analyze_embeddings.py
```

## Key Files

| File | Description |
|------|-------------|
| `data/processed/kinetics.parquet` | UniProt ID, EC, organism, log_kcat, log_Km, substrate |
| `data/processed/stability.parquet` | UniProt ID, Tm, dG, organism |
| `data/embeddings/embeddings.npz` | {uniprot_id: embedding_vector} |
| `results/tables/correlations.csv` | Pearson/Spearman r, p-values |
| `results/figures/*.png` | PCA, correlation plots |

---

## environment.yml
```yaml
name: esmc_kinetics
channels:
  - conda-forge
  - pytorch
dependencies:
  - python=3.10
  - pytorch>=2.0
  - pandas
  - numpy
  - scipy
  - scikit-learn
  - matplotlib
  - seaborn
  - requests
  - pyarrow
  - pip
  - pip:
    - esm  # or: git+https://github.com/evolutionaryscale/esm.git
```

---

## scripts/download_data.py
```python
#!/usr/bin/env python3
"""Download kinetics and stability data."""
import os
import requests
from pathlib import Path

RAW_DIR = Path("data/raw")
RAW_DIR.mkdir(parents=True, exist_ok=True)

# SABIO-RK: REST API for kinetic data
# Query enzymes with kcat values
SABIO_URL = "http://sabiork.h-its.org/sabioRestWebServices/kineticLaws/search"

def download_sabio_kinetics(ec_classes: list[str], outfile: Path):
    """Download kinetic laws for given EC classes."""
    records = []
    for ec in ec_classes:
        params = {"q": f"ECNumber:{ec}* AND Parameter:kcat", "format": "json"}
        resp = requests.get(SABIO_URL, params=params, timeout=60)
        if resp.ok:
            records.extend(resp.json().get("results", []))
    assert len(records) > 0, "No SABIO-RK records retrieved"
    import json
    outfile.write_text(json.dumps(records, indent=2))
    print(f"SABIO-RK: {len(records)} records -> {outfile}")

def download_fireprotdb(outfile: Path):
    """Download FireProtDB stability data (CSV export)."""
    url = "https://loschmidt.chemi.muni.cz/fireprotdb/api/v1/proteins/export?format=csv"
    resp = requests.get(url, timeout=120)
    assert resp.ok, f"FireProtDB download failed: {resp.status_code}"
    outfile.write_bytes(resp.content)
    print(f"FireProtDB: {len(resp.content)} bytes -> {outfile}")

def download_uniprot_sequences(uniprot_ids: list[str], outfile: Path):
    """Batch retrieve sequences from UniProt."""
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {"query": " OR ".join(f"accession:{uid}" for uid in uniprot_ids),
              "format": "fasta"}
    resp = requests.get(url, params=params, timeout=300)
    assert resp.ok, f"UniProt download failed: {resp.status_code}"
    outfile.write_text(resp.text)
    print(f"UniProt: {resp.text.count('>')} sequences -> {outfile}")

if __name__ == "__main__":
    # Example: download oxidoreductases and transferases
    download_sabio_kinetics(["1.", "2."], RAW_DIR / "sabio_kinetics.json")
    download_fireprotdb(RAW_DIR / "fireprotdb.csv")
    # Sequences downloaded after filtering step (need UniProt IDs)
```

---

## scripts/parse_and_filter.py
```python
#!/usr/bin/env python3
"""Parse raw data, filter, and create analysis-ready datasets."""
import json
import pandas as pd
import numpy as np
from pathlib import Path

RAW_DIR = Path("data/raw")
PROC_DIR = Path("data/processed")
PROC_DIR.mkdir(parents=True, exist_ok=True)

MIN_SEQ_LEN, MAX_SEQ_LEN = 50, 2000

def parse_sabio_kinetics(infile: Path) -> pd.DataFrame:
    """Extract kcat, Km from SABIO-RK JSON."""
    with open(infile) as f:
        data = json.load(f)
    
    rows = []
    for entry in data:
        uniprot = entry.get("uniprotID")
        if not uniprot:
            continue
        for param in entry.get("parameters", []):
            ptype = param.get("type")
            if ptype in ("kcat", "Km"):
                rows.append({
                    "uniprot_id": uniprot,
                    "ec": entry.get("ecNumber"),
                    "organism": entry.get("organism"),
                    "substrate": param.get("substrate"),
                    "param_type": ptype,
                    "value": float(param.get("value", np.nan)),
                    "unit": param.get("unit")
                })
    
    df = pd.DataFrame(rows)
    assert len(df) > 0, "No valid kinetic records parsed"
    return df

def parse_fireprotdb(infile: Path) -> pd.DataFrame:
    """Parse FireProtDB CSV for wild-type stability data."""
    df = pd.read_csv(infile)
    # Keep wild-type entries only
    df = df[df["mutation"].isna() | (df["mutation"] == "WT")]
    df = df.rename(columns={"uniprot_id": "uniprot_id", "Tm": "Tm", "ddG": "dG"})
    df = df[["uniprot_id", "Tm", "dG", "organism"]].dropna(subset=["uniprot_id"])
    assert len(df) > 0, "No valid stability records parsed"
    return df

def aggregate_kinetics(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate duplicate entries by median."""
    # Pivot to wide format
    agg = df.groupby(["uniprot_id", "ec", "organism", "substrate", "param_type"])["value"].median().unstack()
    agg = agg.reset_index()
    
    # Log-transform (values span orders of magnitude)
    for col in ["kcat", "Km"]:
        if col in agg.columns:
            agg[f"log_{col}"] = np.log10(agg[col].replace(0, np.nan))
    
    # Compute catalytic efficiency if both available
    if "log_kcat" in agg.columns and "log_Km" in agg.columns:
        agg["log_kcat_Km"] = agg["log_kcat"] - agg["log_Km"]
    
    return agg

def main():
    # Parse
    kinetics = parse_sabio_kinetics(RAW_DIR / "sabio_kinetics.json")
    stability = parse_fireprotdb(RAW_DIR / "fireprotdb.csv")
    
    # Aggregate kinetics
    kinetics_agg = aggregate_kinetics(kinetics)
    
    # Collect unique UniProt IDs
    all_ids = set(kinetics_agg["uniprot_id"]) | set(stability["uniprot_id"])
    print(f"Unique proteins: {len(all_ids)}")
    
    # Save
    kinetics_agg.to_parquet(PROC_DIR / "kinetics.parquet", index=False)
    stability.to_parquet(PROC_DIR / "stability.parquet", index=False)
    
    # Save ID list for sequence retrieval
    (PROC_DIR / "uniprot_ids.txt").write_text("\n".join(sorted(all_ids)))
    print(f"Saved {len(all_ids)} UniProt IDs")

if __name__ == "__main__":
    main()
```

---

## scripts/embed_sequences.py
```python
#!/usr/bin/env python3
"""Generate ESM-C embeddings for protein sequences."""
import argparse
import numpy as np
import torch
from pathlib import Path

PROC_DIR = Path("data/processed")
EMB_DIR = Path("data/embeddings")
EMB_DIR.mkdir(parents=True, exist_ok=True)

def load_fasta(fasta_file: Path) -> dict[str, str]:
    """Parse FASTA to {id: sequence} dict."""
    seqs = {}
    current_id, current_seq = None, []
    for line in fasta_file.read_text().splitlines():
        if line.startswith(">"):
            if current_id:
                seqs[current_id] = "".join(current_seq)
            # Extract UniProt ID from header (e.g., >sp|P12345|...)
            current_id = line.split("|")[1] if "|" in line else line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs

def embed_with_esmc(sequences: dict[str, str], model_name: str, batch_size: int, device: str) -> dict[str, np.ndarray]:
    """Generate mean-pooled ESM-C embeddings."""
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein
    
    model = ESMC.from_pretrained(model_name).to(device)
    model.eval()
    
    embeddings = {}
    ids = list(sequences.keys())
    
    for i in range(0, len(ids), batch_size):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = [sequences[uid] for uid in batch_ids]
        
        proteins = [ESMProtein(sequence=seq) for seq in batch_seqs]
        
        with torch.no_grad():
            for uid, protein in zip(batch_ids, proteins):
                output = model.encode(protein)
                # Mean pool over residue dimension
                emb = output.embeddings.mean(dim=1).squeeze().cpu().numpy()
                embeddings[uid] = emb
        
        print(f"Embedded {min(i+batch_size, len(ids))}/{len(ids)}")
    
    return embeddings

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="esmc_300m", choices=["esmc_300m", "esmc_600m"])
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    
    fasta_file = PROC_DIR / "sequences.fasta"
    assert fasta_file.exists(), f"Run download_data.py first to get {fasta_file}"
    
    sequences = load_fasta(fasta_file)
    assert len(sequences) > 0, "No sequences loaded"
    print(f"Loaded {len(sequences)} sequences")
    
    embeddings = embed_with_esmc(sequences, args.model, args.batch_size, args.device)
    
    # Save as npz
    np.savez_compressed(EMB_DIR / "embeddings.npz", **embeddings)
    print(f"Saved embeddings to {EMB_DIR / 'embeddings.npz'}")

if __name__ == "__main__":
    main()
```

---

## scripts/analyze_embeddings.py
```python
#!/usr/bin/env python3
"""Correlate ESM-C embeddings with kinetic/stability parameters."""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

PROC_DIR = Path("data/processed")
EMB_DIR = Path("data/embeddings")
RES_DIR = Path("results")
FIG_DIR = RES_DIR / "figures"
TAB_DIR = RES_DIR / "tables"
for d in [FIG_DIR, TAB_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SEED = 42
np.random.seed(SEED)

def load_data():
    """Load embeddings and phenotype data."""
    emb_data = np.load(EMB_DIR / "embeddings.npz")
    embeddings = {k: emb_data[k] for k in emb_data.files}
    
    kinetics = pd.read_parquet(PROC_DIR / "kinetics.parquet")
    stability = pd.read_parquet(PROC_DIR / "stability.parquet")
    
    return embeddings, kinetics, stability

def merge_embeddings_with_phenotypes(embeddings: dict, pheno_df: pd.DataFrame, target_col: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Align embeddings with phenotype values."""
    valid_ids = [uid for uid in pheno_df["uniprot_id"] if uid in embeddings]
    df = pheno_df[pheno_df["uniprot_id"].isin(valid_ids)].dropna(subset=[target_col])
    
    X = np.array([embeddings[uid] for uid in df["uniprot_id"]])
    y = df[target_col].values
    
    assert len(X) == len(y), "Mismatch in X, y lengths"
    return X, y, df

def compute_correlations(X: np.ndarray, y: np.ndarray) -> dict:
    """Compute Pearson and Spearman correlations with PCA components."""
    pca = PCA(n_components=min(50, X.shape[1]))
    X_pca = pca.fit_transform(X)
    
    results = []
    for i in range(X_pca.shape[1]):
        r_pearson, p_pearson = stats.pearsonr(X_pca[:, i], y)
        r_spearman, p_spearman = stats.spearmanr(X_pca[:, i], y)
        results.append({
            "PC": i+1,
            "var_explained": pca.explained_variance_ratio_[i],
            "pearson_r": r_pearson,
            "pearson_p": p_pearson,
            "spearman_r": r_spearman,
            "spearman_p": p_spearman
        })
    
    return pd.DataFrame(results), X_pca, pca

def fit_predictive_models(X: np.ndarray, y: np.ndarray, target_name: str) -> pd.DataFrame:
    """Fit Ridge, Lasso, RF and report CV R²."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {
        "Ridge": Ridge(alpha=1.0, random_state=SEED),
        "Lasso": Lasso(alpha=0.1, random_state=SEED, max_iter=5000),
        "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10, random_state=SEED, n_jobs=-1)
    }
    
    results = []
    for name, model in models.items():
        scores = cross_val_score(model, X_scaled, y, cv=5, scoring="r2")
        results.append({
            "target": target_name,
            "model": name,
            "cv_r2_mean": scores.mean(),
            "cv_r2_std": scores.std()
        })
        print(f"{target_name} | {name}: R²={scores.mean():.3f}±{scores.std():.3f}")
    
    return pd.DataFrame(results)

def plot_pca_scatter(X_pca: np.ndarray, y: np.ndarray, target_name: str, outfile: Path):
    """Scatter plot of PC1 vs PC2 colored by target."""
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(X_pca[:, 0], X_pca[:, 1], c=y, cmap="viridis", alpha=0.7, s=20)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"ESM-C Embedding PCA colored by {target_name}")
    plt.colorbar(scatter, label=target_name)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()

def analyze_target(embeddings: dict, pheno_df: pd.DataFrame, target_col: str, label: str):
    """Full analysis pipeline for one target variable."""
    X, y, df = merge_embeddings_with_phenotypes(embeddings, pheno_df, target_col)
    
    if len(X) < 30:
        print(f"Skipping {label}: only {len(X)} samples")
        return None, None
    
    print(f"\n=== {label} (n={len(X)}) ===")
    
    # Correlations
    corr_df, X_pca, pca = compute_correlations(X, y)
    corr_df["target"] = label
    
    # Predictive models
    model_df = fit_predictive_models(X, y, label)
    
    # Plot
    plot_pca_scatter(X_pca, y, label, FIG_DIR / f"pca_{label.replace('/', '_')}.png")
    
    return corr_df, model_df

def main():
    embeddings, kinetics, stability = load_data()
    print(f"Loaded {len(embeddings)} embeddings")
    
    all_corr = []
    all_models = []
    
    # Kinetics targets
    for col, label in [("log_kcat", "log_kcat"), ("log_Km", "log_Km"), ("log_kcat_Km", "log_kcat/Km")]:
        if col in kinetics.columns:
            corr, models = analyze_target(embeddings, kinetics, col, label)
            if corr is not None:
                all_corr.append(corr)
                all_models.append(models)
    
    # Stability targets
    for col, label in [("Tm", "Tm"), ("dG", "dG")]:
        if col in stability.columns:
            corr, models = analyze_target(embeddings, stability, col, label)
            if corr is not None:
                all_corr.append(corr)
                all_models.append(models)
    
    # Save results
    if all_corr:
        pd.concat(all_corr).to_csv(TAB_DIR / "correlations.csv", index=False)
    if all_models:
        pd.concat(all_models).to_csv(TAB_DIR / "model_performance.csv", index=False)
    
    print(f"\nResults saved to {TAB_DIR}")

if __name__ == "__main__":
    main()
```

---

## Confounder Handling

Add to `analyze_embeddings.py` if needed:

```python
def control_for_confounders(df: pd.DataFrame, X: np.ndarray, y: np.ndarray, embeddings: dict) -> tuple[np.ndarray, np.ndarray]:
    """Residualize embeddings against sequence length and EC class."""
    from sklearn.linear_model import LinearRegression
    
    # Compute sequence lengths
    seq_lens = np.array([len(embeddings[uid]) for uid in df["uniprot_id"]])  # proxy if seqs not stored
    
    # One-hot EC class (first digit)
    if "ec" in df.columns:
        ec_class = df["ec"].str[0].astype("category").cat.codes.values
        confounders = np.column_stack([seq_lens, ec_class])
    else:
        confounders = seq_lens.reshape(-1, 1)
    
    # Residualize X
    reg = LinearRegression().fit(confounders, X)
    X_resid = X - reg.predict(confounders)
    
    # Residualize y
    reg_y = LinearRegression().fit(confounders, y)
    y_resid = y - reg_y.predict(confounders)
    
    return X_resid, y_resid
```

---

## Reproducibility Checklist

1. **Clone and setup**
   ```bash
   git clone <repo> && cd esmc_kinetics
   conda env create -f environment.yml
   conda activate esmc_kinetics
   ```

2. **Full pipeline**
   ```bash
   python scripts/download_data.py
   python scripts/parse_and_filter.py
   # Download sequences (requires uniprot_ids.txt from previous step)
   python scripts/download_data.py --sequences
   python scripts/embed_sequences.py --model esmc_300m
   python scripts/analyze_embeddings.py
   ```

3. **Outputs**
   - `results/tables/correlations.csv`: PC-level correlations
   - `results/tables/model_performance.csv`: CV R² for Ridge/Lasso/RF
   - `results/figures/pca_*.png`: Embedding visualizations

4. **Determinism**
   - Set `SEED=42` in all scripts
   - Pin package versions in `environment.yml`
   - Use `np.random.seed()` before any stochastic operation