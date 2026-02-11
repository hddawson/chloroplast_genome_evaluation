import pandas as pd
import numpy as np
from umap import UMAP
import matplotlib.pyplot as plt

# Load data
pca = pd.read_parquet("results/embeddings_full_pca.parquet")
data = pd.read_parquet("data/processed_data.parquet")

# Extract PC columns
pc_cols = [c for c in pca.columns if c.startswith("PC") and c[2:].isdigit()]
assert len(pc_cols) > 0, "No PC columns found"

# Run UMAP
umap_model = UMAP(n_neighbors=15, min_dist=0.1, random_state=42, n_jobs=-1)
umap_coords = umap_model.fit_transform(pca[pc_cols].values)

pca["UMAP1"] = umap_coords[:, 0]
pca["UMAP2"] = umap_coords[:, 1]

# Merge metadata
pca = pca.merge(data[["ID", "Organism", "Order", "pheno_wc2.1_2.5m_bio_8_p50"]], on="ID", how="left")

# Plot by Gene
fig, ax = plt.subplots(figsize=(10, 8))
for gene in pca["Gene"].unique():
    subset = pca[pca["Gene"] == gene]
    ax.scatter(subset["UMAP1"], subset["UMAP2"], s=0.5, alpha=0.3, label=gene)
ax.set_title("UMAP colored by Gene")
ax.legend(markerscale=5)
plt.savefig("plots/umap_by_gene.png", dpi=150)

# Plot by Residue
fig, ax = plt.subplots(figsize=(10, 8))
for res in pca["Residue"].unique():
    subset = pca[pca["Residue"] == res]
    ax.scatter(subset["UMAP1"], subset["UMAP2"], s=0.5, alpha=0.3, label=res)
ax.set_title("UMAP colored by Residue")
plt.savefig("plots/umap_by_residue.png", dpi=150)

# Plot by phenotype (temperature)
fig, ax = plt.subplots(figsize=(10, 8))
sc = ax.scatter(
    pca["UMAP1"], pca["UMAP2"],
    c=pca["pheno_wc2.1_2.5m_bio_8_p50"],
    s=0.5, alpha=0.3, cmap="viridis"
)
plt.colorbar(sc, label="Temp (BIO8 p50)")
ax.set_title("UMAP colored by Temperature")
plt.savefig("plots/umap_by_temperature.png", dpi=150)