import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import hashlib

def prepare_siamese_dataset(clean_seqs_file, pheno_file, min_temp_diff=500):
    seqs = pd.read_parquet(clean_seqs_file)
    pheno = pd.read_parquet(pheno_file)
    df = seqs.merge(pheno, on="ID").drop_duplicates("ID")

    # Pre-filter sequences once
    valid_mask = df["AA_seq"].str.len().ge(20) & df["AA_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY*X-]+")
    df = df[valid_mask].reset_index(drop=True)

    ids = df["ID"].values
    temps = df["pheno_Topt_site_p50"].values
    seqs = df["AA_seq"].values
    genes = df["Gene"].values

    diffs = temps[:, None] - temps
    mask = np.triu(np.abs(diffs) >= min_temp_diff, k=1)
    before = mask.sum()
    mask &= genes[:, None] == genes
    after = mask.sum()

    print("Pairs removed:", before - after)
    print("Remaining pairs:", after)  # keep only within-gene pairs

    i, j = np.where(mask)

    hashes = np.array([hash(s) for s in seqs], dtype=np.int64)
    mask_diff = hashes[i] != hashes[j]
    i, j = i[mask_diff], j[mask_diff]

    print("computed pairs:", len(i))

    out = pd.DataFrame({
        "Gene1": genes[i],
        "Gene2": genes[j],
        "Species1": ids[i],
        "Species2": ids[j],
        "seq1": seqs[i],
        "seq2": seqs[j],
        "label": np.where(diffs[i, j] > 0, 1, 0)
    })

    out = out[out["Gene1"] == out["Gene2"]].reset_index(drop=True)

    print("computed pairs:", len(out))

    return out

def split_and_save(df, train_ratio=0.8, output_prefix="siamese"):
    """Split data and save train/validation TSVs"""
    
    # Stratified split to ensure temperature range coverage
    temp_bins = pd.cut(df['label'], bins=10, labels=False)
    
    train_df, val_df = train_test_split(
        df, 
        test_size=1-train_ratio, 
        stratify=temp_bins,
        random_state=42
    )
    
    print(f"Train set: {len(train_df)} pairs")
    print(f"Validation set: {len(val_df)} pairs")
    
    # Save files
    train_file = f"{output_prefix}_train.tsv"
    val_file = f"{output_prefix}_valid.tsv"
    
    train_df.to_csv(train_file, sep='\t', index=False)
    val_df.to_csv(val_file, sep='\t', index=False)
    
    print(f"Saved {train_file} and {val_file}")
    
    return train_file, val_file

# Example usage:
if __name__ == "__main__":
    # Your file paths
    clean_seqs_file = "data/clean_AA_seqs.parquet"
    pheno_file = "data/pheno_topt_clean.parquet"
    
    print("Processing chloroplast genome data for temperature modeling...")
    print(f"Sequences file: {clean_seqs_file}")
    print(f"Phenotype file: {pheno_file}")
    
    # Prepare siamese dataset
    siamese_df = prepare_siamese_dataset(clean_seqs_file, pheno_file, min_temp_diff=500)
    
    # Split and save
    train_file, val_file = split_and_save(siamese_df, output_prefix="chloroplast_temp")
    
    print(f"\nDataset ready!")
    print(f"Train command: python siamese_esm.py -train {train_file} -valid {val_file} -pretrained facebook/esm2_t33_650M_UR50D -batch_size 8 -epochs 10 -output ./models -project chloroplast_temp -id test_run -lr 1e-5")