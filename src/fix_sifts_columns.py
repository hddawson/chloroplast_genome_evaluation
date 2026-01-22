#!/usr/bin/env python
import pandas as pd
import numpy as np
from pathlib import Path

# Read the SIFTS file
sifts_file = Path("data/SIFTS/pdb_chain_uniprot_plus_original.csv")
output_file = Path("data/SIFTS/pdb_chain_uniprot_plus.csv")

print("Reading SIFTS file...")
df = pd.read_csv(sifts_file, comment='#', low_memory=False)

print(f"Original columns: {df.columns.tolist()}")
print(f"Total rows: {len(df)}")

# Rename columns to match what evCouplings expects
column_mapping = {
    'PDB': 'pdb_id',
    'CHAIN': 'pdb_chain',
    'SP_PRIMARY': 'uniprot_ac',
    'RES_BEG': 'resseq_start',
    'RES_END': 'resseq_end',
    'PDB_BEG': 'pdb_start',
    'PDB_END': 'pdb_end',
    'SP_BEG': 'uniprot_start',
    'SP_END': 'uniprot_end'
}

df = df.rename(columns=column_mapping)

# Add uniprot_id as duplicate of uniprot_ac
df['uniprot_id'] = df['uniprot_ac']

# Convert numeric columns to proper types
# Handle NaN and convert to appropriate numeric types
numeric_cols = ['resseq_start', 'resseq_end', 'pdb_start', 'pdb_end', 
                'uniprot_start', 'uniprot_end']

for col in numeric_cols:
    # Convert to numeric, coercing errors to NaN
    df[col] = pd.to_numeric(df[col], errors='coerce')
    # Convert to Int64 (nullable integer type) to handle NaN
    df[col] = df[col].astype('Int64')

print(f"New columns: {df.columns.tolist()}")
print(f"\nData types after conversion:")
for col in df.columns:
    print(f"  {col}: {df[col].dtype}")

# Save the fixed file
print(f"\nSaving to {output_file}...")
df.to_csv(output_file, index=False)

print("Done! SIFTS file has been fixed.")
print(f"\nFirst few rows of fixed file:")
df_check = pd.read_csv(output_file, nrows=3)
print(df_check)
print(f"\nSample of rows with NaN values:")
print(df[df['pdb_end'].isna()].head(3))
