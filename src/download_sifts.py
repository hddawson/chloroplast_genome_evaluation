
#!/usr/bin/env python

import pandas as pd

import urllib.request

import gzip

import shutil

from pathlib import Path

# Create directory

sifts_dir = Path("data/SIFTS")

sifts_dir.mkdir(parents=True, exist_ok=True)

print("Downloading SIFTS mapping file...")

url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"

gz_file = sifts_dir / "pdb_chain_uniprot.csv.gz"

csv_file = sifts_dir / "pdb_chain_uniprot_plus.csv"

# Download

urllib.request.urlretrieve(url, gz_file)

print("Extracting...")

with gzip.open(gz_file, 'rb') as f_in:

    with open(csv_file, 'wb') as f_out:

        shutil.copyfileobj(f_in, f_out)

print(f"SIFTS mapping saved to {csv_file}")

# Create a dummy FASTA file (empty or minimal) since you're using existing alignments

fasta_file = sifts_dir / "pdb_chain_uniprot_plus.fasta"

fasta_file.touch()

print(f"Created placeholder FASTA file at {fasta_file}")

print("Done! SIFTS database updated.")

