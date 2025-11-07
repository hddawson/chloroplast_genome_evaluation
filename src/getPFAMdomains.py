import os
import subprocess
import pandas as pd
from Bio import SeqIO

# Paths
#fasta_dir = "data/tmp/seqSmall"
fasta_dir = "data/AA_seqs"
pfam_db = "data/Pfam-A.hmm"
output_dir = "pfam_results"
os.makedirs(output_dir, exist_ok=True)

results = []

# Iterate through all fasta files
for fasta_file in os.listdir(fasta_dir):
    if not fasta_file.endswith("AA.fasta"):#, ".fasta", ".faa")):
        continue
    
    gene_name = os.path.splitext(fasta_file)[0]
    fasta_path = os.path.join(fasta_dir, fasta_file)
    
    # Run hmmscan (search each sequence vs Pfam)
    # need to make sure it in the environment export PATH=/programs/hmmer/bin:$PATH
    tblout_path = os.path.join(output_dir, f"{gene_name}_pfam.tblout")
    subprocess.run([
        "hmmscan",
        "--noali",                # no alignment output
        "--tblout", tblout_path,  # tabular output
        pfam_db,
        fasta_path
    ], check=True)
    
    # Parse results
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 18:
                continue
            
            seq_id = parts[2]
            pfam_id = parts[1]
            try:
                evalue = float(parts[4])
                start = int(parts[-5])
                end = int(parts[-4])
            except (ValueError, IndexError):
                print(f"Skipping malformed line in {tblout_path}: {line.strip()}")
                continue
            results.append({
                "Gene": gene_name,
                "Sequence_ID": seq_id,
                "Pfam_ID": pfam_id,
                "Evalue": evalue,
                "Start": start,
                "End": end
            })

# Convert to DataFrame and save
df = pd.DataFrame(results)
df.to_csv("pfam_summary.csv", index=False)
print(df.head())
