"""
Parse all GBF files to extract CDS coordinates.
Output: data/gbf_gene_coords.csv
"""

from Bio import SeqIO
from pathlib import Path
import pandas as pd
from tqdm import tqdm

GBF_DIR = Path("data/gbfs")

records = []

for gbf_file in tqdm(GBF_DIR.glob("*.gbf")):
    # Extract ID from filename (e.g., AP0004231fa.gbf -> AP000423.1)
    stem = gbf_file.stem.replace("fa", "")
    # Insert dot before last digit
    genome_id = stem[:-1] + "." + stem[-1]
    
    try:
        for record in SeqIO.parse(gbf_file, "genbank"):
            genome_size = len(record.seq)
            
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                
                gene = feature.qualifiers.get("gene", [None])[0]
                if gene is None:
                    continue
                
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                
                records.append({
                    "genome_id": genome_id,
                    "genome_size": genome_size,
                    "gene": gene,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
    except Exception as e:
        print(f"Error parsing {gbf_file}: {e}")

df = pd.DataFrame(records)
df.to_csv("data/gbf_gene_coords.csv", index=False)

print(f"Parsed {df['genome_id'].nunique()} genomes, {len(df)} CDS features")
print(f"Genes: {df['gene'].unique().tolist()}")