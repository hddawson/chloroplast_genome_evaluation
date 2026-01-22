"""
Map GWAS results to Arabidopsis chloroplast genomic coordinates.

Usage:
1. Run export_gwas_results.R first to create results/gwas_results_clean.csv
2. Run this script: python map_to_arabidopsis.py
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
from pathlib import Path

# ---- CONFIG ----
ARABIDOPSIS_ID = "AP000423.1"
GBF_FILE = "data/gbfs/AP0004231fa.gbf"  # adjust if needed
GWAS_FILE = "results/gwas_results_clean.csv"
ALN_DIR = Path("data/tmp/alignedGenes")

# ---- PARSE GENBANK FILE ----
def parse_gbf(gbf_path):
    """Extract CDS gene names and coordinates from GenBank file."""
    gene_coords = {}
    
    for record in SeqIO.parse(gbf_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", [None])[0]
                if gene_name is None:
                    gene_name = feature.qualifiers.get("product", ["unknown"])[0]
                
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                
                gene_coords[gene_name] = {
                    "start": start,
                    "end": end,
                    "strand": strand
                }
    
    return gene_coords

# ---- MAP ALIGNMENT POSITIONS ----
def get_alignment_mapping(aln_file, target_id):
    """
    Map alignment positions to ungapped positions for target sequence.
    Returns dict: {aligned_pos (1-based): ungapped_pos (1-based)}
    """
    from Bio import AlignIO
    
    try:
        alignment = AlignIO.read(aln_file, "fasta")
    except:
        return None
    
    # Find target sequence
    target_seq = None
    for record in alignment:
        clean_id = record.id.split("|")[0]
        if clean_id == target_id:
            target_seq = str(record.seq)
            break
    
    if target_seq is None:
        return None
    
    # Build mapping
    mapping = {}
    ungapped_pos = 0
    
    for i, char in enumerate(target_seq):
        aln_pos = i + 1  # 1-based
        if char != "-":
            ungapped_pos += 1
            mapping[aln_pos] = ungapped_pos
        else:
            mapping[aln_pos] = None  # gap in target
    
    return mapping

# ---- MAIN ----
def main():
    # Load GWAS results
    gwas = pd.read_csv(GWAS_FILE)
    assert len(gwas) > 0, "No GWAS results found"
    print(f"Loaded {len(gwas)} GWAS positions")
    
    # Parse GenBank
    gene_coords = parse_gbf(GBF_FILE)
    print(f"Parsed {len(gene_coords)} genes from GenBank")
    print(f"Genes: {list(gene_coords.keys())}")
    
    # Get alignments and build mappings
    genes = gwas["Gene"].unique()
    all_mappings = {}
    
    for gene in genes:
        aln_file = ALN_DIR / f"{gene}_AA_aligned.fasta"
        if not aln_file.exists():
            print(f"Warning: alignment not found for {gene}")
            continue
        
        mapping = get_alignment_mapping(aln_file, ARABIDOPSIS_ID)
        if mapping is None:
            print(f"Warning: {ARABIDOPSIS_ID} not found in {gene} alignment")
            continue
        
        all_mappings[gene] = mapping
        print(f"  {gene}: mapped {sum(v is not None for v in mapping.values())} positions")
    
    # Map GWAS positions
    results = []
    
    for _, row in gwas.iterrows():
        gene = row["Gene"]
        aln_pos = int(row["Aligned_Position"])
        
        # Get AA position in Arabidopsis
        if gene not in all_mappings:
            continue
        
        at_aa_pos = all_mappings[gene].get(aln_pos)
        if at_aa_pos is None:
            continue  # gap in Arabidopsis
        
        # Get genomic coordinates
        if gene not in gene_coords:
            print(f"Warning: {gene} not in GenBank annotations")
            continue
        
        coords = gene_coords[gene]
        
        if coords["strand"] == "+":
            genomic_pos = coords["start"] + (at_aa_pos - 1) * 3
        else:
            genomic_pos = coords["end"] - (at_aa_pos - 1) * 3 - 2
        
        results.append({
            "Gene": gene,
            "Aligned_Position": aln_pos,
            "At_AA_Position": at_aa_pos,
            "Genomic_Position": genomic_pos,
            "Strand": coords["strand"],
            "P_aa_with_pcs": row["P_aa_with_pcs"],
            "P_aa_only": row["P_aa_only"],
            "gwas_hit": row["gwas_hit"],
            "neg_log10_p": row["neg_log10_p"]
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("Genomic_Position")
    
    out_path = "results/gwas_arabidopsis_mapped.csv"
    results_df.to_csv(out_path, index=False)
    
    print(f"\nMapped {len(results_df)} positions to Arabidopsis genome")
    print(f"GWAS hits: {results_df['gwas_hit'].sum()}")
    print(f"Saved to: {out_path}")
    
    # Quick stats
    print(f"\nGenomic range: {results_df['Genomic_Position'].min()} - {results_df['Genomic_Position'].max()}")
    
    return results_df

if __name__ == "__main__":
    main()