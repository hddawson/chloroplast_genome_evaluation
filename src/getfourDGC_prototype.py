import pandas as pd
import os
import glob
from tqdm import tqdm
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import Counter

# --- config ---
indir = "data/tmp"   # where *_pal2nal.fa live
genetic_code = 11
fourfold_prefixes = {"GC","CG","CC","GG","CT","TC","GT","AC"}

def is_unambiguous_codon(c):
    return set(c.upper()) <= set("ACGT")

def process_gene(in_fa, gene):
    """Process one gene alignment and return GC summary DataFrame"""
    records = list(SeqIO.parse(in_fa, "fasta"))
    L = len(records[0].seq)
    assert all(len(r.seq)==L for r in records), "All sequences must have same length"
    assert L % 3 == 0, "Alignment length must be a multiple of 3"

    stop_codons = set(CodonTable.unambiguous_dna_by_id[genetic_code].stop_codons)
    ncod = L // 3

    # consensus valid codons
    valid_k = []
    for k in range(ncod):
        codons_here = []
        for rec in records:
            codon = str(rec.seq[3*k:3*k+3]).upper()
            if "-" in codon or not is_unambiguous_codon(codon) or codon in stop_codons:
                continue
            codons_here.append(codon)
        if not codons_here:
            continue

        prefix_counts = Counter(c[:2] for c in codons_here)
        prefix, count = prefix_counts.most_common(1)[0]
        freq = count / sum(prefix_counts.values())
        if prefix in fourfold_prefixes and freq >= 0.99:
            valid_k.append((k, prefix))

    # mask + count GC
    gc_summary = []
    for rec in records:
        s = list(str(rec.seq).upper())
        counts = Counter()
        for k in range(ncod):
            codon = "".join(s[3*k:3*k+3])
            if "-" in codon or not is_unambiguous_codon(codon) or codon in stop_codons:
                continue
            valid_prefix = next((prefix for vk, prefix in valid_k if k == vk), None)
            if valid_prefix and codon[:2] == valid_prefix:
                base3 = codon[2]
                if base3 in "ACGT":
                    counts[base3] += 1

        total = counts["A"] + counts["T"] + counts["G"] + counts["C"]
        gc_summary.append({
            "taxon": rec.id,
            "ID": rec.id.split("|")[0],
            f"{gene}_nuc_G": counts["G"],
            f"{gene}_nuc_C": counts["C"],
            f"{gene}_nuc_A": counts["A"],
            f"{gene}_nuc_T": counts["T"],
            f"{gene}_GC_content": (counts["G"] + counts["C"]) / total if total > 0 else None
        })

    return pd.DataFrame(gc_summary)

def main():
    files = glob.glob(os.path.join(indir, "**/*_pal2nal.fa"), recursive=True)
    all_dfs = []

    for f in tqdm(files, desc="Processing genes"):
        gene = os.path.basename(f).split("_")[0]  # assumes format gene_pal2nal.fa
        df = process_gene(f, gene)
        all_dfs.append(df)

    # Merge all gene dfs
    merged = all_dfs[0]
    for df in all_dfs[1:]:
        merged = pd.merge(merged, df, on=["taxon","ID"], how="outer")

    merged.to_csv("data/all_genes_gc_summary.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
