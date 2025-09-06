import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import re 
from tqdm import tqdm
import seaborn as sns
from Bio import SeqIO
from Bio import SeqUtils
import numpy as np
import glob
from Bio.Data import CodonTable
from collections import defaultdict, Counter
from pandarallel import pandarallel
from scipy.stats import pearsonr


# do this subcommand 
#/programs/pal2nal/pal2nal.pl data/tmp/alignedGenes/rbcL_AA_aligned.fasta data/tmp/genesToAlign/rbcL_CDS.fasta -codontable 11 -output fasta > data/rbcL_pal2nal.fa

# --- config ---
in_fa = "data/tmp/rbcL_aln/rbcL_pal2nal.fa"     # codon-aware DNA alignment
genetic_code = 11  # 1 for Standard, 11 for Bacterial/Archaeal/Plant Plastid
mode = "strict"    
mask_char = "N"

fourfold_prefixes = { "GC","CG","CC","GG","CT","TC","GT","AC" }

records = list(SeqIO.parse(in_fa, "fasta"))
L = len(records[0].seq)
assert all(len(r.seq)==L for r in records), "All sequences must have same length"
assert L % 3 == 0, "Alignment length must be a multiple of 3"

def is_unambiguous_codon(c):
    return set(c.upper()) <= set("ACGT")

# Precompute which codon indices are 4D per taxon
per_taxon_4d = [set() for _ in records]  # set of codon indices (0..ncodons-1)
stop_codons = set(CodonTable.unambiguous_dna_by_id[genetic_code].stop_codons)

ncod = L // 3

for i, rec in tqdm(enumerate(records)):
    s = str(rec.seq).upper()
    for k in range(ncod):
        codon = s[3*k:3*k+3]
        if "-" in codon or not is_unambiguous_codon(codon):
            continue
        if codon in stop_codons:
            continue
        if codon[:2] in fourfold_prefixes:
            per_taxon_4d[i].add(k)

valid_k = []
codon_distributions = {}  # store codon counts per site
prefix_distributions = {}

for k in tqdm(range(ncod)):
    codons_here = []
    for rec in records:
        codon = str(rec.seq[3*k:3*k+3]).upper()
        if ("-" in codon) or (not is_unambiguous_codon(codon)) or (codon in stop_codons):
            continue
        codons_here.append(codon)
    
    if not codons_here:
        continue
    
    # Count codons and prefixes
    codon_counts = Counter(codons_here)
    prefix_counts = Counter(c[:2] for c in codons_here)
    
    # Majority prefix
    prefix, count = prefix_counts.most_common(1)[0]
    freq = count / sum(prefix_counts.values())
    
    codon_distributions[k] = codon_counts
    prefix_distributions[k] = prefix_counts
    
    #95% consensus gives 207 sites
    if prefix in fourfold_prefixes and freq >= 0.99:
        valid_k.append((k, prefix))

print(f"Relaxed 4-fold codons (95% consensus): {len(valid_k)}")

# Example: show distribution at one site
site = valid_k[0]

masked_records = []
gc_summary = []

for rec in tqdm(records):
    s = list(str(rec.seq).upper())
    counts = Counter()

    for k in range(ncod):  # iterate over codons
        codon = "".join(s[3*k:3*k+3])

        # skip ambiguous/gap/stop codons
        if "-" in codon or not is_unambiguous_codon(codon) or codon in stop_codons:
            s[3*k:3*k+3] = ["N"]*3
            continue

        # check if codon is valid 4D site
        valid_prefix = next((prefix for vk, prefix in valid_k if k == vk), None)

        if valid_prefix and codon[:2] == valid_prefix:
            # keep only 3rd base
            base3 = codon[2]
            s[3*k:3*k+2] = ["N","N"]
            s[3*k+2] = base3
            if base3 in "ACGT":
                counts[base3] += 1
        else:
            # mask non-4D or deviant codon
            s[3*k:3*k+3] = ["N"]*3

    # update sequence with masked bases
    rec.seq = rec.seq.__class__("".join(s))
    masked_records.append(rec)

    # compute GC summaries
    total = counts["A"] + counts["T"] + counts["G"] + counts["C"]
    gc_summary.append({
        "taxon": rec.id,
        "rbcL_nuc_G": counts["G"],
        "rbcL_nuc_C": counts["C"],
        "rbcL_nuc_A": counts["A"],
        "rbcL_nuc_T": counts["T"],
        "GC_content": (counts["G"] + counts["C"]) / total if total > 0 else None
    })

# Save masked alignment
#SeqIO.write(masked_records, "data/tmp/rbcL_aln/aligned_cds.4D.masked.fasta", "fasta")

# Save GC summary
df_gc = pd.DataFrame(gc_summary)

df_gc["ID"] = df_gc["taxon"].apply(lambda x: x.split("|")[0])
df_gc.to_csv("data/rbcL_gc_summary_per_taxon.tsv", sep="\t", index=False)