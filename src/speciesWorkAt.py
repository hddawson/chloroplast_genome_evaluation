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
import shutil

import os
from Bio import SeqIO
from collections import defaultdict

def write_proteins_by_gene():
    """
    Reads protein fasta files from a specified directory, groups sequences by protein name,
    and writes them into separate fasta files for each protein.
    """
    # Dictionary to store sequences by protein name
    proteins = defaultdict(list)

    # Read all protein fasta files
    fasta_dir = "data/speciesWork/At/proteinFastas/"
    for fasta_file in os.listdir(fasta_dir):
        if not fasta_file.endswith(".fasta"):
            continue
        
        # Extract sample name (e.g., AP0004231fa from AP0004231fa.Protein.fasta)
        sample_name = fasta_file.replace("fa.Protein.fasta", "")
        
        # Parse sequences
        for record in SeqIO.parse(os.path.join(fasta_dir, fasta_file), "fasta"):
            # Extract protein name (e.g., rps12 from rps12_join{...})
            protein_name = record.id.split("_")[0]
            
            # Create new ID with sample name
            record.id = f"{sample_name}_{record.id}"
            record.description = ""
            
            proteins[protein_name].append(record)

    # Write protein-specific fastas
    out_dir = "data/speciesWork/At/proteinsByGene/"
    os.makedirs(out_dir, exist_ok=True)

    for protein_name, records in proteins.items():
        out_file = os.path.join(out_dir, f"{protein_name}.fasta")
        SeqIO.write(records, out_file, "fasta")
        print(f"Wrote {len(records)} sequences to {out_file}")


def count_polymorphisms_per_gene(aln_file):
    """Count polymorphic sites in an alignment"""
    aln = list(SeqIO.parse(aln_file, 'fasta'))
    aln_len = len(aln[0].seq)
    
    assert all(len(r.seq) == aln_len for r in aln), f'Sequences have different lengths in {aln_file}'
    
    aln_matrix = np.array([list(str(r.seq).upper()) for r in aln])
    
    # Count sites with >1 non-gap allele
    n_polymorphic = 0
    for i in range(aln_len):
        site = aln_matrix[:, i]
        non_gap = site[site != '-']
        if len(np.unique(non_gap)) > 1:
            n_polymorphic += 1
    
    return n_polymorphic, aln_len

def get_gene_category(gene_name):
    """Assign gene to functional category based on prefix"""
    prefixes = ['psb', 'psa', 'pet', 'atp', 'ndh', 'rbc']
    for prefix in prefixes:
        if gene_name.startswith(prefix):
            return prefix
    return 'other'


def count_polymorphisms():

    # Process all alignments
    aln_dir = "data/speciesWork/At/alignedProteins/"
    results = []

    for aln_file in tqdm(sorted(os.listdir(aln_dir))):
        if not aln_file.endswith(".fasta"):
            continue

        if aln_file.startswith("ycf"):
            continue
        
        gene_name = aln_file.replace(".fasta", "")
        n_poly, aln_len = count_polymorphisms_per_gene(os.path.join(aln_dir, aln_file))
        poly_rate = n_poly / aln_len if aln_len > 0 else 0
        category = get_gene_category(gene_name)
        results.append({'gene': gene_name, 'polymorphic_sites': n_poly, 
                        'length': aln_len, 'poly_rate': poly_rate, 'category': category})


    results_df = pd.DataFrame(results)

    # Color mapping
    colors = {'psb': 'blue', 'psa': 'green', 'pet': 'red', 
            'atp': 'orange', 'ndh': 'purple', 'rbc': 'brown', 'other': 'gray'}

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Polymorphic sites per gene
    for category in results_df['category'].unique():
        mask = results_df['category'] == category
        axes[0, 0].bar(np.where(mask)[0], results_df[mask]['polymorphic_sites'], 
                    color=colors[category], label=category, alpha=0.7)
    axes[0, 0].set_xlabel('Gene')
    axes[0, 0].set_ylabel('Polymorphic residues per gene')
    axes[0, 0].set_title('Residue Polymorphisms across 11 A. thaliana accessions')
    axes[0, 0].legend()

    # Polymorphism rate by category
    category_means = results_df.groupby('category')['poly_rate'].mean()
    axes[0, 1].bar(range(len(category_means)), category_means.values, 
                color=[colors[c] for c in category_means.index])
    axes[0, 1].set_xticks(range(len(category_means)))
    axes[0, 1].set_xticklabels(category_means.index, rotation=45)
    axes[0, 1].set_ylabel('Mean polymorphism rate')
    axes[0, 1].set_title('Polymorphism rate by category')

    # Distribution by category
    for category in results_df['category'].unique():
        data = results_df[results_df['category'] == category]['polymorphic_sites']
        axes[1, 0].hist(data, bins=15, alpha=0.5, label=category, color=colors[category])
    axes[1, 0].set_xlabel('Polymorphic sites')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of polymorphic sites')
    axes[1, 0].legend()

    # Length vs polymorphisms
    for category in results_df['category'].unique():
        mask = results_df['category'] == category
        axes[1, 1].scatter(results_df[mask]['length'], results_df[mask]['polymorphic_sites'], 
                        alpha=0.6, label=category, color=colors[category])
    axes[1, 1].set_xlabel('Gene length (aa)')
    axes[1, 1].set_ylabel('Polymorphic sites')
    axes[1, 1].set_title('Length vs polymorphisms')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig('data/speciesWork/At/gene_polymorphisms.png', dpi=300)
    plt.show()

    print(f"\nTotal genes: {len(results_df)}")
    print(f"Mean polymorphic sites: {results_df['polymorphic_sites'].mean():.2f}")
    print(f"Mean polymorphism rate: {results_df['poly_rate'].mean():.4f}")
    print(f"\nBy category:")
    print(results_df.groupby('category')[['polymorphic_sites', 'poly_rate']].mean())

if __name__ == "__main__":
    #write_proteins_by_gene()
    count_polymorphisms()