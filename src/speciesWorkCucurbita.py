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
    fasta_dir = "data/speciesWork/Cucurbita/proteinFastas/"
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
    out_dir = "data/speciesWork/Cucurbita/proteinsByGene/"
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
    aln_dir = "data/speciesWork/Cucurbita/alignedProteins/"
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
    axes[0, 0].set_title('Residue Polymorphisms across 18 Cucurbita accessions')
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
    plt.savefig('data/speciesWork/Cucurbita/gene_polymorphisms.png', dpi=300)
    plt.show()

    print(f"\nTotal genes: {len(results_df)}")
    print(f"Mean polymorphic sites: {results_df['polymorphic_sites'].mean():.2f}")
    print(f"Mean polymorphism rate: {results_df['poly_rate'].mean():.4f}")
    print(f"\nBy category:")
    print(results_df.groupby('category')[['polymorphic_sites', 'poly_rate']].mean())

def copy_genomes():
    for index, row in cucurbita_data.iterrows():
        file_name = row["FileBasename"] + ".fa"
        src = os.path.join("data/genomes", file_name)
        dst_dir = "data/speciesWork/Cucurbita/genomes"
        os.makedirs(dst_dir, exist_ok=True)
        dst = os.path.join(dst_dir, file_name)
        
        assert os.path.exists(src), f"File {src} does not exist"
        shutil.copy2(src, dst)

def extract_variants_from_alignments():
    """Extract polymorphic sites from Cucurbita protein alignments"""
    aln_dir = "data/speciesWork/Cucurbita/alignedProteins/"
    
    variants = []
    
    for aln_file in tqdm(sorted(os.listdir(aln_dir))):
        if not aln_file.endswith(".fasta") or aln_file.startswith("ycf"):
            continue
        
        gene_name = aln_file.replace("_aligned.fasta", "")
        aln = list(SeqIO.parse(os.path.join(aln_dir, aln_file), 'fasta'))
        aln_len = len(aln[0].seq)
        
        assert all(len(r.seq) == aln_len for r in aln), f'Unequal lengths in {aln_file}'
        
        aln_matrix = np.array([list(str(r.seq).upper()) for r in aln])
        
        # Find polymorphic sites
        for pos in range(aln_len):
            site = aln_matrix[:, pos]
            non_gap = site[site != '-']
            
            if len(np.unique(non_gap)) > 1:
                # Calculate residue index (non-gap position)
                residue_idx = np.sum(site[0] != '-' for _ in range(pos) if aln_matrix[0, _] != '-') + 1
                
                variants.append({
                    'gene': gene_name,
                    'aligned_position': pos + 1,
                    'residue_index': residue_idx,
                    'n_alleles': len(np.unique(non_gap)),
                    'alleles': ','.join(np.unique(non_gap))
                })
    
    return pd.DataFrame(variants)


def map_association_scores(variants_df, results_file='results/reference_mapped_results.csv'):
    """Map association P-values from modeling results to Cucurbita variants"""
    
    assert os.path.exists(results_file), f"Results file {results_file} not found"
    
    results_df = pd.read_csv(results_file)
    
    # Merge on gene and residue index
    mapped = variants_df.merge(
        results_df[['Gene', 'Residue_Index', 'P_res', 'R2_partial']],
        left_on=['gene', 'residue_index'],
        right_on=['Gene', 'Residue_Index'],
        how='left'
    )
    
    # Calculate significance threshold
    p_threshold = results_df['P_res'].quantile(0.1)
    mapped['is_significant'] = mapped['P_res'] <= p_threshold
    
    return mapped, p_threshold


def summarize_significant_variants(mapped_df, p_threshold):
    """Summarize variants at significant association sites"""
    
    sig = mapped_df[mapped_df['is_significant'] == True].copy()
    
    print(f"\n=== CUCURBITA VARIANT ANALYSIS ===")
    print(f"Total polymorphic sites: {len(mapped_df)}")
    print(f"Sites with association data: {mapped_df['P_res'].notna().sum()}")
    print(f"Significant sites (P_res < {p_threshold:.2e}): {len(sig)}")
    
    if len(sig) > 0:
        print(f"\n=== TOP SIGNIFICANT VARIANTS ===")
        top = sig.nsmallest(20, 'P_res')[['gene', 'residue_index', 'P_res', 'R2_partial', 'n_alleles', 'alleles']]
        print(top.to_string(index=False))
        
        print(f"\n=== VARIANTS BY GENE ===")
        gene_counts = sig.groupby('gene').size().sort_values(ascending=False)
        print(gene_counts)
    
    return sig



if __name__ == "__main__":

    tax_data = pd.read_csv("data/taxonomy_info.csv")

    #find all entries with Cucurbita in the organism name
    cucurbita_data = tax_data[tax_data["Organism"].str.contains("Cucurbita")]

    #copy_genomes()

    #ls data/speciesWork/Cucurbita/genomes/ > data/speciesWork/Cucurbita/annotationList.txt
    #src/2_annotate.sh data/speciesWork/Cucurbita/annotationList.txt
    #mkdir data/speciesWork/Cucurbita/proteinFastas
    #cp  data/speciesWork/Cucurbita/annotationResults/*/*/*Protein.fasta  data/speciesWork/Cucurbita/proteinFastas/
    #write_proteins_by_gene()
    #count_polymorphisms()

    variants_df = extract_variants_from_alignments()
    variants_df.to_csv('data/speciesWork/Cucurbita/variants.csv', index=False)

    # Map association scores
    mapped_df, p_threshold = map_association_scores(variants_df, results_file="data/speciesWork/Cucurbita/cucurbitaPepo_reference_mapped_results.csv")
    mapped_df.to_csv('data/speciesWork/Cucurbita/variants_with_scores.csv', index=False)

    # Summarize significant variants
    sig_variants = summarize_significant_variants(mapped_df, p_threshold)