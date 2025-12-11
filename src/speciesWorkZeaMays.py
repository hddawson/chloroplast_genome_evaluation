import pandas as pd
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import numpy as np
import shutil
from itertools import combinations

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
    fasta_dir = "data/speciesWork/ZeaMays/proteinFastas/"
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
    out_dir = "data/speciesWork/ZeaMays/proteinsByGene/"
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
    aln_dir = "data/speciesWork/ZeaMays/alignedProteins/"
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
    axes[0, 0].set_title('Residue Polymorphisms across 18 ZeaMays accessions')
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
    plt.savefig('data/speciesWork/ZeaMays/gene_polymorphisms.png', dpi=300)
    plt.show()

    print(f"\nTotal genes: {len(results_df)}")
    print(f"Mean polymorphic sites: {results_df['polymorphic_sites'].mean():.2f}")
    print(f"Mean polymorphism rate: {results_df['poly_rate'].mean():.4f}")
    print(f"\nBy category:")
    print(results_df.groupby('category')[['polymorphic_sites', 'poly_rate']].mean())

def copy_genomes():
    for index, row in ZeaMays_data.iterrows():
        file_name = row["FileBasename"] + ".fa"
        src = os.path.join("data/genomes", file_name)
        dst_dir = "data/speciesWork/ZeaMays/genomes"
        os.makedirs(dst_dir, exist_ok=True)
        dst = os.path.join(dst_dir, file_name)
        
        assert os.path.exists(src), f"File {src} does not exist"
        shutil.copy2(src, dst)

def extract_variants_from_alignments():
    """Extract polymorphic sites from ZeaMays protein alignments"""
    aln_dir = "data/speciesWork/ZeaMays/alignedProteins/"
    
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
    """Map association P-values from modeling results to ZeaMays variants"""
    
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
    
    print(f"\n=== ZeaMays VARIANT ANALYSIS ===")
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

def identify_optimal_crosses(mapped_df, p_threshold, n_crosses=10):
    """
    Identify sample pairs that maximize differences at significant sites.
    Returns ranked cross combinations based on variant differences.
    """
    
    # Filter to significant sites only
    sig_sites = mapped_df[mapped_df['is_significant'] == True].copy()
    
    assert len(sig_sites) > 0, "No significant sites found"
    
    # Reload alignments to get per-sample genotypes at significant sites
    aln_dir = "data/speciesWork/ZeaMays/alignedProteins/"
    
    genotypes = {}  # {(gene, residue_idx): {sample: allele}}
    site_info = {}  # {(gene, residue_idx): {'p_value': x, 'R2': y}}
    samples = set()
    
    for _, row in tqdm(sig_sites.iterrows(), total=len(sig_sites), desc="Loading genotypes"):
        gene = row['gene']
        res_idx = row['residue_index']
        
        site_info[(gene, res_idx)] = {
            'p_value': row['P_res'],
            'R2_partial': row['R2_partial']
        }
        
        # Try multiple possible alignment file names
        possible_files = [f"{gene}.fasta", f"{gene}_aligned.fasta"]
        aln_path = None
        
        for fname in possible_files:
            test_path = os.path.join(aln_dir, fname)
            if os.path.exists(test_path):
                aln_path = test_path
                break
        
        if not aln_path:
            continue
            
        aln = list(SeqIO.parse(aln_path, 'fasta'))
        
        for record in aln:
            sample = record.id.split('_')[0]
            samples.add(sample)
            
            seq = str(record.seq).upper()
            non_gap_pos = [i for i, c in enumerate(seq) if c != '-']
            
            if res_idx <= len(non_gap_pos):
                allele = seq[non_gap_pos[res_idx - 1]]
                
                if (gene, res_idx) not in genotypes:
                    genotypes[(gene, res_idx)] = {}
                genotypes[(gene, res_idx)][sample] = allele
    
    samples = sorted(samples)
    print(f"Found {len(samples)} samples across {len(genotypes)} significant sites")
    
    assert len(samples) > 1, "Need at least 2 samples for crosses"
    assert len(genotypes) > 0, "No genotype data found"
    
    # Calculate pairwise differences with details
    from itertools import combinations
    
    cross_scores = []
    
    for sample1, sample2 in tqdm(list(combinations(samples, 2)), desc="Evaluating crosses"):
        differences = []
        
        for site, alleles in genotypes.items():
            if sample1 in alleles and sample2 in alleles:
                if alleles[sample1] != alleles[sample2]:
                    gene, res_idx = site
                    differences.append({
                        'gene': gene,
                        'residue': res_idx,
                        'parent1_allele': alleles[sample1],
                        'parent2_allele': alleles[sample2],
                        'p_value': site_info[site]['p_value']
                    })
        
        if len(differences) > 0:
            # Create formatted string of differences
            diff_details = '; '.join([
                f"{d['gene']}:{d['residue']}({d['parent1_allele']}/{d['parent2_allele']},P={d['p_value']:.2e})"
                for d in sorted(differences, key=lambda x: x['p_value'])
            ])
            
            cross_scores.append({
                'parent1': sample1,
                'parent2': sample2,
                'n_differences': len(differences),
                'pct_different': 100 * len(differences) / len(genotypes),
                'differences': diff_details
            })
    
    assert len(cross_scores) > 0, "No valid crosses found"
    
    cross_df = pd.DataFrame(cross_scores).sort_values('n_differences', ascending=False)
    
    print(f"\n=== TOP {n_crosses} CROSSES ===")
    for i, row in cross_df.head(n_crosses).iterrows():
        print(f"\n{row['parent1']} x {row['parent2']}: {row['n_differences']} differences ({row['pct_different']:.1f}%)")
        print(f"  {row['differences']}")
    
    return cross_df

def plot_allele_heatmap(crosses_df, mapped_df, n_crosses=10):
    """
    Create an alignment-style plot showing alleles at significant sites for top crosses
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    
    # Get unique parents from top crosses
    top_crosses = crosses_df.head(n_crosses)
    parents = set()
    for _, row in top_crosses.iterrows():
        parents.add(row['parent1'])
        parents.add(row['parent2'])
    parents = sorted(parents)
    
    # Get significant sites sorted by P-value
    sig_sites = mapped_df[mapped_df['is_significant'] == True].copy()
    sig_sites = sig_sites.sort_values('P_res')
    
    # Reload genotypes for significant sites
    aln_dir = "data/speciesWork/ZeaMays/alignedProteins/"
    genotypes = {}
    site_info = []
    
    for _, row in tqdm(sig_sites.iterrows(), total=len(sig_sites), desc="Loading genotypes"):
        gene = row['gene']
        res_idx = row['residue_index']
        
        possible_files = [f"{gene}.fasta", f"{gene}_aligned.fasta"]
        aln_path = None
        
        for fname in possible_files:
            test_path = os.path.join(aln_dir, fname)
            if os.path.exists(test_path):
                aln_path = test_path
                break
        
        if not aln_path:
            continue
        
        site_info.append({
            'gene': gene,
            'residue': res_idx,
            'p_value': row['P_res']
        })
        
        aln = list(SeqIO.parse(aln_path, 'fasta'))
        
        for record in aln:
            sample = record.id.split('_')[0]
            if sample not in parents:
                continue
                
            seq = str(record.seq).upper()
            non_gap_pos = [i for i, c in enumerate(seq) if c != '-']
            
            if res_idx <= len(non_gap_pos):
                allele = seq[non_gap_pos[res_idx - 1]]
                if sample not in genotypes:
                    genotypes[sample] = []
                genotypes[sample].append(allele)
    
    # Get organism names
    id_to_organism = ZeaMays_data.set_index('FileBasename')['Organism'].to_dict()
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(len(site_info)*0.5, len(parents)*0.4 + 2),
                                   gridspec_kw={'height_ratios': [1, len(parents)]})
    
    # Top panel: P-values
    p_values = [s['p_value'] for s in site_info]
    ax1.bar(range(len(p_values)), -np.log10(p_values), color='steelblue', alpha=0.7)
    ax1.set_ylabel('-log10(P)', fontsize=10)
    ax1.set_xlim(-0.5, len(site_info) - 0.5)
    ax1.set_xticks([])
    ax1.grid(axis='y', alpha=0.3)
    
    # Bottom panel: Allele display
    ax2.set_xlim(-0.5, len(site_info) - 0.5)
    ax2.set_ylim(-0.5, len(parents) - 0.5)
    
    # Add alleles as text
    for i, sample in enumerate(parents):
        if sample in genotypes:
            for j, allele in enumerate(genotypes[sample]):
                ax2.text(j, i, allele, ha='center', va='center', 
                        fontsize=8, family='monospace', weight='bold')
    
    # Set labels
    site_labels = [f"{s['gene']}\n{s['residue']}" for s in site_info]
    ax2.set_xticks(range(len(site_labels)))
    ax2.set_xticklabels(site_labels, fontsize=8)
    
    organism_labels = [id_to_organism.get(p, p).replace('ZeaMays ', '') for p in parents]
    ax2.set_yticks(range(len(organism_labels)))
    ax2.set_yticklabels(organism_labels, fontsize=9)
    
    ax2.set_xlabel('Gene:Residue', fontsize=10)
    ax2.grid(True, alpha=0.2)
    
    plt.suptitle('Alleles at significant sites in ZeaMays accessions', fontsize=12, y=0.995)
    plt.tight_layout()
    plt.savefig('data/speciesWork/ZeaMays/allele_alignment.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Saved allele alignment plot")


if __name__ == "__main__":

    tax_data = pd.read_csv("data/taxonomy_info.csv")

    #find all entries with ZeaMays in the organism name
    ZeaMays_data = tax_data[tax_data["Organism"].str.contains("Zea mays", case=False, na=False)].copy()

    #copy_genomes()

    #ls data/speciesWork/ZeaMays/genomes/ > data/speciesWork/ZeaMays/annotationList.txt
    #src/2_annotate.sh data/speciesWork/ZeaMays/annotationList.txt
    #mkdir data/speciesWork/ZeaMays/proteinFastas
    #cp  data/speciesWork/ZeaMays/annotationResults/*/*/*Protein.fasta  data/speciesWork/ZeaMays/proteinFastas/
    
    #src/alignerIntraspecific.sh (make sure to modify the script to point to the correct directories)
    write_proteins_by_gene()
    #count_polymorphisms()

    variants_df = extract_variants_from_alignments()
    variants_df.to_csv('data/speciesWork/ZeaMays/variants.csv', index=False)

    # Map association scores
    mapped_df, p_threshold = map_association_scores(variants_df, results_file="data/speciesWork/ZeaMays/ZeaMays_reference_mapped_results.csv")
    mapped_df.to_csv('data/speciesWork/ZeaMays/variants_with_scores.csv', index=False)

    # Summarize significant variants
    sig_variants = summarize_significant_variants(mapped_df, p_threshold)

    crosses_df = identify_optimal_crosses(mapped_df, p_threshold, n_crosses=10)

    id_to_organism = ZeaMays_data.set_index('FileBasename')['Organism'].to_dict()
    print(id_to_organism)

    crosses_df['parent1_organism'] = crosses_df['parent1'].map(id_to_organism)
    crosses_df['parent2_organism'] = crosses_df['parent2'].map(id_to_organism)

    crosses_df.to_csv('data/speciesWork/ZeaMays/optimal_crosses.csv', index=False)

    # Plot cross comparisons
    plot_allele_heatmap(crosses_df, mapped_df, n_crosses=10)