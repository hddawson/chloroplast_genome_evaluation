#!/usr/bin/env python3
"""
Analyze chloroplast variants in context of protein residues and statistical results.
Identifies variants in CDS, predicts synonymous/non-synonymous changes, and cross-references
with significant association sites from R analysis.
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

# Configuration
GBF_FILE = 'data/speciesWork/At/gbfs/NC0009321fa.gbf'
VARIANT_FILE = 'data/chloroplast_mitochondria_variant_calls.txt'
RESULTS_FILE = 'results/reference_mapped_results.csv'
OUTPUT_FILE = 'results/variant_analysis.csv'

# Load data
print("Loading variant calls...")
variant_df = pd.read_csv(VARIANT_FILE, sep='\t')
c_variants = variant_df[variant_df['CHROM'] == 'chloroplast'].copy()
assert len(c_variants) > 0, "No chloroplast variants found"

print("Loading GBF annotations...")
record = SeqIO.read(GBF_FILE, 'genbank')
assert len(record.seq) > 0, "Empty sequence in GBF file"

print("Loading association results...")
results_df = pd.read_csv(RESULTS_FILE)
assert len(results_df) > 0, "Empty results file"

# Get lower decile threshold for P_res
p_res_threshold = results_df['P_res'].quantile(0.1)
print(f"Lower quartile P_res threshold: {p_res_threshold:.2e}")

# Build CDS lookup table
print("\nBuilding CDS annotation table...")
cds_annotations = []

for feature in record.features:
    if feature.type != 'CDS':
        continue
    
    gene_name = feature.qualifiers.get('gene', [''])[0]
    if not gene_name:
        continue
    
    translation = feature.qualifiers.get('translation', [''])[0]
    strand = feature.location.strand
    
    # Handle multi-part features (e.g., genes split by introns)
    for part in feature.location.parts:
        start = int(part.start) + 1  # Convert to 1-based
        end = int(part.end)
        
        cds_annotations.append({
            'gene': gene_name,
            'start': start,
            'end': end,
            'strand': strand,
            'translation': translation
        })

cds_df = pd.DataFrame(cds_annotations)
print(f"Found {len(cds_df)} CDS regions across {cds_df['gene'].nunique()} genes")

# Analyze each variant
print("\nAnalyzing variants...")
variant_results = []

for idx, variant in c_variants.iterrows():
    pos = variant['POS']
    
    result = {
        'POS': pos,
        'in_CDS': False,
        'gene': None,
        'residue_index': None,
        'ref_codon': None,
        'ref_aa': None,
        'is_nonsynonymous': False,
        'num_nonsynonymous': 0,
        'total_alleles': 0,
        'in_significant_site': False,
        'P_res': None
    }
    
    # Check if variant is in any CDS
    cds_match = cds_df[(cds_df['start'] <= pos) & (cds_df['end'] >= pos)]
    
    if len(cds_match) > 0:
        result['in_CDS'] = True
        cds = cds_match.iloc[0]
        result['gene'] = cds['gene']
        
        # Calculate residue index within CDS
        if cds['strand'] == 1:
            cds_pos = pos - cds['start']
        else:
            cds_pos = cds['end'] - pos
        
        codon_index = cds_pos // 3
        codon_pos = cds_pos % 3
        result['residue_index'] = codon_index + 1
        
        # Extract reference codon
        if cds['strand'] == 1:
            codon_start = cds['start'] - 1 + (codon_index * 3)
            ref_codon = str(record.seq[codon_start:codon_start + 3])
        else:
            codon_start = cds['end'] - (codon_index * 3) - 3
            ref_codon = str(record.seq[codon_start:codon_start + 3].reverse_complement())
        
        result['ref_codon'] = ref_codon
        
        if len(ref_codon) == 3:
            result['ref_aa'] = str(Seq(ref_codon).translate())
            
            # Check all alternate alleles for non-synonymous changes
            sample_cols = c_variants.columns[2:].tolist()
            alleles = [variant[col] for col in sample_cols]
            ref_allele = alleles[0]
            
            result['total_alleles'] = len(set(alleles))
            
            for allele in set(alleles):
                if allele == ref_allele:
                    continue
                
                # Create alternate codon
                alt_codon = list(ref_codon)
                alt_codon[codon_pos] = allele
                alt_codon = ''.join(alt_codon)
                
                if len(alt_codon) == 3:
                    alt_aa = str(Seq(alt_codon).translate(table=11))
                    
                    if alt_aa != result['ref_aa']:
                        result['is_nonsynonymous'] = True
                        result['num_nonsynonymous'] += 1
        
        # Check if this position is significant in results
        gene_results = results_df[
            (results_df['Gene'] == result['gene']) & 
            (results_df['Residue_Index'] == result['residue_index'])
        ]
        
        if len(gene_results) > 0:
            p_res = gene_results.iloc[0]['P_res']
            result['P_res'] = p_res
            
            if p_res <= p_res_threshold:
                result['in_significant_site'] = True
    
    variant_results.append(result)

# Convert to DataFrame and save
variant_results_df = pd.DataFrame(variant_results)
variant_results_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nSaved detailed results to: {OUTPUT_FILE}")

# Summary statistics
print("\n=== VARIANT ANALYSIS SUMMARY ===")
print(f"1. Total variants: {len(variant_results_df)}")
print(f"2. Variants in CDS: {variant_results_df['in_CDS'].sum()}")
print(f"3. Non-synonymous variants: {variant_results_df['is_nonsynonymous'].sum()}")

significant_nonsyn = variant_results_df[
    variant_results_df['in_significant_site'] & 
    variant_results_df['is_nonsynonymous']
]
print(f"4. Non-synonymous variants at significant sites (lower quartile P_res): {len(significant_nonsyn)}")

# Additional breakdowns
if len(significant_nonsyn) > 0:
    print("\n=== SIGNIFICANT NON-SYNONYMOUS VARIANTS ===")
    for _, row in significant_nonsyn.iterrows():
        print(f"  {row['gene']} residue {row['residue_index']} (pos {row['POS']}): "
              f"P_res = {row['P_res']:.2e}")

# Gene-level summary
print("\n=== VARIANTS BY GENE ===")
gene_summary = variant_results_df[variant_results_df['in_CDS']].groupby('gene').agg({
    'POS': 'count',
    'is_nonsynonymous': 'sum',
    'in_significant_site': 'sum'
}).rename(columns={
    'POS': 'total_variants',
    'is_nonsynonymous': 'nonsynonymous',
    'in_significant_site': 'in_significant_sites'
})

print(gene_summary.sort_values('nonsynonymous', ascending=False))
