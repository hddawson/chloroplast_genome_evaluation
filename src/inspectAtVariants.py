#!/usr/bin/env python3
"""
Convert GBF annotations, genome FASTA, and variant calls to IGV-compatible formats.
Outputs: GFF3 annotation file, FASTA genome, and VCF variant file.
"""

import pandas as pd
from Bio import SeqIO
import sys
import os

# Configuration
GBF_FILE = 'data/speciesWork/At/gbfs/NC0009321fa.gbf'
VARIANT_FILE = 'data/chloroplast_mitochondria_variant_calls.txt'
OUTPUT_DIR = 'data/speciesWork/At/igv_files'

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Extract genome FASTA from GBF
print("Extracting genome sequence...")
fasta_output = os.path.join(OUTPUT_DIR, 'NC0009321.fasta')
record = SeqIO.read(GBF_FILE, 'genbank')
assert len(record.seq) > 0, "Empty sequence in GBF file"

with open(fasta_output, 'w') as f:
    f.write(f">{record.id}\n")
    f.write(str(record.seq) + "\n")
print(f"Wrote FASTA: {fasta_output}")

# 2. Convert annotations to GFF3
# 2. Convert annotations to GFF3
print("Converting annotations to GFF3...")
gff_output = os.path.join(OUTPUT_DIR, 'NC0009321.gff3')

with open(gff_output, 'w') as gff:
    gff.write("##gff-version 3\n")

    for feature in record.features:
        if feature.type == 'gene':
            gene_name = feature.qualifiers.get('gene', [''])[0] or f"gene_{feature.location.start}"
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'
            gff.write(f"{record.id}\tGenBank\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_name};Name={gene_name}\n")

        elif feature.type == 'CDS':
            gene_name = feature.qualifiers.get('gene', [''])[0] or f"gene_{feature.location.start}"
            product = feature.qualifiers.get('product', [''])[0]
            translation = feature.qualifiers.get('translation', [''])[0]  # <— include AA sequence if available

            # IGV prefers CDS to have a Parent (usually an mRNA)
            parent_id = f"{gene_name}_mRNA"
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'

            attrs = [f"ID={gene_name}_CDS", f"Parent={parent_id}", f"Name={gene_name}"]
            if product:
                attrs.append(f"product={product}")
            if translation:
                attrs.append(f"translation={translation}")

            gff.write(f"{record.id}\tGenBank\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={parent_id};Parent={gene_name}\n")
            gff.write(f"{record.id}\tGenBank\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{';'.join(attrs)}\n")

print(f"Wrote GFF3: {gff_output}")

# 3. Convert variants to VCF
print("Converting variants to VCF...")
variant_df = pd.read_csv(VARIANT_FILE, sep='\t')
c_variant_df = variant_df[variant_df['CHROM'] == 'chloroplast'].copy()
assert len(c_variant_df) > 0, "No chloroplast variants found"

vcf_output = os.path.join(OUTPUT_DIR, 'chloroplast_variants.vcf')
sample_cols = c_variant_df.columns[2:].tolist()

with open(vcf_output, 'w') as vcf:
    # Header
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write(f"##reference={record.id}\n")
    vcf.write(f"##contig=<ID={record.id},length={len(record.seq)}>\n")
    vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    
    # Column header
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_cols) + "\n")
    
    # Variant rows
    for idx, row in c_variant_df.iterrows():
        pos = row['POS']
        alleles = [row[col] for col in sample_cols]
        
        # Determine REF and ALT
        ref = alleles[0]  # First sample as reference
        alt_alleles = list(set(alleles) - {ref})
        
        if not alt_alleles:
            continue
            
        alt = ','.join(alt_alleles)
        
        # Genotypes
        allele_map = {ref: '0'}
        for i, a in enumerate(alt_alleles, 1):
            allele_map[a] = str(i)
        
        genotypes = [allele_map.get(a, '.') for a in alleles]
        gt_string = '\t'.join(genotypes)
        
        vcf.write(f"{record.id}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt_string}\n")

print(f"Wrote VCF: {vcf_output}")

# 4. Summary statistics
print("\n=== Summary ===")
print(f"Genome length: {len(record.seq)} bp")
print(f"Total features: {len(record.features)}")
print(f"CDS features: {sum(1 for f in record.features if f.type == 'CDS')}")
print(f"Total variants: {len(c_variant_df)}")

# Variants in CDS
variants_in_cds = 0
for _, variant in c_variant_df.iterrows():
    pos = variant['POS']
    for feature in record.features:
        if feature.type == 'CDS':
            for part in feature.location.parts:
                if int(part.start) <= pos <= int(part.end):
                    variants_in_cds += 1
                    break

print(f"Variants in CDS: {variants_in_cds} ({100*variants_in_cds/len(c_variant_df):.1f}%)")

print("\n=== IGV Files Ready ===")
print(f"Load these files into IGV:")
print(f"1. Genome: {fasta_output}")
print(f"2. Annotations: {gff_output}")
print(f"3. Variants: {vcf_output}")