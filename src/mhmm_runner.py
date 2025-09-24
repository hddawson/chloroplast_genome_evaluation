#!/usr/bin/env python3
"""
TMHMM2.0 Analysis for Chloroplast Proteins
Run transmembrane topology prediction on FASTA sequences
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import re
from concurrent.futures import ProcessPoolExecutor
import tempfile
from tqdm import tqdm

def parse_fasta_headers(fasta_file):
    """Extract sample, gene info from FASTA headers"""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Parse header: >PX136589.1|Gene_atpA|Taxonomy:...
        header_parts = record.id.split('|')
        sample_id = header_parts[0]
        gene_name = header_parts[1].replace('Gene_', '') if len(header_parts) > 1 else 'unknown'
        
        sequences.append({
            'sample': sample_id,
            'gene': gene_name,
            'sequence': str(record.seq),
            'seq_length': len(record.seq)
        })
    return sequences

def run_tmhmm_batch(sequences, output_dir, batch_size=1000):
    """Run TMHMM on batches of sequences"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    results = []
    
    # Process in batches to avoid memory issues
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch = sequences[i:i+batch_size]
        batch_file = output_dir / f"batch_{i//batch_size}.fasta"
        batch_results_file = output_dir / f"batch_{i//batch_size}_tmhmm.txt"
        
        # Write batch FASTA
        with open(batch_file, 'w') as f:
            for j, seq_data in enumerate(batch):
                f.write(f">{seq_data['sample']}|{seq_data['gene']}\n")
                f.write(f"{seq_data['sequence']}\n")
        
        print(f"Processing batch {i//batch_size + 1} ({len(batch)} sequences)...")
        
        # Run TMHMM (using --short format as per HPC docs)
        try:
            cmd = f"tmhmm --short < {batch_file} > {batch_results_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Parse results
            batch_results = parse_tmhmm_output(batch_results_file, batch)
            results.extend(batch_results)
            
            # Clean up batch files
            batch_file.unlink()
            
        except subprocess.CalledProcessError as e:
            print(f"Error running TMHMM on batch {i//batch_size}: {e}")
            continue
    
    return results

def parse_tmhmm_output(tmhmm_file, sequence_batch):
    """Parse TMHMM --short output file"""
    results = []
    seq_lookup = {f"{s['sample']}|{s['gene']}": s for s in sequence_batch}
    
    with open(tmhmm_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            # Short format: seq_id len=123 ExpAA=45.67 First60=0.00 PredHel=2 Topology=o18-40i52-74o
            parts = line.split('\t')
            if len(parts) < 2:
                continue
                
            seq_id = parts[0]
            info = parts[1]
            
            # Parse the info string
            tm_helices = 0
            topology = ""
            
            # Extract PredHel (predicted helices)
            pred_match = re.search(r'PredHel=(\d+)', info)
            if pred_match:
                tm_helices = int(pred_match.group(1))
            
            # Extract Topology
            topo_match = re.search(r'Topology=(\S+)', info)
            if topo_match:
                topology = topo_match.group(1)
            
            if seq_id in seq_lookup:
                seq_data = seq_lookup[seq_id]
                results.append({
                    'sample': seq_data['sample'],
                    'gene': seq_data['gene'],
                    'seq_length': seq_data['seq_length'],
                    'tm_helices': tm_helices,
                    'topology': topology,
                    'is_membrane_protein': tm_helices > 0
                })
    
    return results

def analyze_tm_by_gene(results_df):
    """Analyze TM predictions by gene family"""
    
    gene_summary = results_df.groupby('gene').agg({
        'tm_helices': ['count', 'mean', 'std', 'max'],
        'is_membrane_protein': 'mean',
        'seq_length': 'mean'
    }).round(3)
    
    gene_summary.columns = ['n_sequences', 'mean_tm_helices', 'std_tm_helices', 
                           'max_tm_helices', 'prop_membrane', 'mean_length']
    
    return gene_summary.sort_values('mean_tm_helices', ascending=False)

def main():
    """Main analysis pipeline"""
    
    # Configuration
    fasta_dir = "data/AA_seqs"
    output_dir = "data/tmhmm_results"
    
    print("=== TMHMM Analysis for Chloroplast Proteins ===\n")
    
    # Check if TMHMM is installed
    try:
        result = subprocess.run(["which", "tmhmm"], capture_output=True, text=True)
        if result.returncode != 0:
            raise FileNotFoundError
        print(f"TMHMM found at: {result.stdout.strip()}")
    except FileNotFoundError:
        print("ERROR: TMHMM not found. Please ensure it's in PATH:")
        print("export PATH=/programs/TMHMM2.0c/bin:$PATH")
        return
    
    # Find all FASTA files
    fasta_files = list(Path(fasta_dir).glob("*_AA.fasta"))
    print(f"Found {len(fasta_files)} FASTA files")
    
    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        return
    
    # Collect all sequences
    all_sequences = []
    for fasta_file in fasta_files:
        print(f"Reading {fasta_file.name}...")
        sequences = parse_fasta_headers(fasta_file)
        all_sequences.extend(sequences)
    
    print(f"\nTotal sequences: {len(all_sequences)}")
    
    # Run TMHMM analysis
    print("\nRunning TMHMM predictions...")
    results = run_tmhmm_batch(all_sequences, output_dir)
    
    if not results:
        print("No results generated. Check TMHMM installation.")
        return
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Save raw results
    results_df.to_csv(f"{output_dir}/tmhmm_results.csv", index=False)
    print(f"\nResults saved to {output_dir}/tmhmm_results.csv")
    
    # Gene-level analysis
    print("\n=== Gene-Level Analysis ===")
    gene_summary = analyze_tm_by_gene(results_df)
    print(gene_summary.head(20))
    
    gene_summary.to_csv(f"{output_dir}/gene_tm_summary.csv")
    
    # Basic statistics
    print(f"\n=== Summary Statistics ===")
    print(f"Total sequences analyzed: {len(results_df)}")
    print(f"Membrane proteins (TM > 0): {results_df['is_membrane_protein'].sum()}")
    print(f"Proportion membrane proteins: {results_df['is_membrane_protein'].mean():.3f}")
    print(f"Mean TM helices: {results_df['tm_helices'].mean():.2f}")
    print(f"Max TM helices: {results_df['tm_helices'].max()}")
    
    # Top membrane protein genes
    print(f"\n=== Top Membrane Protein Genes ===")
    membrane_genes = gene_summary[gene_summary['prop_membrane'] > 0.5].head(10)
    print(membrane_genes[['mean_tm_helices', 'prop_membrane']])

if __name__ == "__main__":
    main()
