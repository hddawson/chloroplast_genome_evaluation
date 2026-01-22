#!/usr/bin/env python3
"""
Calculate tRNA stem GC content using tRNAscan-SE structure predictions.

Usage:
    python trna_stem_gc.py --indir data/tmp/alignedGenes/ --outfile stem_gc_results.tsv

Output columns:
    accession, gene, stem_gc, loop_gc, stem_length, loop_length
"""

import argparse
import subprocess
import tempfile
import os
import re
from pathlib import Path
from collections import defaultdict


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def degap(seq):
    """Remove gap characters from sequence."""
    return seq.replace('-', '').replace('.', '')


def run_trnascan(seq, trnascan_path="tRNAscan-SE"):
    """
    Run tRNAscan-SE on a single sequence, return (detected_seq, structure) tuple.
    Returns (None, None) if prediction fails.
    
    tRNAscan-SE output format:
        Seq: TGGGCGTGGCCAAG...
        Str: >>>>>>>..>>>...<<<..>>>>
    
    Structure notation:
        > < = stem (base pairs)
        . = loop/unpaired
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        infile = os.path.join(tmpdir, "input.fa")
        ssfile = os.path.join(tmpdir, "output.ss")
        
        # Write sequence to temp file
        with open(infile, 'w') as f:
            f.write(">query\n")
            f.write(seq + "\n")
        
        # Run tRNAscan-SE with secondary structure output
        cmd = [
            trnascan_path,
            "-O",  # organellar (chloroplast) mode
            "-f", ssfile,  # secondary structure output
            "-q",  # quiet
            infile
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        except subprocess.TimeoutExpired:
            return None, None
        
        # Parse secondary structure file
        if not os.path.exists(ssfile):
            return None, None
            
        with open(ssfile) as f:
            content = f.read()
        
        # Extract sequence and structure lines
        seq_match = re.search(r'Seq:\s*([A-Za-z]+)', content)
        str_match = re.search(r'Str:\s*([><.]+)', content)
        
        if seq_match and str_match:
            return seq_match.group(1), str_match.group(1)
        
        return None, None


def gc_content(bases):
    """Calculate GC content from a list of bases."""
    if not bases:
        return None
    gc = sum(1 for b in bases if b.upper() in 'GC')
    return gc / len(bases)


def parse_trna_regions(seq, structure):
    """
    Parse tRNA into distinct structural regions using stack-based pairing.
    
    tRNAscan-SE notation:
        > = 5' side of stem (pairs with <)
        < = 3' side of stem  
        . = unpaired (loop)
    
    Returns dict with bases for each region.
    """
    assert len(seq) == len(structure), f"Length mismatch: seq={len(seq)}, struct={len(structure)}"
    
    n = len(structure)
    
    # Match > with < using a stack
    stack = []
    pairs = []
    
    for i, char in enumerate(structure):
        if char == '>':
            stack.append(i)
        elif char == '<':
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    
    pairs.sort(key=lambda x: x[0])
    
    # Group consecutive pairs into stems
    stems = []
    current_stem = []
    
    for pair in pairs:
        if not current_stem or pair[0] == current_stem[-1][0] + 1:
            current_stem.append(pair)
        else:
            if current_stem:
                stems.append(current_stem)
            current_stem = [pair]
    if current_stem:
        stems.append(current_stem)
    
    # Find loop regions (consecutive . characters)
    loop_ranges = []
    i = 0
    while i < n:
        if structure[i] == '.':
            start = i
            while i < n and structure[i] == '.':
                i += 1
            loop_ranges.append((start, i))
        else:
            i += 1
    
    # Initialize regions
    regions = {
        'acceptor_stem': [],
        'd_stem': [],
        'd_loop': [],
        'anticodon_stem': [],
        'anticodon_loop': [],
        'variable_loop': [],
        't_stem': [],
        't_loop': [],
    }
    
    # Assign stems (canonical order: acceptor, D, anticodon, T)
    stem_names = ['acceptor_stem', 'd_stem', 'anticodon_stem', 't_stem']
    
    for idx, stem in enumerate(stems[:4]):
        if idx < len(stem_names):
            name = stem_names[idx]
            for (i, j) in stem:
                regions[name].append(seq[i])
                regions[name].append(seq[j])
    
    # Assign loops based on position relative to stems
    if len(stems) >= 4:
        # D-loop: enclosed by D-stem
        d_5_end = max(p[0] for p in stems[1])
        d_3_start = min(p[1] for p in stems[1])
        for start, end in loop_ranges:
            if start > d_5_end and end <= d_3_start:
                regions['d_loop'].extend(seq[start:end])
                break
        
        # Anticodon loop: enclosed by anticodon stem
        ac_5_end = max(p[0] for p in stems[2])
        ac_3_start = min(p[1] for p in stems[2])
        for start, end in loop_ranges:
            if start > ac_5_end and end <= ac_3_start:
                regions['anticodon_loop'].extend(seq[start:end])
                break
        
        # T-loop: enclosed by T stem
        t_5_end = max(p[0] for p in stems[3])
        t_3_start = min(p[1] for p in stems[3])
        for start, end in loop_ranges:
            if start > t_5_end and end <= t_3_start:
                regions['t_loop'].extend(seq[start:end])
                break
        
        # Variable loop: between anticodon 3' and T stem 5'
        ac_3_end = max(p[1] for p in stems[2])
        t_5_start = min(p[0] for p in stems[3])
        for start, end in loop_ranges:
            if start >= ac_3_end and end <= t_5_start:
                regions['variable_loop'].extend(seq[start:end])
    
    return regions


def calc_gc_by_region(seq, structure):
    """
    Calculate GC content for each tRNA structural region.
    
    Returns dict with GC and length for each region, plus totals.
    """
    regions = parse_trna_regions(seq, structure)
    
    result = {}
    all_stem_bases = []
    all_loop_bases = []
    
    for name, bases in regions.items():
        result[f'{name}_gc'] = gc_content(bases)
        result[f'{name}_len'] = len(bases)
        
        if 'stem' in name:
            all_stem_bases.extend(bases)
        elif 'loop' in name:
            all_loop_bases.extend(bases)
    
    result['total_stem_gc'] = gc_content(all_stem_bases)
    result['total_stem_len'] = len(all_stem_bases)
    result['total_loop_gc'] = gc_content(all_loop_bases)
    result['total_loop_len'] = len(all_loop_bases)
    
    return result


def extract_accession(header):
    """Extract accession from FASTA header."""
    # Format: PX136589.1|Gene_trnQ-UUG|Taxonomy:...
    return header.split('|')[0]


def extract_gene(filepath):
    """Extract gene name from filename."""
    # Format: trnA-UGC_CDS_aligned.fasta
    name = os.path.basename(filepath)
    match = re.match(r'(trn[A-Z]-[A-Z]{3})', name)
    return match.group(1) if match else name.replace('_CDS_aligned.fasta', '')


def main():
    parser = argparse.ArgumentParser(description='Calculate tRNA stem GC content')
    parser.add_argument('--indir', required=True, help='Directory with aligned FASTA files')
    parser.add_argument('--outfile', required=True, help='Output TSV file')
    parser.add_argument('--trnascan', default='/programs/tRNAscan-SE-2.0.8/bin/tRNAscan-SE',
                        help='Path to tRNAscan-SE binary')
    parser.add_argument('--pattern', default='trn*_aligned.fasta', 
                        help='Glob pattern for input files')
    parser.add_argument('--resume', action='store_true',
                        help='Resume from existing output file')
    args = parser.parse_args()
    
    indir = Path(args.indir)
    assert indir.is_dir(), f"Input directory not found: {indir}"
    
    # Find all tRNA alignment files
    fasta_files = sorted(indir.glob(args.pattern))
    assert len(fasta_files) > 0, f"No files matching {args.pattern} in {indir}"
    print(f"Found {len(fasta_files)} alignment files")
    
    # Column order
    columns = ['accession', 'gene',
               'acceptor_stem_gc', 'acceptor_stem_len',
               'd_stem_gc', 'd_stem_len', 'd_loop_gc', 'd_loop_len',
               'anticodon_stem_gc', 'anticodon_stem_len', 'anticodon_loop_gc', 'anticodon_loop_len',
               'variable_loop_gc', 'variable_loop_len',
               't_stem_gc', 't_stem_len', 't_loop_gc', 't_loop_len',
               'total_stem_gc', 'total_stem_len', 'total_loop_gc', 'total_loop_len']
    
    # Load existing results if resuming
    completed = set()
    if args.resume and os.path.exists(args.outfile):
        with open(args.outfile) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    completed.add((parts[0], parts[1]))  # (accession, gene)
        print(f"Resuming: {len(completed)} entries already completed")
    
    # Open output file (append if resuming, else write fresh)
    if args.resume and os.path.exists(args.outfile):
        outf = open(args.outfile, 'a')
    else:
        outf = open(args.outfile, 'w')
        outf.write('\t'.join(columns) + '\n')
    
    # Cache: degapped_seq -> (detected_seq, structure)
    structure_cache = {}
    
    failed = defaultdict(int)
    written = 0
    
    try:
        from tqdm import tqdm
        fasta_iter = tqdm(fasta_files)
    except ImportError:
        fasta_iter = fasta_files
    
    for fasta_file in fasta_iter:
        gene = extract_gene(fasta_file)
        sequences = parse_fasta(fasta_file)
        
        for header, aligned_seq in sequences:
            accession = extract_accession(header)
            
            # Skip if already completed
            if (accession, gene) in completed:
                continue
            
            seq = degap(aligned_seq)
            
            # Skip very short sequences (likely truncated)
            if len(seq) < 50:
                failed['too_short'] += 1
                continue
            
            # Get structure (from cache or run tRNAscan-SE)
            if seq not in structure_cache:
                try:
                    detected_seq, structure = run_trnascan(seq, args.trnascan)
                except Exception as e:
                    failed['trnascan_error'] += 1
                    print(f"Warning: tRNAscan error for {accession}: {e}")
                    continue
                structure_cache[seq] = (detected_seq, structure)
            else:
                detected_seq, structure = structure_cache[seq]
            
            if structure is None:
                failed['trnascan_failed'] += 1
                continue
            
            assert detected_seq is not None, "Structure exists but detected_seq is None"
            
            # Calculate GC for each region
            try:
                gc_results = calc_gc_by_region(detected_seq, structure)
            except AssertionError as e:
                failed['length_mismatch'] += 1
                continue
            
            if gc_results['total_stem_gc'] is None:
                failed['no_stem'] += 1
                continue
            
            row = {'accession': accession, 'gene': gene, **gc_results}
            
            # Write immediately (no buffering entire dataset)
            row_vals = []
            for col in columns:
                val = row.get(col, '')
                if val is None:
                    row_vals.append('NA')
                elif isinstance(val, float):
                    row_vals.append(f'{val:.4f}')
                else:
                    row_vals.append(str(val))
            outf.write('\t'.join(row_vals) + '\n')
            written += 1
            
            # Flush periodically
            if written % 100 == 0:
                outf.flush()
    
    outf.close()
    
    print(f"\nDone. Wrote {written} new results to {args.outfile}")
    print(f"Unique sequences cached: {len(structure_cache)}")
    if failed:
        print(f"Failed: {dict(failed)}")


if __name__ == '__main__':
    main()

