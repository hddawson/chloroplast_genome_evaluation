import pandas as pd
import os
from Bio import SeqIO
from tqdm import tqdm
import numpy as np
from sklearn.decomposition import PCA
import csv
import re 
from collections import Counter

def get_gc(dna_sequence):
    """Calculate the GC content of a DNA sequence."""
    if not isinstance(dna_sequence, str):
        raise ValueError("Input must be a string representing a DNA sequence.")
    dna_sequence = dna_sequence.upper()  # Ensure the sequence is in uppercase
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    return gc_count / len(dna_sequence) if len(dna_sequence) > 0 else 0

def get_dinucleotide_frequencies(dna_sequence):
    """Calculate the frequency of each dinucleotide in a DNA sequence."""
    if not isinstance(dna_sequence, str):
        raise ValueError("Input must be a string representing a DNA sequence.")
    dna_sequence = dna_sequence.upper()  # Ensure the sequence is in uppercase
    dinucleotides = [dna_sequence[i:i+2] for i in range(len(dna_sequence) - 1)]
    frequencies = {f"{a}{b}": dinucleotides.count(f"{a}{b}") / (len(dinucleotides) or 1) for a in 'ACGT' for b in 'ACGT'}
    return frequencies

def get_dipeptide_frequencies(aa_sequence):
    """Calculate the frequency of each dipeptide in an amino acid sequence."""
    if not isinstance(aa_sequence, str):
        raise ValueError("Input must be a string representing an amino acid sequence.")
    aa_sequence = aa_sequence.upper()  # Ensure the sequence is in uppercase
    dipeptides = [aa_sequence[i:i+2] for i in range(len(aa_sequence) - 1)]
    frequencies = {f"{a}{b}": dipeptides.count(f"{a}{b}") / (len(dipeptides) or 1) for a in 'ACDEFGHIKLMNPQRSTVWY' for b in 'ACDEFGHIKLMNPQRSTVWY'}
    return frequencies


def get_codon_usage(cds):
    """Calculate the frequency of each codon in a CDS sequence."""
    if not isinstance(cds, str):
        raise ValueError("Input must be a string representing a CDS sequence.")
    check  = len(cds) % 3 == 0 #"CDS length must be a multiple of 3"
    cds = cds.upper()  # Ensure the sequence is in uppercase
    codons = [cds[i:i+3] for i in range(0, len(cds) - 2, 3)]
    frequencies = {codon: codons.count(codon) / (len(codons) or 1) for codon in set(codons)}
    #if check is False: replace all the values with NA
    if not check:
        for codon in frequencies:
            frequencies[codon] = np.nan
    else:
        # Normalize frequencies by dividing by the total count
        total_codons = sum(frequencies.values())
        if total_codons > 0:
            frequencies = {codon: count / total_codons for codon, count in frequencies.items()}

    return frequencies


def get_AA_freq_and_total(proteome):
    """
    Calculate the frequencies of the 20 canonical amino acids and the total count.

    This function only considers the 20 standard amino acids and will include
    a frequency of 0 for any canonical amino acid not found in the input.

    Parameters:
        proteome (str): A string representing a proteome (concatenated amino acid sequences).

    Returns:
        dict: A dictionary containing the frequency of each of the 20 canonical 
              amino acids, plus a 'Total_Amino_Acids' key with the total count.
    """
    # Define the 20 canonical amino acids
    CANONICAL_AA = "ACDEFGHIKLMNPQRSTVWY"

    # Initialize the frequency dictionary with all canonical AAs set to 0.0
    AA_freq = {aa: 0.0 for aa in CANONICAL_AA}

    # Count all characters in the input proteome to get the correct total
    all_counts = Counter(proteome)
    total_AA = sum(all_counts.values())
    
    # Handle the edge case of an empty proteome
    if total_AA == 0:
        AA_freq['Total_Amino_Acids'] = 0
        return AA_freq

    # Calculate frequencies for the canonical amino acids that are present
    for aa, count in all_counts.items():
        if aa in CANONICAL_AA:
            AA_freq[aa] = count / total_AA
    
    # Add the total number of amino acids to the final dictionary
    AA_freq['Total_Amino_Acids'] = total_AA
    
    return AA_freq

def get_proteome_from_gbf(gbf_file):
    """
    Extract the proteome (all translations) from a GenBank file, removing non-amino acid characters.
    
    Parameters:
        gbf_file (str): Path to the GenBank file.

    Returns:
        str: Concatenated and cleaned protein sequences from all CDS features with translations.
    """
    proteome = []
    valid_amino_acids = re.compile(r'[ACDEFGHIKLMNPQRSTVWY]')

    for record in SeqIO.parse(gbf_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                #convert amino acid string to upper just in case 
                AAs = feature.qualifiers["translation"][0]
                AAs = AAs.upper()
                # Clean translation by keeping only valid amino acids
                cleaned_translation = "".join(valid_amino_acids.findall(AAs))
                proteome.append(cleaned_translation)
    
    return "".join(proteome)

def get_gene_group_concatenated_sequences(
    gene_group: str,
    gbf_file: str,
    feature_types: list = ['CDS', 'tRNA', 'rRNA'],
    dna_sequence: bool = True
) -> str:
    """
    Extracts and concatenates sequences from a GenBank file based on a gene prefix
    and specified feature types (e.g., CDS, tRNA, rRNA).

    Args:
        gene_group (str): The gene name prefix to search for (e.g., "psb", "trn").
        gbf_file (str): Path to the GenBank file.
        feature_types (list): A list of feature types to search within 
                              (e.g., ['CDS', 'tRNA', 'rRNA']).
        dna_sequence (bool): If True, returns concatenated DNA sequences. 
                             If False, returns concatenated protein sequences 
                             (only from features with a '/translation' qualifier, like CDS).

    Returns:
        str: A single string of the concatenated sequences.
    """
    if not isinstance(gene_group, str) or len(gene_group) < 2:
        raise ValueError("gene_group must be a string with at least two characters.")
    if not isinstance(feature_types, list):
        raise ValueError("feature_types must be a list or tuple of strings.")

    sequences = []
    target_features = set(feature_types) # Use a set for faster lookups

    for record in SeqIO.parse(gbf_file, "genbank"):
        for feature in record.features:
            # Check if the feature is one of the types we're interested in
            if feature.type not in target_features:
                continue

            # Check if the feature has a 'gene' qualifier that matches the prefix
            gene_name = feature.qualifiers.get("gene", [""])[0]
            if gene_name.startswith(gene_group):
                if dna_sequence:
                    # For any matching feature, extract its DNA sequence
                    sequences.append(str(feature.extract(record.seq)))
                else:
                    # Only extract protein if a translation is available
                    if 'translation' in feature.qualifiers:
                        sequences.append(feature.qualifiers['translation'][0])
    
    return ''.join(sequences)


def get_genome_from_gbf(gbf_file):
    """
    Extract the genome sequence from a GenBank file.
    
    Parameters:
        gbf_file (str): Path to the GenBank file.

    Returns:
        str: Genome sequence.
    """
    for record in SeqIO.parse(gbf_file, "genbank"):
        return str(record.seq)
    