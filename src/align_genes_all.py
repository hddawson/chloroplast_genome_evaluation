import os
import shutil
import subprocess
import warnings
import random
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO, BiopythonParserWarning
from tqdm import tqdm


def align_genes(input_fasta, output_fasta):
    """
    Align sequences in a FASTA file using MAFFT.
    """
    subprocess.run(["/programs/mafft/bin/mafft", "--auto", input_fasta], stdout=open(output_fasta, 'w'))

def align_genes_with_threads(input_fasta, output_fasta, threads=4):
    """
    Align sequences in a FASTA file using MAFFT with multithreading.
    
    :param input_fasta: Path to input FASTA file
    :param output_fasta: Path to output aligned FASTA file
    :param threads: Number of threads to use (default: 4)
    """
    for VROOM in range(threads):
      print("VROOM")
    with open(output_fasta, 'w') as out_f:
        subprocess.run(
            ["/programs/mafft/bin/mafft", "--auto","--treeout","--distout", "--adjustdirection",
            "--thread", str(threads),
            "--threadtb", str(threads),
            "--threadit", str(threads),
            input_fasta], 
            stdout=out_f,
            check=True
        )
        
def get_gene_DNA_from_gbf(gene_name,gbf_file):
    """
    Get the DNA sequence of a gene from the GenBank file.
    """
    ret = ""
    with open(gbf_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "gene" in feature.qualifiers:
                    if gene_name in feature.qualifiers["gene"]:
                        # Extract the DNA sequence
                        dna_sequence = str(feature.location.extract(record.seq))
                        ret = dna_sequence
    if ret:
        return ret
    else:
        print(f"Gene {gene_name} not found in {gbf_file}.")
        return None

def get_gene_DNA_from_gbf_directory(gene_name, gbf_dir):
    """
    Get the DNA sequences of a gene from the GenBank files in the dir, returning a dictionary.
    """
    gene_sequences = {}
    gbf_files = [f for f in os.listdir(gbf_dir) if f.endswith(".gbf.fixed")]
    for gbf_file in gbf_files:
        gbf_path = os.path.join(gbf_dir, gbf_file)
        sequence = get_gene_DNA_from_gbf(gene_name, gbf_path)
        if sequence:
            header = f"{gbf_file.split('.')[0]}_{gene_name}"
            gene_sequences[header] = sequence
    return gene_sequences

def extract_all_genes_AA(lui_list, gbf_dir):
    """
    Extract all protein-coding gene sequences (CDS translations) for given LUIs.
    """
    all_genes = {}
    gbf_files = [f for f in os.listdir(gbf_dir) if f.endswith(".gbf.fixed") and f.split(".")[0] in lui_list]
    print(gbf_files)
    for gbf_file in tqdm(gbf_files):
        gbf_path = os.path.join(gbf_dir, gbf_file)
        with open(gbf_path, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]

                        protein_sequence = feature.qualifiers["translation"][0]
                        LUI = gbf_file.split(".")[0]
                        header = f"{LUI}_{gene_name}"

                        if gene_name not in all_genes:
                            all_genes[gene_name] = {}

                        all_genes[gene_name][header] = protein_sequence

    return all_genes

def extract_all_genes_DNA(lui_list, gbf_dir):
    """
    Extract all protein-coding gene sequences (CDS translations) for given LUIs.
    """
    all_genes = {}
    gbf_files = [f for f in os.listdir(gbf_dir) if f.endswith(".gbf.fixed") and f.split(".")[0] in lui_list]
    print(gbf_files)
    for gbf_file in tqdm(gbf_files):
        gbf_path = os.path.join(gbf_dir, gbf_file)
        with open(gbf_path, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "gene":
                        gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                        dna_sequence = str(feature.location.extract(record.seq))
                        
                        LUI = gbf_file.split(".")[0]
                        header = f"{LUI}_{gene_name}"

                        if gene_name not in all_genes:
                            all_genes[gene_name] = {}

                        all_genes[gene_name][header] = dna_sequence

    #chains into, 
    return all_genes

def write_fasta(gene_dict, output_file):
    """
    Write a dictionary of gene sequences to a FASTA file.
    """
    with open(output_file, "w") as f:
        for header, sequence in gene_dict.items():
            f.write(f">{header}\n{sequence}\n")

# Use ProcessPoolExecutor to parallelize MAFFT alignment
def process_gene(gene_name, gene_dict):
    """
    Write FASTA file, align sequences, and save output.
    """
    if not gene_dict:
        return None

    temp_fasta = os.path.join(temp_dir, f"{gene_name}.fasta")
    aligned_fasta = os.path.join(aligned_dir, f"{gene_name}_aligned.fasta")

    write_fasta(gene_dict, temp_fasta)  # Write to FASTA
    align_genes(temp_fasta, aligned_fasta)  # Align with MAFFT

    return f"Aligned sequences for {gene_name} saved to: {aligned_fasta}"

def fasta_to_gene_df(fasta_file):
    """
    Convert a FASTA file to a DataFrame with gene names and sequences.
    """
    gene_data = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            gene_name = record.id
            sequence = str(record.seq)
            gene_data.append({"gene": gene_name, "sequence": sequence})

    return pd.DataFrame(gene_data)

def remove_unwanted_species(df, metadata):
    """
    remove species from the DF that do not have associated metadata
    """
    pre = len(df)
    df["LUI"] = df["gene"].apply(lambda x: x.split("_")[0])
    df = df[df["LUI"].isin(metadata["LUI"])]

    #print(f"Removed {pre - len(df)} sequences without metadata.")

    return df

def filter_gene_no_metadata(file,outfile):
    df = fasta_to_gene_df(file)
    df["length"] = df["sequence"].apply(len)
    mean = df["length"].mean()
    std = df["length"].std()
    #filter the df
    filtered_df = df[(df["length"] >= mean - 3*std -1) & (df["length"] <= mean + 3*std + 1)]
    print(f"Filtered {len(df) - len(filtered_df)} sequences based on length.")
    #save the filtered df to a new fasta file in the filtered directory
    with open(outfile, "w") as handle:
        for index, row in filtered_df.iterrows():
            handle.write(f">{row['gene']}\n{row['sequence']}\n")

    print(f"Aligned sequences for {file} saved to: {outfile}")
    
def make_random_small_fasta(n, inf, out):
    """
    :param n: number of sequences to select
    :param inf: input fasta file
    :param out: output fasta file
    
    Make a small fasta file with random sequences from the input fasta file.
    """
    with open(inf, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        random_records = random.sample(records, n)  # Select 10 random records

    with open(out, "w") as handle:
        SeqIO.write(random_records, handle, "fasta")

    print(f"Random small fasta file saved to: {out}")

def filter_and_align_gene(file,fastas_dir, filtered_dir, aligned_dir,metadata):
    fasta_file = os.path.join(fastas_dir, file)
    df = fasta_to_gene_df(fasta_file)
    print(df.head())
    df = remove_unwanted_species(df,metadata)
    #add length to the df 
    print(len(df))
    df["length"] = df["sequence"].apply(len)
    mean = df["length"].mean()
    std = df["length"].std()
    #filter the df
    filtered_df = df[(df["length"] >= mean - 3*std -1) & (df["length"] <= mean + 3*std + 1)]
    print(f"Filtered {len(df) - len(filtered_df)} sequences based on length.")
    #save the filtered df to a new fasta file in the filtered directory
    filtered_fasta = os.path.join(filtered_dir, file)
    with open(filtered_fasta, "w") as handle:
        for index, row in filtered_df.iterrows():
            handle.write(f">{row['gene']}\n{row['sequence']}\n")

    #align the filtered fasta file
    aligned_fasta = os.path.join(aligned_dir, file)
    align_genes(filtered_fasta, aligned_fasta)
    print(f"Aligned sequences for {file} saved to: {aligned_fasta}")
    
def align_gene(file,fastas_dir, aligned_dir,metadata):
    fasta_file = os.path.join(fastas_dir, file)
    aligned_fasta = os.path.join(aligned_dir, file)
    align_genes(fasta_file, aligned_fasta)
    print(f"Aligned sequences for {file} saved to: {aligned_fasta}")

if __name__ == "__main__":
    gbf = "/workdir/hdd29/chloro_env_adaptation/data/gbfs/NC0079421Glycinemax.gbf.fixed"
    gene_ex = "psbA"
    genes = get_gene_DNA_from_gbf_directory(gene_ex, "data/gbfs/")
    print(genes)

    write_fasta(genes, "data/psbA.fasta")

    filter_and_align_gene("psbA.fasta", "data/", "data/filtered_fastas", "data/aligned_filtered_genes", pd.read_csv("data/processed_data.csv"))

"""
#just aligning genes 
fastas_dir = "data/gene_fastas"
aligned_dir = "data/aligned_genes"

if not os.path.exists(aligned_dir):
os.makedirs(aligned_dir)

envdata = pd.read_csv("data/processed_data.csv")

num_workers = os.cpu_count() - 5  # Use all but 2 cores

# Run alignment in parallel
with ProcessPoolExecutor(max_workers=num_workers) as executor:
futures = {
executor.submit(align_gene, file, fastas_dir, aligned_dir, metadata=envdata): file
for file in os.listdir(fastas_dir)
}
for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning genes"):
result = future.result()
if result:
print(result)  # Print status messages


"""
"""
#filtering and aligning genes 
fastas_dir = "data/gene_fastas"
filtered_dir = "data/filtered_fastas"
aligned_dir = "data/aligned_filtered_genes"
if not os.path.exists(filtered_dir):
os.makedirs(filtered_dir)
if not os.path.exists(aligned_dir):
os.makedirs(aligned_dir)

envdata = pd.read_csv("data/processed_data.csv")

num_workers = os.cpu_count() - 5  # Use all but 2 cores

# Run alignment in parallel
with ProcessPoolExecutor(max_workers=num_workers) as executor:
futures = {
executor.submit(filter_and_align_gene, file, fastas_dir, filtered_dir, aligned_dir, metadata=envdata): file
for file in os.listdir(fastas_dir)
}
for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning genes"):
result = future.result()
if result:
print(result)  # Print status messages

#for file in tqdm(os.listdir(fastas_dir)[:3]):
#    filter_and_align_gene(file, fastas_dir, filtered_dir, aligned_dir,metadata=envdata)
"""
"""
# Suppress Biopython warnings
warnings.simplefilter("ignore", BiopythonParserWarning)

# Define directories
gbf_dir = "data/gbfs"
temp_dir = "data/temp"
aligned_dir = "data/aligned"

# Ensure the necessary directories exist
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(aligned_dir, exist_ok=True)
# Extract all genes for the clade
with open("data/poaceae_LUIs.txt", "r") as file:
lui_list = [line.strip() for line in file]

print(lui_list)

print("runnin runnin runnin")
gene_data = extract_all_genes(lui_list, gbf_dir)

# Number of parallel workers (adjust for your system)
num_workers = os.cpu_count() - 5  # Use all but 2 cores

# Run alignment in parallel
with ProcessPoolExecutor(max_workers=num_workers) as executor:
futures = {executor.submit(process_gene, gene_name, gene_dict): gene_name for gene_name, gene_dict in gene_data.items()}

for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning genes"):
result = future.result()
if result:
print(result)  # Print status messages

print("All alignments completed.")"""
