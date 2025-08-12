from Bio import Entrez, SeqIO
import pandas as pd
from tqdm import tqdm
import time
import os
import datetime
import re

def detect_technology(comment):
    comment = comment.lower()
    if re.search(r'\billumina\b|\bmiseq\b|\bnextseq\b|\bhisep\b|\bnovaseq\b', comment):
        return "Illumina"
    elif re.search(r'\bpacbio\b|\bpacific\s+biosciences\b|\bsequel\b|\brs ii\b', comment):
        return "PacBio"
    elif re.search(r'\bnanopore\b|\boxford\s+nanopore\b|\bminion\b|\bpromethion\b', comment):
        return "Oxford Nanopore"
    elif re.search(r'\bsanger\b|\bcapillary\b', comment):
        return "Sanger"
    elif re.search(r'\b454\b|\broche\b', comment):
        return "454/Roche"
    elif re.search(r'\bion[\s_-]?torrent\b|\bpgm\b', comment):
        return "Ion Torrent"
    elif re.search(r'\bsolid\b', comment):
        return "SOLiD"
    else:
        return "Unknown"

def sanitize_filename(name):
    """
    Sanitize a filename by removing all non-alphanumeric characters.
    """
    return ''.join(c for c in name if c.isalnum())

def detect_technology_from_record(record):
    structured = record.annotations.get("structured_comment", {})
    if isinstance(structured, dict):
        for section in structured.values():
            for key, value in section.items():
                if isinstance(value, str):
                    tech = detect_technology(value)
                    if tech != "Unknown":
                        return tech
    comment = record.annotations.get("comment", "")
    return detect_technology(comment)

def fetch_chloroplast_ids(query, total_records, retmax=500):
    print("Fetching chloroplast genome IDs...")
    all_ids = []
    for start in tqdm(range(0, total_records, retmax), desc="Fetching IDs"):
        handle = Entrez.esearch(db="nuccore", term=query, retstart=start, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        all_ids.extend(record["IdList"])
        time.sleep(0.5)
    return all_ids

def extract_metadata_and_download_genomes(id_list, batch_size=100, output_dir="data/downloaded_genomes"):
    print("Processing batches: extracting metadata and downloading genomes...")
    os.makedirs(output_dir, exist_ok=True)
    taxonomy_data = []
    for i in tqdm(range(0, len(id_list), batch_size), desc="Processing batches"):
        batch_ids = id_list[i:i+batch_size]
        try:
            # Fetch GenBank records (metadata)
            handle = Entrez.efetch(db="nuccore", id=",".join(batch_ids), rettype="gb", retmode="text")
            records = list(SeqIO.parse(handle, "genbank"))
            handle.close()

            # Extract metadata and prepare FASTA sequences
            fasta_seqs = []
            for record in records:
                organism = record.annotations.get("organism", "Unknown")
                taxonomy = record.annotations.get("taxonomy", [])
                sequencing_date = record.annotations.get("date", "Unknown")
                year = None
                try:
                    year = datetime.datetime.strptime(sequencing_date, "%d-%b-%Y").year
                except:
                    pass
                tech = detect_technology_from_record(record)
                taxonomy_data.append({
                    "ID": record.id,
                    "FileBasename": sanitize_filename(record.id),
                    "Organism": organism,
                    "Taxonomy": "; ".join(taxonomy),
                    "Year": year,
                    "SequencingTech": tech
                })
                fasta_seqs.append(record.format("fasta"))

            # Write batch fasta file
            batch_file = os.path.join(output_dir, f"chloroplast_genomes_batch_{i//batch_size + 1}.fasta")
            with open(batch_file, "w") as f:
                f.writelines(fasta_seqs)

            time.sleep(0.5)
        except Exception as e:
            print(f"Error in batch {i//batch_size + 1}: {e}")
    return pd.DataFrame(taxonomy_data)

def split_fasta_files(input_dir, output_dir):
    """
    Split multi-sequence FASTA files into individual files.

    Parameters:
        input_dir (str): Directory containing multi-sequence FASTA files.
        output_dir (str): Directory to save individual genome files.
    """
    os.makedirs(output_dir, exist_ok=True)
    for batch_file in os.listdir(input_dir):
        if batch_file.endswith(".fasta"):
            batch_path = os.path.join(input_dir, batch_file)
            print(f"Processing {batch_path}...")

            # Read the multi-sequence FASTA file
            with open(batch_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    # Extract accession and Latin name from the header
                    header_parts = record.description.split()
                    if len(header_parts) < 3:
                        print(f"Warning: Header has fewer parts: {record.description}")
                        continue
                    
                    accession = sanitize_filename(header_parts[0])  # have to remove non-alphanumeric characters

                    # Save the sequence into a new file
                    output_file = os.path.join(output_dir, f"{accession}.fa")
                    with open(output_file, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")

                    print(f"Saved: {output_file}")

if __name__ == "__main__":
    Entrez.email = "hdd29@cornell.edu"

    query = 'chloroplast[Title] AND "complete genome"[Title]'
    print("Retrieving total record count...")
    handle = Entrez.esearch(db="nuccore", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    total_count = int(record["Count"])
    print(f"Total matching records: {total_count}")

    ids = fetch_chloroplast_ids(query, total_count, retmax=500)  # Limiting to first 500 for testing
    if ids:
        print(f"Retrieved {len(ids)} IDs")
        metadata_df = extract_metadata_and_download_genomes(ids, batch_size=100, output_dir="data/downloaded_genomes")
        metadata_file = "data/taxonomy_info.csv"
        os.makedirs(os.path.dirname(metadata_file), exist_ok=True)
        metadata_df.to_csv(metadata_file, index=False)
        print(f"Metadata saved to {metadata_file}")
        print(f"Genomes saved to data/genomes")
    else:
        print("No IDs found.")

    input_dir = "data/downloaded_genomes/"

    # Output directory for individual genome FASTAs
    output_dir = "data/genomes/"
    os.makedirs(output_dir, exist_ok=True)

    split_fasta_files(input_dir, output_dir)

    print(f"Done! Individual genomes saved to {output_dir}")



    
