from Bio import Entrez
import pandas as pd
from tqdm import tqdm
import time
import os

# Function to fetch all chloroplast genome IDs in batches
def fetch_chloroplast_ids(query, total_records, retmax=500):
    """
    Fetch all chloroplast genome IDs from NCBI.

    Parameters:
        query (str): NCBI search query.
        total_records (int): Total number of records.
        retmax (int): Number of records to fetch per request.

    Returns:
        list: List of NCBI IDs.
    """
    print("Fetching chloroplast genome IDs...")
    all_ids = []

    try:
        for start in tqdm(range(0, total_records, retmax), desc="Fetching IDs"):
            search_handle = Entrez.esearch(db="nuccore", term=query, retstart=start, retmax=retmax)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            all_ids.extend(search_record["IdList"])
            time.sleep(1)  # Rate limit safety
    except Exception as e:
        print(f"Error fetching IDs: {e}")

    return all_ids

# Function to download genomes in batches
def download_genomes(id_list, batch_size=100, output_dir="data/genomes"):
    """
    Download chloroplast genomes using their NCBI IDs.

    Parameters:
        id_list (list): List of NCBI IDs to fetch.
        batch_size (int): Number of genomes to fetch per request.
    """
    print("Downloading chloroplast genomes...")
    for i in tqdm(range(0, len(id_list), batch_size), desc="Downloading batches"):
        batch_ids = id_list[i:i+batch_size]
        try:
            # Fetch the batch of genomes
            fetch_handle = Entrez.efetch(db="nuccore", id=",".join(batch_ids), rettype="fasta", retmode="text")
            sequences = fetch_handle.read()
            fetch_handle.close()

            # Save to a file
            batch_file = os.path.join(output_dir, f"chloroplast_genomes_batch_{i//batch_size + 1}.fasta")
            #print("--------------------------------- writing to batch file ---------------------------------")
            with open(batch_file, "w") as f:
                f.write(sequences)
            
            print(f"Batch {i//batch_size + 1} saved to {batch_file}")
            time.sleep(1)  # Rate limit safety

        except Exception as e:
            print(f"Error downloading batch {i//batch_size + 1}: {e}")

if __name__ == "__main__":
    # Set your email for NCBI's Entrez
    Entrez.email = "hdd29@cornell.edu"

    # Output directory for genomes
    output_dir = "data/downloaded_genomes"
    os.makedirs(output_dir, exist_ok=True)

    # Define the NCBI search query
    search_query = 'chloroplast[Title] AND "complete genome"[Title]'

    # Step 1: Get the total number of matching records
    print("Retrieving total record count...")
    search_handle = Entrez.esearch(db="nuccore", term=search_query, retmax=1)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    total_count = int(search_record["Count"])
    print(f"Total records matching query: {total_count}")

    # Step 2: Fetch all IDs
    chloroplast_ids = fetch_chloroplast_ids(search_query, total_records=total_count, retmax=500)

    # Step 3: Download genomes using the retrieved IDs
    if chloroplast_ids:
        print(f"Total IDs retrieved: {len(chloroplast_ids)}")
        download_genomes(chloroplast_ids, batch_size=100, output_dir=output_dir)
        print(f"All genomes downloaded to '{output_dir}'.")
    else:
        print("No IDs retrieved. Exiting.")
