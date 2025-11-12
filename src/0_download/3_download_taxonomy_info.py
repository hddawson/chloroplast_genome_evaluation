from Bio import Entrez
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import time
import os

# Set your email for NCBI's Entrez
Entrez.email = "hdd29@cornell.edu"

# Output directories
output_dir = "data/chloroplast_gbs/"
taxonomy_file = "data/taxonomy_info.csv"
os.makedirs(output_dir, exist_ok=True)

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

# Function to download genomes in batches and extract taxonomic info
def download_genomes_and_extract_taxonomy(id_list, batch_size=100):
    """
    Download chloroplast genomes as GenBank files and extract taxonomy.

    Parameters:
        id_list (list): List of NCBI IDs to fetch.
        batch_size (int): Number of genomes to fetch per request.
    """
    print("Downloading chloroplast genomes and extracting taxonomy...")
    taxonomy_data = []

    for i in tqdm(range(0, len(id_list), batch_size), desc="Downloading batches"):
        batch_ids = id_list[i:i+batch_size]
        try:
            # Fetch the batch of genomes
            fetch_handle = Entrez.efetch(db="nuccore", id=",".join(batch_ids), rettype="gb", retmode="text")
            records = SeqIO.parse(fetch_handle, "genbank")
            
            for record in records:
                # Save GenBank file
                gb_file = os.path.join(output_dir, f"{record.id}.gb")
                with open(gb_file, "w") as f:
                    SeqIO.write(record, f, "genbank")
                
                # Extract taxonomy information
                organism = record.annotations.get("organism", "Unknown")
                taxonomy = record.annotations.get("taxonomy", [])
                taxonomy_data.append({"ID": record.id, "Organism": organism, "Taxonomy": "; ".join(taxonomy)})

            fetch_handle.close()
            time.sleep(1)  # Rate limit safety

        except Exception as e:
            print(f"Error downloading batch {i//batch_size + 1}: {e}")

    # Save taxonomy information to a CSV file
    taxonomy_df = pd.DataFrame(taxonomy_data)
    taxonomy_df.to_csv(taxonomy_file, index=False)
    print(f"Taxonomy information saved to '{taxonomy_file}'.")

# Define the NCBI search query
search_query = '(chloroplast[All Fields] AND "complete genome"[All Fields] AND ("Arabidopsis thaliana"[Organism] OR ("Arabidopsis thaliana"[Organism] OR Arabidopsis thaliana[All Fields]))) AND "Arabidopsis thaliana"[Primary Organism] AND chloroplast[filter]'

# Step 1: Get the total number of matching records
print("Retrieving total record count...")
search_handle = Entrez.esearch(db="nuccore", term=search_query, retmax=1)
search_record = Entrez.read(search_handle)
search_handle.close()
total_count = int(search_record["Count"])
print(f"Total records matching query: {total_count}")

# Step 2: Fetch all IDs
chloroplast_ids = fetch_chloroplast_ids(search_query, total_records=total_count, retmax=500)

# Step 3: Download genomes and extract taxonomy information
if chloroplast_ids:
    print(f"Total IDs retrieved: {len(chloroplast_ids)}")
    download_genomes_and_extract_taxonomy(chloroplast_ids, batch_size=100)
    print(f"All genomes downloaded to '{output_dir}' and taxonomy info saved to '{taxonomy_file}'.")
else:
    print("No IDs retrieved. Exiting.")
