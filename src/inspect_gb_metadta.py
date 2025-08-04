from Bio import Entrez, SeqIO
import os
import random

# Setup
Entrez.email = "hdd29@cornell.edu"
output_dir = "data/sample_chloroplast_gbs_full/"
os.makedirs(output_dir, exist_ok=True)

def fetch_random_chloroplast_ids(query, sample_size=5):
    # Get up to 5000 results to sample from
    search_handle = Entrez.esearch(db="nuccore", term=query, retmax=5000)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    all_ids = search_record["IdList"]
    return random.sample(all_ids, min(sample_size, len(all_ids)))

def download_and_dump_full_metadata(id_list):
    fetch_handle = Entrez.efetch(db="nuccore", id=",".join(id_list), rettype="gb", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "genbank"))

    for record in records:
        record_id = record.id
        filename = os.path.join(output_dir, f"{record_id}.gb")
        SeqIO.write(record, filename, "genbank")
        print("="*80)
        print(f"Record ID: {record_id}")
        print(f"Full description: {record.description}")
        print("\nAnnotations:")
        for key, value in record.annotations.items():
            print(f"  {key}: {value}")
        print("\nFeatures:")
        for feature in record.features:
            print(f"  {feature.type}: {feature.qualifiers}")
        print("\nRaw record saved to:", filename)
        print("="*80, "\n")

    fetch_handle.close()

if __name__ == "__main__":
    query = 'chloroplast[Title] AND "complete genome"[Title]'
    ids = fetch_random_chloroplast_ids(query, sample_size=3)
    if ids:
        download_and_dump_full_metadata(ids)
    else:
        print("No sequences found.")
