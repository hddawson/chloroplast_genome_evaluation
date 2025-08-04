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

def extract_metadata_from_genomes(id_list, batch_size=100):
    print("Extracting metadata from GenBank records...")
    taxonomy_data = []

    for i in tqdm(range(0, len(id_list), batch_size), desc="Processing batches"):
        batch_ids = id_list[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="nuccore", id=",".join(batch_ids), rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "genbank")

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
                    "Organism": organism,
                    "Taxonomy": "; ".join(taxonomy),
                    "Year": year,
                    "SequencingTech": tech
                })
            handle.close()
            time.sleep(0.5)
        except Exception as e:
            print(f"Error in batch {i//batch_size + 1}: {e}")

    return pd.DataFrame(taxonomy_data)

if __name__ == "__main__":
    Entrez.email = "hdd29@cornell.edu"

    taxonomy_file = "data/taxonomy_info.csv"
    os.makedirs(os.path.dirname(taxonomy_file), exist_ok=True)

    query = 'chloroplast[Title] AND "complete genome"[Title]'
    print("Retrieving total record count...")
    handle = Entrez.esearch(db="nuccore", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    total_count = int(record["Count"])
    print(f"Total matching records: {total_count}")

    ids = fetch_chloroplast_ids(query, total_count, retmax=500)
    if ids:
        print(f"Retrieved {len(ids)} IDs")
        df = extract_metadata_from_genomes(ids, batch_size=300)
        df.to_csv(taxonomy_file, index=False)
        print(f"Metadata saved to {taxonomy_file}")
    else:
        print("No IDs found.")
