from Bio import Entrez
import pandas as pd
from tqdm import tqdm
import time
from Bio import SeqIO

def fetch_species_count(query, description):
    """Fetch and count unique species for a given query."""
    print(f"\n{description}")
    
    # Get total count
    handle = Entrez.esearch(db="nuccore", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    total_count = int(record["Count"])
    print(f"Total records: {total_count}")
    
    # Fetch all IDs
    print("Fetching IDs...")
    all_ids = []
    retmax = 500
    for start in tqdm(range(0, total_count, retmax)):
        handle = Entrez.esearch(db="nuccore", term=query, retstart=start, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        all_ids.extend(record["IdList"])
        time.sleep(0.5)
    
    # Fetch organisms in batches
    print("Fetching organism names...")
    organisms = []
    batch_size = 1000
    for i in tqdm(range(0, len(all_ids), batch_size)):
        batch_ids = all_ids[i:i+batch_size]
        handle = Entrez.efetch(db="nuccore", id=",".join(batch_ids), 
                               rettype="gb", retmode="text")
        # Parse only annotations, not full sequence
        records = SeqIO.parse(handle, "genbank")
        for rec in records:
            organism = rec.annotations.get("organism", "Unknown")
            organisms.append(organism)
        handle.close()
        #time.sleep(0.5)
    
    # Count unique species
    unique_species = set(organisms)
    print(f"Unique species: {len(unique_species)}")
    
    return {
        "total_records": total_count,
        "unique_species": len(unique_species),
        "species_list": sorted(unique_species)
    }

if __name__ == "__main__":
    Entrez.email = "hdd29@cornell.edu"
    
    # Mammalian mitochondria
    mammal_query = '(mitochondrion[filter] AND "complete genome"[All Fields] AND Mammalia[Organism])'
    mammal_results = fetch_species_count(mammal_query, "MAMMALIAN MITOCHONDRIA")
    
    # Plant mitochondria
    plant_query = '(mitochondrion[filter] AND "complete genome"[All Fields] AND Viridiplantae[Organism])'
    plant_results = fetch_species_count(plant_query, "PLANT MITOCHONDRIA")
    
    # Save summary
    summary = pd.DataFrame([
        {"Group": "Mammals", "Total_Records": mammal_results["total_records"], 
         "Unique_Species": mammal_results["unique_species"]},
        {"Group": "Plants", "Total_Records": plant_results["total_records"], 
         "Unique_Species": plant_results["unique_species"]}
    ])
    
    summary.to_csv("mitochondrial_diversity_summary.csv", index=False)
    print("\n=== SUMMARY ===")
    print(summary)
    
    # Save species lists
    pd.DataFrame({"Species": mammal_results["species_list"]}).to_csv(
        "mammal_species_list.csv", index=False)
    pd.DataFrame({"Species": plant_results["species_list"]}).to_csv(
        "plant_species_list.csv", index=False)
    
    print("\nFiles saved: mitochondrial_diversity_summary.csv, mammal_species_list.csv, plant_species_list.csv")