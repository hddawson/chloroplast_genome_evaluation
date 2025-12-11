from Bio import Entrez
import pandas as pd
from tqdm import tqdm
import time

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
    
    # Fetch TaxIds from docsums
    print("Fetching taxonomy IDs...")
    taxids = []
    batch_size = 200
    for i in tqdm(range(0, len(all_ids), batch_size)):
        batch_ids = all_ids[i:i+batch_size]
        handle = Entrez.esummary(db="nuccore", id=",".join(batch_ids))
        records = Entrez.read(handle)
        handle.close()
        for rec in records:
            taxid = rec.get('TaxId', None)
            if taxid:
                taxids.append(str(taxid))
        time.sleep(0.5)

    # Get unique taxids
    unique_taxids = list(set(taxids))
    print(f"Unique taxonomy IDs: {len(unique_taxids)}")
    
    # Fetch organism names from taxonomy database
    print("Fetching organism names from taxonomy database...")
    taxid_to_organism = {}
    batch_size = 200
    for i in tqdm(range(0, len(unique_taxids), batch_size)):
        batch_taxids = unique_taxids[i:i+batch_size]
        handle = Entrez.efetch(db="taxonomy", id=",".join(batch_taxids), retmode="xml")
        tax_records = Entrez.read(handle)
        handle.close()
        for tax_rec in tax_records:
            taxid = tax_rec.get('TaxId', 'Unknown')
            organism = tax_rec.get('ScientificName', 'Unknown')
            rank = tax_rec.get('Rank', 'Unknown')
            lineage = tax_rec.get('Lineage', 'Unknown')
            taxid_to_organism[taxid] = {
                'organism': organism,
                'rank': rank,
                'lineage': lineage
            }
        time.sleep(0.5)
    
    organisms = [taxid_to_organism[tid]['Organism'] for tid in unique_taxids if tid in taxid_to_organism]
    
    # Count unique species
    unique_species = set(organisms)
    print(f"Unique species: {len(unique_species)}")
    
    return {
        "total_records": total_count,
        "unique_species": len(unique_species),
        "species_list": sorted(unique_species),
        "taxid_info": taxid_to_organism
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
    out_dir = "results/mitoDiversity"
    
    # Save detailed taxonomic information
    mammal_tax_df = pd.DataFrame([
        {"TaxId": tid, "Organism": info['organism'], "Rank": info['rank'], "Lineage": info['lineage']}
        for tid, info in mammal_results["taxid_info"].items()
    ])
    mammal_tax_df.to_csv(f"{out_dir}/mammal_taxonomy_details.csv", index=False)

    plant_tax_df = pd.DataFrame([
        {"TaxId": tid, "Organism": info['organism'], "Rank": info['rank'], "Lineage": info['lineage']}
        for tid, info in plant_results["taxid_info"].items()
    ])
    plant_tax_df.to_csv(f"{out_dir}/plant_taxonomy_details.csv", index=False)
    
    # Save species lists
    pd.DataFrame({"Species": mammal_results["species_list"]}).to_csv(
        f"{out_dir}/mammal_species_list.csv", index=False)
    pd.DataFrame({"Species": plant_results["species_list"]}).to_csv(
        f"{out_dir}/plant_species_list.csv", index=False)

    print("\nFiles saved: mitochondrial_diversity_summary.csv, mammal/plant_taxonomy_details.csv, mammal/plant_species_list.csv")