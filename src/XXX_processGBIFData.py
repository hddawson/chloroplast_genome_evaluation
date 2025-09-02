import pandas as pd
from tqdm import tqdm
import os
import gc
from typing import Set

#download the occurrence set from GBIF website
#large file, but more reproducible than querying species by species
#since individual queries often fail

def load_taxonomy_organisms(taxonomy_file: str) -> Set[str]:
    """Load the organism names from taxonomy_info.csv into a set for fast lookup."""
    print("Loading taxonomy organisms...")
    taxonomy_df = pd.read_csv(taxonomy_file)
    organisms = set(taxonomy_df['Organism'].dropna().unique())
    print(f"Loaded {len(organisms)} unique organisms from taxonomy file")
    #process all terms in organisms to remove punctuation, and keep only alphanumeric
    organisms = {"".join(ch for ch in org if ch.isalnum()) for org in organisms}
    return organisms

def process_gbif_chunks(gbif_file: str, taxonomy_organisms: Set[str], 
                       chunk_size: int = 1000000, output_file: str = "filtered_gbif.csv"):
    """
    Process GBIF data in chunks, filtering out iNaturalist records and 
    keeping only records that match taxonomy organisms.
    """
    INAT_KEY = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"
    
    # Process file in chunks without counting total lines first
    first_chunk = True
    total_kept = 0
    total_processed = 0
    chunk_num = 0
    
    print(f"Processing in chunks of {chunk_size:,} records...")
    
    try:
        chunk_iterator = pd.read_csv(
            gbif_file, 
            chunksize=chunk_size,
            low_memory=False,
            sep='\t',  # GBIF files are often tab-separated
            quoting=3,  # QUOTE_NONE - handle embedded quotes better
            on_bad_lines='skip',  # Skip malformed lines
            dtype={
                'gbifID': 'str',
                'datasetKey': 'str',
                'scientificName': 'str',
                'decimalLatitude': 'float32',
                'decimalLongitude': 'float32',
                'year': 'Int16',
                'month': 'Int8',
                'day': 'Int8'
            }
        )
    except Exception as e:
        print(f"Error with tab separator, trying comma separator: {e}")
        try:
            chunk_iterator = pd.read_csv(
                gbif_file, 
                chunksize=chunk_size,
                low_memory=False,
                sep=',',
                quoting=3,  # QUOTE_NONE
                on_bad_lines='skip',
                dtype={
                    'gbifID': 'str',
                    'datasetKey': 'str',
                    'scientificName': 'str',
                    'decimalLatitude': 'float32',
                    'decimalLongitude': 'float32',
                    'year': 'Int16',
                    'month': 'Int8',
                    'day': 'Int8'
                }
            )
        except Exception as e2:
            print(f"Error with comma separator too: {e2}")
            print("Trying with python engine and minimal options...")
            chunk_iterator = pd.read_csv(
                gbif_file, 
                chunksize=chunk_size,
                engine='python',  # More robust but slower
                on_bad_lines='skip',
                dtype={
                    'gbifID': 'str',
                    'datasetKey': 'str',
                    'scientificName': 'str'
                }
            )
    
    # Use tqdm without total - it will show rate and count
    with tqdm(desc="Processing chunks", unit="chunks") as pbar:
        for chunk in chunk_iterator:
            chunk_num += 1
            # Filter out iNaturalist records
            chunk_filtered = chunk[chunk['datasetKey'] != INAT_KEY]
            
            # Filter to only records with scientificName in taxonomy organisms
            # Use .isin() for vectorized operation - much faster than loops
            chunk_filtered['species_clean'] = chunk_filtered['species'].str.replace(r'[^0-9a-zA-Z]', '', regex=True)
            chunk_filtered['scientificName_clean'] = chunk_filtered['scientificName'].str.replace(r'[^0-9a-zA-Z]', '', regex=True)

            # now compare against taxonomy_organisms
            mask = (
                chunk_filtered['species_clean'].isin(taxonomy_organisms) |
                chunk_filtered['scientificName_clean'].isin(taxonomy_organisms)
            )
            final_chunk = chunk_filtered[mask]
            
            # Track statistics
            chunk_kept = len(final_chunk)
            total_kept += chunk_kept
            total_processed += len(chunk)
            
            # Write to output file
            if chunk_kept > 0:
                final_chunk.to_csv(
                    output_file, 
                    mode='a' if not first_chunk else 'w',
                    header=first_chunk,
                    index=False
                )
                first_chunk = False
            
            # Update progress bar with current statistics
            pbar.set_postfix({
                'chunk': chunk_num,
                'kept': f"{total_kept:,}",
                'processed': f"{total_processed:,}",
                'rate': f"{total_kept/total_processed*100:.2f}%" if total_processed > 0 else "0%"
            })
            pbar.update(1)
            
            # Force garbage collection to manage memory
            del chunk, chunk_filtered, final_chunk
            gc.collect()
    
    print(f"\nProcessing complete!")
    print(f"Total records processed: {total_processed:,}")
    print(f"Total records kept: {total_kept:,}")
    print(f"Retention rate: {total_kept/total_processed*100:.2f}%")
    print(f"Output saved to: {output_file}")

def get_summary_stats(output_file: str):
    """Get basic statistics about the filtered dataset."""
    if not os.path.exists(output_file):
        print("Output file not found. Run processing first.")
        return
    
    print("\nGenerating summary statistics...")
    
    # Read in smaller chunks for summary to avoid memory issues
    chunk_size = 500000
    total_records = 0
    unique_species = set()
    countries = {}
    years = {}
    
    for chunk in tqdm(pd.read_csv(output_file, chunksize=chunk_size), 
                      desc="Analyzing output"):
        total_records += len(chunk)
        
        # Collect unique species
        unique_species.update(chunk['scientificName'].dropna().unique())
        
        # Count countries
        for country in chunk['countryCode'].dropna():
            countries[country] = countries.get(country, 0) + 1
        
        # Count years
        for year in chunk['year'].dropna():
            years[year] = years.get(year, 0) + 1
    
    print(f"\nSummary Statistics:")
    print(f"Total filtered records: {total_records:,}")
    print(f"Unique species: {len(unique_species):,}")
    print(f"Countries represented: {len(countries)}")
    print(f"Year range: {min(years.keys()) if years else 'N/A'} - {max(years.keys()) if years else 'N/A'}")

def main():
    """Main execution function."""
    # File paths
    gbif_file = "data/0011215-250827131500795.csv"
    taxonomy_file = "data/taxonomy_info.csv"
    output_file = "data/filtered_gbif_data.csv"
    
    # Check if input files exist
    if not os.path.exists(gbif_file):
        print(f"ERROR: GBIF file not found: {gbif_file}")
        return
    
    if not os.path.exists(taxonomy_file):
        print(f"ERROR: Taxonomy file not found: {taxonomy_file}")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    print("Starting GBIF data processing...")
    print(f"Input GBIF file: {gbif_file}")
    print(f"Taxonomy file: {taxonomy_file}")
    print(f"Output file: {output_file}")
    print("-" * 50)
    
    # Step 1: Load taxonomy organisms into memory
    taxonomy_organisms = load_taxonomy_organisms(taxonomy_file)
    
    # Step 2: Process GBIF data in chunks
    process_gbif_chunks(gbif_file, taxonomy_organisms, 
                       chunk_size=1000000, output_file=output_file)
    
    # Step 3: Generate summary statistics
    get_summary_stats(output_file)
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()