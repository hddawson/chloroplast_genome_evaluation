#!/usr/bin/env python3
"""
Script to generate LUIs from taxonomy_info.csv
Makes ID column alphanumeric, takes first two words of Organism column,
makes them alphanumeric, and joins them together.
"""

import pandas as pd
import re
import sys

def alphanum_only(s):
    """Remove all non-alphanumeric characters from a string"""
    return re.sub(r'[^A-Za-z0-9]', '', str(s))

def clean_species_name(name):
    """
    Extract first two words from organism name, excluding taxonomic specifiers
    """
    specifiers = {"x", "sp", "var", "subsp", "f", "spp", "cf", "aff", "ex"}
    
    # Split by whitespace and clean
    words = name.split()
    cleaned_words = [w for w in words if w.lower().strip(".") not in specifiers]
    
    # Take first two words
    if len(cleaned_words) >= 2:
        return " ".join(cleaned_words[:2])
    else:
        return " ".join(cleaned_words)

def generate_lui(id_val, organism_name):
    """
    Generate LUI by combining alphanumeric ID with alphanumeric species name
    """
    # Make ID alphanumeric
    clean_id = alphanum_only(id_val)
    
    # Get first two words of organism name and make alphanumeric
    species = clean_species_name(organism_name)
    clean_species = alphanum_only(species)
    
    # Join them together
    lui = f"{clean_id}{clean_species}"
    
    return lui

def main():
    # Read the taxonomy info file
    print("Reading taxonomy_info.csv...")
    df = pd.read_csv('data/taxonomy_info.csv')
    
    print(f"Processing {len(df)} rows...")
    
    # Generate LUIs
    df['LUI'] = df.apply(lambda row: generate_lui(row['ID'], row['Organism']), axis=1)
    
    # Show some examples
    print("\nFirst 10 LUIs generated:")
    for i, row in df.head(10).iterrows():
        print(f"ID: {row['ID']} | Organism: {row['Organism']} | LUI: {row['LUI']}")
    
    # Save to file
    output_file = 'data/generated_luis.txt'
    df['LUI'].to_csv(output_file, index=False, header=False)
    
    print(f"\nGenerated {len(df)} LUIs")
    print(f"Saved to: {output_file}")
    
    # Show unique count
    unique_count = df['LUI'].nunique()
    print(f"Unique LUIs: {unique_count}")
    
    if len(df) != unique_count:
        print(f"Warning: {len(df) - unique_count} duplicate LUIs found")
        
        # Show some duplicates
        duplicates = df[df['LUI'].duplicated(keep=False)].sort_values('LUI')
        if len(duplicates) > 0:
            print("\nFirst few duplicate LUIs:")
            for lui in duplicates['LUI'].unique()[:5]:
                dupes = duplicates[duplicates['LUI'] == lui]
                print(f"LUI: {lui}")
                for _, row in dupes.iterrows():
                    print(f"  - {row['ID']} | {row['Organism']}")

if __name__ == "__main__":
    main() 