#!/usr/bin/env python3
import pandas as pd

# File paths
occurrences_file = "data/combinedOccurrences.csv"
taxonomy_file = "data/taxonomy_info.csv"
output_file = "data/taxonomy_missing_or_under10.csv"
output_file2 = "data/taxonomy_found_try1.csv"

# Read the data
occ_df = pd.read_csv(occurrences_file)
tax_df = pd.read_csv(taxonomy_file)

# Count occurrences by queryTerm
occ_counts = occ_df['queryTerm'].value_counts()

# Get set of species with < 10 occurrences
rare_species = {species for species, count in occ_counts.items() if count < 10}

# Get set of species in taxonomy
tax_species = set(tax_df['Organism'])

# Get set of species in occurrences
occ_species = set(occ_df['queryTerm'])

# Species in taxonomy but NOT in occurrences
missing_species = tax_species - occ_species

# Combine: missing OR rare
target_species = rare_species | missing_species

# Filter taxonomy_info for these species
filtered_tax = tax_df[tax_df['Organism'].isin(target_species)]

#get the complement of the filtered taxonomy
filtered_tax_complement = tax_df[~tax_df['Organism'].isin(target_species)]

# Save to CSV
filtered_tax.to_csv(output_file, index=False)
filtered_tax_complement.to_csv(output_file2, index=False)

print(f"Found {len(missing_species)} species missing entirely from occurrences.")
print(f"Found {len(rare_species)} species with fewer than 10 occurrences.")
print(f"Total filtered species: {len(filtered_tax)}")

#how many unique species in the filtered taxonomy
unique_species_count = filtered_tax['Organism'].nunique()
print(f"Unique species in filtered taxonomy: {unique_species_count}")
print(f"Saved {len(filtered_tax)} rows to {output_file}")

print(f"Saved {len(filtered_tax_complement)} rows to {output_file2}")
print(f"Unique species in filtered taxonomy complement: {filtered_tax_complement['Organism'].nunique()}")
print(f"Total unique species in taxonomy: {tax_df['Organism'].nunique()}")
