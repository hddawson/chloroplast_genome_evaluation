import pandas as pd
import os
import csv

#THE LUI is my universal indicator, which is the accession ID and the first two words of the organism name

#this script will take in a taxonomy file, and add the columns LUI and species to it

# species name will be determined by taking the organism name, which may have some common specifiers in it that prevent simply taking the first two words
#such as "x.". "sp.", "var.", "subsp.", "f.", "spp.", "cf.", "aff.", "ex", etc.
#the LUI will be the accession ID and the first two words of the species name,
#and the species name will be the first two words of the species name, with the common specifiers removed

def get_lui_and_species(organism_name, accession_id):
    common_specifiers = ["x.", "sp.", "var.", "subsp.", "f.", "spp.", "cf.", "aff.", "ex"]
    # Remove common specifiers
    for specifier in common_specifiers:
        organism_name = organism_name.replace(specifier, "")
    # Get the first two words
    species_name = " ".join(organism_name.split()[:2])
    lui = f"{accession_id}_{species_name}"
    return lui, species_name


