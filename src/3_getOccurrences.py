#!/usr/bin/env python3
"""
Script to process species occurrence queries using multiprocessing
takes in a taxonomy_info.csv file with ID and Organism names

Prepares a list of organism names and their corresponding ID, and feeds this into a bash script for processing
the occurrences 
"""

import subprocess
import time
import sys
from pathlib import Path
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 3_getOccurrences.py <taxonomy_info.csv> <output_dir>")
        sys.exit(1)

    taxonomy_file = sys.argv[1]
    output_dir = sys.argv[2]
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    tax_data = pd.read_csv(taxonomy_file)
    assert 'ID' in tax_data.columns and 'Organism' in tax_data.columns, "CSV must contain 'ID' and 'Organism' columns"

    #the list of queries are the unique organism names
    queries = tax_data[['Organism']].drop_duplicates().values.tolist()

    print(f"Processing {len(queries)} queries...")

    #write the organism queries to a file
    queries_file = Path(output_dir) / "queries.txt"
    with open(queries_file, 'w') as f:
        for organism in queries:
            f.write(f"{organism[0]}\n")

    #call the bash script to process the occurrences
    occurenceCaller = "/workdir/hdd29/chloroplast_genome_evaluation/src/occurrenceCaller.sh"

    command = f"bash {occurenceCaller} {queries_file} {output_dir}"
    print(f"Executing command: {command}")
    start_time = time.time()
    subprocess.run(command, shell=True, check=True)
    end_time = time.time()
    print(f"Occurrences processed in {end_time - start_time:.2f} seconds.")