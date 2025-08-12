#!/usr/bin/env python3
"""
Script to process species occurrence queries using multiprocessing
Handles quoted species names properly and avoids GNU parallel issues
"""

import subprocess
import multiprocessing as mp
import time
import sys
from pathlib import Path

def process_query(args):
    """Process a single query"""
    lui, species_name = args
    
    try:
        # Run the R script
        result = subprocess.run([
            'pixi', 'run', 'Rscript', 
            'src/pull_species_occurrences.r', 
            lui, species_name
        ], capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            print(f"✓ Success: {lui} - {species_name}")
            return True
        else:
            print(f"✗ Failed: {lui} - {species_name}")
            print(f"  Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"⏰ Timeout: {lui} - {species_name}")
        return False
    except Exception as e:
        print(f"💥 Error: {lui} - {species_name}: {e}")
        return False

def parse_queries_file(filename):
    """Parse the queries file and return list of (lui, species_name) tuples"""
    queries = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Find the first space (separates LUI from species name)
            space_idx = line.find(' ')
            if space_idx == -1:
                continue
                
            lui = line[:space_idx]
            species_name = line[space_idx + 1:]
            
            # Remove surrounding quotes from species name
            if species_name.startswith('"') and species_name.endswith('"'):
                species_name = species_name[1:-1]
                
            queries.append((lui, species_name))
    
    return queries

def main():
    if len(sys.argv) != 2:
        print("Usage: python src/process_queries.py <queries_file>")
        sys.exit(1)
    
    queries_file = sys.argv[1]
    
    if not Path(queries_file).exists():
        print(f"Error: File {queries_file} not found")
        sys.exit(1)
    
    # Parse queries
    print("Parsing queries file...")
    queries = parse_queries_file(queries_file)
    print(f"Found {len(queries)} queries to process")
    
    # Set up multiprocessing
    num_cores = min(mp.cpu_count(), 20)  # Cap at 6 for API safety
    print(f"Using {num_cores} processes")
    
    # Process queries
    start_time = time.time()
    
    with mp.Pool(num_cores) as pool:
        results = list(pool.imap(process_query, queries, chunksize=10))
    
    # Summary
    successful = sum(results)
    failed = len(results) - successful
    elapsed = time.time() - start_time
    
    print(f"\n=== Summary ===")
    print(f"Total queries: {len(queries)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Time elapsed: {elapsed:.1f} seconds")
    
    if failed > 0:
        print(f"\nFailed queries will be logged. Check the output above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main() 