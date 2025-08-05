#!/usr/bin/env python3
"""
Script to check that each generated LUI has a corresponding genome file
"""

import pandas as pd
import os
from pathlib import Path

def check_genome_files():
    """Check that each LUI has a corresponding genome file"""
    
    # Read the generated LUIs
    print("Reading generated LUIs...")
    try:
        with open('data/generated_luis.txt', 'r') as f:
            luis = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print("Error: data/generated_luis.txt not found. Run generate_luis.py first.")
        return
    
    print(f"Found {len(luis)} LUIs to check")
    
    # Genome directory path
    genome_dir = Path("../theRefseqening/theRefseqening/data/genomes")
    
    if not genome_dir.exists():
        print(f"Error: Genome directory not found: {genome_dir}")
        return
    
    # Check each LUI
    missing_files = []
    found_files = []
    
    for lui in luis:
        genome_file = genome_dir / f"{lui}.fa"
        if genome_file.exists():
            found_files.append(lui)
        else:
            missing_files.append(lui)
    
    # Report results
    print(f"\n=== Results ===")
    print(f"Total LUIs: {len(luis)}")
    print(f"Found genome files: {len(found_files)}")
    print(f"Missing genome files: {len(missing_files)}")
    print(f"Success rate: {len(found_files)/len(luis)*100:.1f}%")
    
    if missing_files:
        print(f"\nFirst 10 missing genome files:")
        for lui in missing_files[:10]:
            print(f"  - {lui}.fa")
        
        if len(missing_files) > 10:
            print(f"  ... and {len(missing_files) - 10} more")
    
    # Save lists for further use
    if found_files:
        with open('data/luis_with_genomes.txt', 'w') as f:
            for lui in found_files:
                f.write(f"{lui}\n")
        print(f"\nSaved {len(found_files)} LUIs with genomes to: data/luis_with_genomes.txt")
    
    if missing_files:
        with open('data/luis_missing_genomes.txt', 'w') as f:
            for lui in missing_files:
                f.write(f"{lui}\n")
        print(f"Saved {len(missing_files)} missing LUIs to: data/luis_missing_genomes.txt")

if __name__ == "__main__":
    check_genome_files() 