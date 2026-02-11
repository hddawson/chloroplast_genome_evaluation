#!/usr/bin/env python3
"""
Generate 4x4 grid of hormone-receptor combinations for AF3 GPU step.

Takes existing _data.json files (with MSAs already computed) and swaps
the ligand SMILES to create all 16 combinations.
"""

import json
import shutil
from pathlib import Path
import argparse

# The 4 hormones with their SMILES
HORMONES = { #maybe the backslash escaping is messing w these??
    "auxin_IAA": "C1=CC=C2C(=C1)C(=CN2)CC(=O)O",
    "ethylene": "C=C",
    "gibberellin_GA3": "CC12C(C=CC3(C1C(C45C3CCC(C4)(C(=C)C5)O)C(=O)O)OC2=O)O",
    "abscisic_acid": "CC1=CC(=O)CC([C@]1(/C=C/C(=C\C(=O)O)/C)O)(C)C",  # simplified ABA
    "jasmonic_acid": "CC/C=C\C[C@@H]1[C@H](CCC1=O)CC(=O)O",
}

# Map from source directories to protein names
# Adjust these paths based on your actual directory names
SOURCE_DIRS = {
    "PYR1": "abscisic_acid_PYR1_20260130_140851",
    "ETR1": "ethylene_ETR1_20260130_141735",
    "GID1A": "gibberellin_GA3_GID1A_20260130_144405",
    "COI1": "jasmonic_acid_COI1_20260130_145303",
}

def load_data_json(path: Path) -> dict:
    """Load the _data.json file."""
    with open(path) as f:
        return json.load(f)

def modify_ligand(data: dict, new_smiles: str) -> dict:
    """Replace the ligand SMILES in the data."""
    # Find and modify the ligand entry in sequences
    for seq in data.get("sequences", []):
        if "ligand" in seq:
            seq["ligand"]["smiles"] = new_smiles
            break
    return data

def generate_grid(output_base: Path, source_base: Path):
    """Generate all 16 combinations."""
    
    output_base.mkdir(parents=True, exist_ok=True)
    
    generated = []
    
    for protein_name, source_dir in SOURCE_DIRS.items():
        source_path = source_base / source_dir
        
        # Find the _data.json file
        data_files = list(source_path.glob("*_data.json"))
        if not data_files:
            print(f"Warning: No _data.json found in {source_path}")
            continue
        
        data_file = data_files[0]
        original_data = load_data_json(data_file)
        
        for hormone_name, smiles in HORMONES.items():
            # Create new job name
            job_name = f"{hormone_name}_{protein_name}"
            
            # Create output directory
            job_dir = output_base / job_name
            job_dir.mkdir(exist_ok=True)
            
            # Modify data with new ligand
            new_data = json.loads(json.dumps(original_data))  # deep copy
            new_data["name"] = job_name
            modify_ligand(new_data, smiles)
            
            # Write new _data.json
            output_file = job_dir / f"{job_name}_data.json"
            with open(output_file, "w") as f:
                json.dump(new_data, f, indent=2)
            
            generated.append(job_name)
            print(f"Generated: {output_file}")
    
    print(f"\nTotal: {len(generated)} combinations generated")
    return generated

def main():
    parser = argparse.ArgumentParser(
        description="Generate 4x4 hormone-receptor grid for AF3"
    )
    parser.add_argument("--source", "-s", required=True,
                       help="Source directory with existing AF3 outputs")
    parser.add_argument("--output", "-o", required=True,
                       help="Output directory for new combinations")
    
    args = parser.parse_args()
    
    source_base = Path(args.source)
    output_base = Path(args.output)
    
    assert source_base.exists(), f"Source not found: {source_base}"
    
    generate_grid(output_base, source_base)

if __name__ == "__main__":
    main()