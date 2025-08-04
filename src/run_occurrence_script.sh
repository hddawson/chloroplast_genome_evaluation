#!/bin/bash
# Wrapper script for the R occurrence script
# Usage: ./src/run_occurrence_script.sh <lui> "<species_name>"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <lui> <species_name>"
    exit 1
fi

lui="$1"
# Join all remaining arguments as the species name (in case it has spaces)
shift
species_name="$*"

# Remove surrounding quotes if they exist
species_name=$(echo "$species_name" | sed 's/^"//;s/"$//')

# Run the R script with the arguments
pixi run Rscript src/pull_species_occurrences.r "$lui" "$species_name"

# Optional: Add some logging
echo "$(date): Processed $lui - $species_name" >> parallel_processing.log 