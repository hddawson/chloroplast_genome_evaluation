#!/bin/bash
set -e

# Check if file list is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <file_list>"
    echo "  file_list: Path to a file containing basenames (one per line) to process"
    echo "  Example: $0 data/smallFiles.txt"
    exit 1
fi

file_list="$1"

# Check if file list exists
if [ ! -f "$file_list" ]; then
    echo "Error: File list '$file_list' not found"
    exit 1
fi

#you will want to do this all on a machine with rapid IO. SSD or better, lots of memory.
# Set directories
input_dir="/workdir/hdd29/theRefseqening/data/genomes"  # Directory containing input FASTA files
final_dir="/workdir/hdd29/theRefseqening/data/results_oldAnno"      # Final directory for results
temp_base_dir="/workdir/hdd29/tmp"                       # Temporary base directory
singularity_image="cpgavas2_0.03.sif"          # Singularity image
bind_dir="/workdir/hdd29/theRefseqening/data/genomes"  # Directory to bind in Singularity

# Create directories if they don't exist
mkdir -p "$final_dir"
mkdir -p "$temp_base_dir"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate singularity

# Set TMPDIR for the session
export TMPDIR="$temp_base_dir"

# Function to process each file
process_file() {
  fasta_file="$1"
  base_name=$(basename "$fasta_file" .fa)
  pid_name=$(echo "$base_name" | tr -d '_.')

  echo "Processing $fasta_file..."

  # Temporary working directory
  temp_dir="${temp_base_dir}/dir_${pid_name}"
  mkdir -p "$temp_dir"

  # Run the CLI tool in Singularity
  singularity exec \
    --env TMPDIR="$temp_base_dir" \
    --bind /workdir/hdd29/tmp:/tmp \
    --bind "$bind_dir:/mnt" \
    "$singularity_image" run-cpgavas2 \
    -in "/mnt/$(basename "$fasta_file")" \
    -db 1 \
    -pid "$pid_name" \
    -out "$temp_dir"

  # Check if the process succeeded
  if [ $? -eq 0 ]; then
    echo "Annotation for $base_name completed successfully."
  else
    echo "Error annotating $base_name. Skipping."
    rm -rf "$temp_dir"
    return 1
  fi

  # Copy results from temporary directory to final directory
  final_output_dir="$final_dir/${pid_name}"
  mkdir -p "$final_output_dir"
  cp -r "$temp_dir"/* "$final_output_dir/"
  echo "Copied results for $base_name to $final_output_dir."

  # Clean up temporary working directory
  rm -rf "$temp_dir"
}

export -f process_file
export singularity_image bind_dir final_dir temp_base_dir TMPDIR # Export variables for subprocesses

# Find all .fa files and process them in parallel with limited jobs
cat "$file_list" | parallel --progress -j 50 process_file {}