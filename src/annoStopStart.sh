#!/bin/bash
set -e

# Set directories
input_dir="/workdir/hdd29/theRefseqening/data/genomes"  # Directory containing input FASTA files
final_dir="/workdir/hdd29/theRefseqening/data/results"      # Final directory for results
temp_base_dir="/workdir/hdd29/tmp"                       # Temporary base directory
singularity_image="cpgavas2_0.03.sif"          # Singularity image
bind_dir="/workdir/hdd29/theRefseqening/data/genomes"  # Directory to bind in Singularity

# Create directories if they don't exist
mkdir -p "$final_dir"
mkdir -p "$temp_base_dir"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate singularity

export TMPDIR="$temp_base_dir"

# Function to process each file
process_file() {
  fasta_file="$1"
  base_name=$(basename "$fasta_file" .fa)
  pid_name=$(echo "$base_name" | tr -d '_.')
  
  final_output_dir="$final_dir/${pid_name}"
  
  # ⛔ Skip if output already exists and is non-empty 
  if [[ -d "$final_output_dir" && -n "$(ls -A "$final_output_dir" 2>/dev/null)" ]]; then
    echo "Skipping $base_name: output already exists."
    return 0
  fi

  echo "Processing $fasta_file..."
  
  # Create the final directory from the start
  mkdir -p "$final_output_dir"
  
  # Run Singularity, binding the final directory and writing directly to it
  singularity exec \
    --containall \
    --writable-tmpfs \
    --bind "$bind_dir:/mnt" \
    --bind "$final_output_dir:/output" \
    "$singularity_image" run-cpgavas2 \
    -in "/mnt/$(basename "$fasta_file")" \
    -db 1 \
    -pid "$pid_name" \
    -out "/output" \
     > "$final_output_dir/run.log" 2>&1
  
  # Check for success
  if [ $? -eq 0 ]; then
    echo "Annotation for $base_name completed in $final_output_dir."
  else
    echo "Error annotating $base_name. Log is in $final_output_dir."
    # To start a failed run from scratch, you can uncomment the next line
    # rm -rf "$final_output_dir" 
    return 1
  fi
}

export -f process_file
export singularity_image bind_dir final_dir temp_base_dir TMPDIR # Export variables for subprocesses

# Find all .fa files and process them in parallel with limited jobs
find "$input_dir" -name "*.fa" | parallel --progress -j 20 process_file {}
