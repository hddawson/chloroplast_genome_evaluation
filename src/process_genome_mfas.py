import os
from Bio import SeqIO


def sanitize_filename(name):
    """
    Sanitize a filename by removing all non-alphanumeric characters.
    """
    return ''.join(c for c in name if c.isalnum())

def split_fasta_files(input_dir, output_dir):
    """
    Split multi-sequence FASTA files into individual files.

    Parameters:
        input_dir (str): Directory containing multi-sequence FASTA files.
        output_dir (str): Directory to save individual genome files.
    """
    for batch_file in os.listdir(input_dir):
        if batch_file.endswith(".fasta"):
            batch_path = os.path.join(input_dir, batch_file)
            print(f"Processing {batch_path}...")

            # Read the multi-sequence FASTA file
            with open(batch_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    # Extract accession and Latin name from the header
                    header_parts = record.description.split()
                    if len(header_parts) < 3:
                        print(f"Warning: Header has fewer parts: {record.description}")
                        continue
                    
                    accession = header_parts[0]  # First part is accession
                    latin_name = "_".join(header_parts[1:3])  # Concatenate first two Latin name words
                    sanitized_name = sanitize_filename(accession + latin_name)

                    # Save the sequence into a new file
                    output_file = os.path.join(output_dir, f"{sanitized_name}.fa")
                    with open(output_file, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")

                    print(f"Saved: {output_file}")

if __name__ == "__main__":
    # Ensure the script is run as a standalone program
    print("Starting the FASTA file splitting process...")
    # Input directory containing downloaded FASTA batches
    input_dir = "data/downloaded_genomes/"

    # Output directory for individual genome FASTAs
    output_dir = "data/genomes/"
    os.makedirs(output_dir, exist_ok=True)
    print("Splitting FASTA files into individual genomes...")
    split_fasta_files(input_dir, output_dir)
    print(f"Done! Individual genomes saved to {output_dir}")
