library(ape)
library(seqinr)

# Read tree
tree_file <- "results/raxml_input/aa_treeConstrained_fiveParsimony_ASR.raxml.ancestralTree"
tree <- read.tree(tree_file)
tip_ids <- tree$tip.label

cat("Tree has", length(tip_ids), "tips\n\n")

# Get all fasta files
fasta_dir <- "data/AA_seqs/"
fasta_files <- list.files(fasta_dir, pattern = "_AA\\.fasta$", full.names = TRUE)
stopifnot(length(fasta_files) > 0)

cat("Found", length(fasta_files), "fasta files\n\n")

# Function to extract ID from header
extract_id <- function(header) {
  sub("\\|.*", "", header)
}

# Open output file for writing
output_file <- "esmc_input.fasta"
if (file.exists(output_file)) file.remove(output_file)

total_written <- 0

# Process each fasta file individually
for (idx in seq_along(fasta_files)) {
  fasta_file <- fasta_files[idx]
  
  cat("[", idx, "/", length(fasta_files), "] Processing:", basename(fasta_file), "\n")
  
  # Read sequences
  seqs <- read.fasta(fasta_file, seqtype = "AA", as.string = TRUE)
  n_total <- length(seqs)
  
  # Filter to tree tips
  filtered_seqs <- list()
  filtered_names <- c()
  
  for (i in seq_along(seqs)) {
    seq_id <- extract_id(names(seqs)[i])
    
    if (seq_id %in% tip_ids) {
      filtered_seqs <- c(filtered_seqs, seqs[i])
      filtered_names <- c(filtered_names, seq_id)
    }
  }
  
  n_kept <- length(filtered_seqs)
  cat("  Kept", n_kept, "/", n_total, "sequences\n")
  
  # Append to output file manually
  if (n_kept > 0) {
    names(filtered_seqs) <- filtered_names
    
    # Open connection for appending
    con <- file(output_file, open = if (idx == 1) "w" else "a")
    for (j in seq_along(filtered_seqs)) {
      writeLines(paste0(">", filtered_names[j]), con)
      writeLines(filtered_seqs[[j]], con)
    }
    close(con)
    
    total_written <- total_written + n_kept
  }
  
  # Clean up
  rm(seqs, filtered_seqs, filtered_names)
  gc(verbose = FALSE)
  
  cat("\n")
}

cat("========================================\n")
cat("Total sequences written:", total_written, "\n")
cat("Output file:", output_file, "\n")
cat("========================================\n")

stopifnot(total_written > 0)