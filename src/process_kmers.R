#!/usr/bin/env Rscript

#Script to process kmer files and create presence/absence matrix
#Usage: Rscript process_kmers.R [kmer_dir] [output_file] [num_cores]


# Load required libraries
library(data.table)
library(parallel)
library(Matrix)
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default values
kmer_dir <- "jellyfish_output/k21"
output_file <- "data/kmer_presence_absence_matrix.rds"
num_cores <- 19

# Override defaults with command line arguments
if (length(args) >= 1) kmer_dir <- args[1]
if (length(args) >= 2) output_file <- args[2]
if (length(args) >= 3) num_cores <- as.integer(args[3])

cat("=== Kmer Processing Script ===\n")
cat("Input directory:", kmer_dir, "\n")
cat("Output file:", output_file, "\n")
cat("Number of cores:", num_cores, "\n\n")



# Function to process a single kmer file and return presence/absence
process_kmer_file <- function(file_path) {
  tryCatch({
    # Read the kmer file
    kmers <- fread(file_path, select = 1)  # Only read the kmer column (V1)
    
    # Create a data.table with LUI and presence indicator
    lui <- extract_lui_from_filename(basename(file_path))
    
    # Return a data.table with LUI and all kmers present
    result <- data.table(
      LUI = lui,
      kmer = kmers$V1,
      present = 1
    )
    
    return(result)
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "-", e$message))
    return(NULL)
  })
}

# Check if input directory exists
if (!dir.exists(kmer_dir)) {
  stop("Error: Input directory '", kmer_dir, "' does not exist")
}

# Load the list of LUIs to process
lui_file <- "data/n_status_luis.txt"
if (!file.exists(lui_file)) {
  stop("Error: LUI file '", lui_file, "' does not exist")
}

target_luis <- fread(lui_file, header = FALSE, col.names = "LUI")
cat("Loaded", nrow(target_luis), "target LUIs from", lui_file, "\n")

# Get all kmer count files
all_kmer_files <- list.files(kmer_dir, pattern = "_counts\\.tsv$", full.names = TRUE)

if (length(all_kmer_files) == 0) {
  stop("Error: No kmer files found in '", kmer_dir, "'")
}


# Extract LUIs from filenames and filter to only target LUIs
extract_lui_from_filename <- function(filename) {
  base_name <- gsub("_counts\\.tsv$", "", filename)
  return(base_name)
}

file_luis <- sapply(basename(all_kmer_files), extract_lui_from_filename)
kmer_files <- all_kmer_files[file_luis %in% target_luis$LUI]

if (length(kmer_files) == 0) {
  stop("Error: No kmer files found for target LUIs in '", kmer_dir, "'")
}

cat("Found", length(kmer_files), "kmer files matching target LUIs out of", length(all_kmer_files), "total files\n")

# Show which target LUIs are missing
missing_luis <- target_luis$LUI[!target_luis$LUI %in% file_luis]
if (length(missing_luis) > 0) {
  cat("Warning:", length(missing_luis), "target LUIs not found in kmer directory:\n")
  if (length(missing_luis) <= 10) {
    cat(paste(missing_luis, collapse = ", "), "\n")
  } else {
    cat(paste(head(missing_luis, 10), collapse = ", "), "... and", length(missing_luis) - 10, "more\n")
  }
}

# Adjust number of cores if needed
num_cores <- min(num_cores, length(kmer_files))
cat("Processing files using", num_cores, "cores...\n")

# Process files in parallel
cat("Starting parallel processing...\n")
start_time <- Sys.time()

results_list <- mclapply(kmer_files, process_kmer_file, mc.cores = num_cores)

# Remove NULL results (failed files)
results_list <- results_list[!sapply(results_list, is.null)]

if (length(results_list) == 0) {
  stop("Error: No files were successfully processed")
}

cat("Successfully processed", length(results_list), "files\n")

# Combine all results
cat("Combining results...\n")
all_kmers <- rbindlist(results_list, fill = TRUE)

all_kmers$lui_id <- as.integer(factor(all_kmers$LUI))
all_kmers$kmer_id <- as.integer(factor(all_kmers$kmer))

cat("Processed", length(unique(all_kmers$LUI)), "genomes\n")
cat("Total unique kmers:", length(unique(all_kmers$kmer)), "\n")

cat("Creating presence/absence matrix...\n")
# Create sparse matrix
sparse_mat <- sparseMatrix(
  i = all_kmers$lui_id,
  j = all_kmers$kmer_id,
  x = 1L,
  dims = c(length(unique(all_kmers$lui_id)), length(unique(all_kmers$kmer_id))),
  dimnames = list(
    levels(factor(all_kmers$LUI)),
    levels(factor(all_kmers$kmer))
  )
)


# Create output directory if it doesn't exist
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(sparse_mat, output_file)
#saveRDS(presence_absence, output_file)

cat("Created presence/absence matrix with dimensions:", nrow(sparse_mat), "x", ncol(sparse_mat), "\n")


cat("Output saved to:", output_file, "\n") 