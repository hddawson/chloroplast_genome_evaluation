# Title..........Geographic data from BIEN and GBIF
# Created at..... 11-06-2023
# Updated at..... 08-11-2025
# Author: henry dawson
if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
  install.packages("CoordinateCleaner", repos = "https://cloud.r-project.org")
}
library(CoordinateCleaner)
# --- Load Required Packages ---
require(BIEN)
require(rgbif)
require(tidyverse)
require(dplyr)
require(data.table) # Optimized for fast file I/O
require(CoordinateCleaner) # For geographic cleaning

# --- Argument Parsing and Setup ---
args <- commandArgs(trailingOnly = TRUE)

# Expecting two arguments:
# args[1]: source_term (e.g., "LUI", a unique ID for your workflow)
# args[2]: query_term (e.g., "Quercus rubra", the species scientific name)
if (length(args) < 2) {
  stop("Usage: Rscript your_script.R <query_term> <output_dir>")
}

query_term <- args[1]  # e.g., "Quercus rubra"
output_dir <- args[2]  # e.g., "/path/to/output"

# Clean species name for file/directory paths
data_name <- gsub(query_term, pattern = ' ', replacement = '_')

# Define output directory based on source_term and cleaned query_term
dir.output <- paste0(output_dir, data_name)

# Create output directory if it doesn't exist
if (!dir.exists(dir.output)) {
  dir.create(dir.output, recursive = TRUE) # recursive=TRUE to create parent dirs if needed
}

# Define subdirectories for metadata
path_species_gbif <- paste0(dir.output, '/species_metadata_gbif')
path_species_bien <- paste0(dir.output, '/species_metadata_bien')

if (!dir.exists(path_species_gbif)) dir.create(path_species_gbif)
if (!dir.exists(path_species_bien)) dir.create(path_species_bien)

# --- 1. Fetch GBIF Data ---
message(paste("Fetching GBIF data for:", query_term))

# wait for 1 second to avoid hitting GBIF API too quickly
Sys.sleep(1)

# Initialize an empty dataframe for GBIF results
gbif_data <- data.frame(
  key = character(), scientificName = character(), decimalLatitude = numeric(),
  decimalLongitude = numeric(), datasetKey = character(), occurrenceID = character(),
  countryCode = character(), kingdom = character(), phylum = character(), class = character(),
  order = character(), family = character(), genus = character(), species = character(),
  taxonRank = character(), individualCount = integer(), publishingOrgKey = character(),
  occurrenceStatus = character(), basisOfRecord = character(), isotherm = numeric(),
  institutionCode = character(), collectionCode = character(), catalogNumber = character(),
  recordNumber = character(), recordedBy = character(), `_score` = numeric(),
  preparations = character(),
  stringsAsFactors = FALSE
) # Define columns to prevent issues if initial search is empty

tryCatch({
  gbif_raw_data <- rgbif::occ_search(search = query_term, hasCoordinate = TRUE, hasGeospatialIssue=FALSE, limit = 10000) 
  if (!is.null(gbif_raw_data$data) && nrow(gbif_raw_data$data) > 0) {
    gbif_data <- as.data.frame(gbif_raw_data$data)
    gbif_data$queryTerm <- query_term  # Add query term for tracking
  }
}, error = function(e) {
  message(paste("Error fetching GBIF data for", query_term, ":", e$message))
})

# Handle GBIF data saving and initial output
if (nrow(gbif_data) == 0) {
  message(paste("No GBIF occurrences found for:", query_term))
  gbif_processed_output <- data.frame(queryTerm = query_term, scientificName=NA, decimalLatitude = NA, decimalLongitude = NA, year = NA)
} else {
  message(paste("Found", nrow(gbif_data), "GBIF occurrences for:", query_term))
  # Save all metadata
  data.table::fwrite(file = paste0(path_species_gbif, '/', data_name, '.csv'), gbif_data)
  gbif_data <- gbif_data[gbif_data$datasetKey != "50c9509d-22c7-4a22-a47d-8c48425ef4a7", ]
  year_col <- if ("year" %in% colnames(gbif_data)) gbif_data$year else NA
  #drop iNaturalist observations, which have datasetKey 50c9509d-22c7-4a22-a47d-8c48425ef4a7
  gbif_processed_output <- data.frame(
    queryTerm = gbif_data$queryTerm,
    scientificName = gbif_data$scientificName,
    decimalLatitude = gbif_data$decimalLatitude,
    decimalLongitude = gbif_data$decimalLongitude,
    year = year_col
  )
}
gbif_processed_output$package <- "GBIF"

# --- 2. Fetch BIEN Data ---
message(paste("Fetching BIEN data for:", query_term))

# Initialize an empty dataframe for BIEN results
bien_data <- data.frame(
  scrubbed_species_binomial = character(), latitude = numeric(), longitude = numeric(),
  dataset = character(), datasource = character(), date_collected = character(),
  locality = character(), country = character(), state_province = character(),
  county = character(), elevation_m = numeric(), identified_by = character(),
  collector = character(), taxon_id = integer(),
  stringsAsFactors = FALSE
) # Define columns for consistency

tryCatch({
  bien_raw_data <- BIEN::BIEN_occurrence_species(species = query_term, cultivated = FALSE, only.geovalid = TRUE)
  if (!is.null(bien_raw_data) && nrow(bien_raw_data) > 0) {
    bien_data <- as.data.frame(bien_raw_data)
  }
}, error = function(e) {
  message(paste("Error fetching BIEN data for", query_term, ":", e$message))
})

# Handle BIEN data saving and initial output
if (nrow(bien_data) == 0) {
  message(paste("No BIEN occurrences found for:", query_term))
  bien_processed_output <- data.frame(queryTerm = query_term, scientificName = NA, decimalLatitude = NA, decimalLongitude = NA, year= NA)
} else {
  message(paste("Found", nrow(bien_data), "BIEN occurrences for:", query_term))
  # Save all metadata
  data.table::fwrite(file = paste0(path_species_bien, '/', data_name, '.csv'), bien_data)
  bien_processed_output <- data.frame(
    queryTerm = query_term,
    scientificName = bien_data$scrubbed_species_binomial, # BIEN uses 'scrubbed_species_binomial'
    decimalLatitude = bien_data$latitude,                  # BIEN uses 'latitude'
    decimalLongitude = bien_data$longitude                 # BIEN uses 'longitude'
  )
  # record date collected into year     
  bien_processed_output$year <- as.numeric(format(as.Date(bien_data$date_collected, format = "%Y-%m-%d"), "%Y"))
}
bien_processed_output$package <- "BIEN"


# --- 3. Merge and Initial Cleaning ---
message("Merging and performing initial data cleaning...")

# Combine GBIF and BIEN data
# Ensure column names match for rbind
names(gbif_processed_output) <- c("queryTerm", "scientificName", "decimalLatitude", "decimalLongitude", "year", "package")
names(bien_processed_output) <- c("queryTerm", "scientificName", "decimalLatitude", "decimalLongitude", "year", "package")

merge_data <- rbind(gbif_processed_output, bien_processed_output)

# Filter out rows with NA coordinates
merge_data <- merge_data %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

# Filter out points at (0,0) or identical sum/diff (often indicative of errors)
merge_data <- merge_data %>%
  mutate(sum_coords = decimalLatitude + decimalLongitude) %>%
  mutate(diff_coords = decimalLatitude - decimalLongitude) %>%
  filter(!(sum_coords == 0 & diff_coords == 0)) %>% # Exclude points at (0,0)
  select(-sum_coords, -diff_coords) # Remove temporary columns

# Remove exact duplicates across all selected columns
merge_data <- merge_data %>%
  distinct(decimalLatitude, decimalLongitude, package, .keep_all = TRUE)

# --- 4. Geographic Cleaning with CoordinateCleaner ---
message("Applying CoordinateCleaner for geographic validation...")

# Initial CoordinateCleaner validation (e.g., impossible coordinates)
# cc_val checks for valid coordinate format/range
inval_flags <- CoordinateCleaner::cc_val(
  merge_data,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged",
  verbose = FALSE
)
merge_data <- merge_data[inval_flags, ]

pirate_flags <- CoordinateCleaner::cc_sea(
  merge_data,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged",
  verbose = FALSE
)
merge_data <- merge_data[pirate_flags, ]

# Remove outliers per species
message("Cleaning spatial outliers (cc_outl)...")

cleaned_outliers_data <- data.frame()
if (nrow(merge_data) > 0) {
  cleaned_outliers_data <- tryCatch(
    {
      CoordinateCleaner::cc_outl(
        x = merge_data,
        lon = "decimalLongitude",
        lat = "decimalLatitude",
        species = "scientificName",
        method = 'mad', # Median Absolute Deviation
        mltpl = 30,     # Multiplier for MAD. Adjust as needed.
        verbose = FALSE
      )
    },
    error = function(e) {
      message(paste("Error during cc_outl for some species:", e$message, "Returning original data for those."))
      return(merge_data) # Return original data if cc_outl errors
    }
  )
}


message(paste("Final number of cleaned records:", nrow(merge_data)))

# --- 5. Output ---
message("Saving final cleaned occurrence data...")

# Ensure correct column order and content for final output
final_output_cols <- c("queryTerm", "scientificName", "decimalLatitude", "decimalLongitude", "year", "package")
merge_data_final <- merge_data %>% select(all_of(final_output_cols))

# Save the data
# Using data_name in the filename ensures it's specific to the queried species
data.table::fwrite(file = paste0(dir.output, '/', data_name, '_occurrences_clean.csv'), merge_data_final)

message(paste("Processing complete for", query_term, ". Output saved to:", dir.output))