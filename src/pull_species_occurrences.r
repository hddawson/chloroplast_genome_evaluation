

# Title..........Geographic data from BIEN and GBIF
# Created at..... 11-06-2023
# Updated at..... 06-24-2024 (and further cleaned by Gemini)
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
  stop("Usage: Rscript your_script.R <source_term> <query_term>")
}

source_term <- args[1] # e.g., "LUI"
query_term <- args[2]  # e.g., "Quercus rubra"

# Clean species name for file/directory paths
data_name <- gsub(query_term, pattern = ' ', replacement = '_')

# Define output directory based on source_term and cleaned query_term
dir.output <- paste0(getwd(), '/data/env_data/', source_term, "_", data_name)

# Check if cleaned output file already exists
output_file_path <- file.path(dir.output, paste0(data_name, '_occurrences_clean.csv'))

if (file.exists(output_file_path)) {
  message(paste("Cleaned output already exists for", query_term, "- skipping."))
  quit(save = "no", status = 0)  # Exit script successfully
}

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
  gbif_raw_data <- rgbif::occ_search(scientificName = query_term, hasCoordinate = TRUE, limit = 100000) # Increased limit
  if (!is.null(gbif_raw_data$data) && nrow(gbif_raw_data$data) > 0) {
    gbif_data <- as.data.frame(gbif_raw_data$data)
  }
}, error = function(e) {
  message(paste("Error fetching GBIF data for", query_term, ":", e$message))
})

# Handle GBIF data saving and initial output
if (nrow(gbif_data) == 0) {
  message(paste("No GBIF occurrences found for:", query_term))
  data.table::fwrite(file = paste0(path_species_gbif, '/missing_species.csv'), data.frame(species = query_term), append = TRUE)
  gbif_processed_output <- data.frame(scientificName = query_term, decimalLatitude = NA, decimalLongitude = NA, year = NA)
} else {
  message(paste("Found", nrow(gbif_data), "GBIF occurrences for:", query_term))
  # Save all metadata
  data.table::fwrite(file = paste0(path_species_gbif, '/', data_name, '.csv'), gbif_data)
  gbif_processed_output <- data.frame(
    scientificName = gbif_data$scientificName,
    decimalLatitude = gbif_data$decimalLatitude,
    decimalLongitude = gbif_data$decimalLongitude,
    year = gbif_data$year
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
  data.table::fwrite(file = paste0(path_species_bien, '/missing_species.csv'), data.frame(species = query_term), append = TRUE)
  bien_processed_output <- data.frame(scientificName = query_term, decimalLatitude = NA, decimalLongitude = NA, year= NA)
} else {
  message(paste("Found", nrow(bien_data), "BIEN occurrences for:", query_term))
  # Save all metadata
  data.table::fwrite(file = paste0(path_species_bien, '/', data_name, '.csv'), bien_data)
  bien_processed_output <- data.frame(
    scientificName = bien_data$scrubbed_species_binomial, # BIEN uses 'scrubbed_species_binomial'
    decimalLatitude = bien_data$latitude,                  # BIEN uses 'latitude'
    decimalLongitude = bien_data$longitude                 # BIEN uses 'longitude'
  )
  # Add year column if it exists, otherwise use NA
  if ("year" %in% colnames(bien_data)) {
    bien_processed_output$year <- bien_data$year
  } else {
    bien_processed_output$year <- NA
  }
}
bien_processed_output$package <- "BIEN"


# --- 3. Merge and Initial Cleaning ---
message("Merging and performing initial data cleaning...")

# Combine GBIF and BIEN data
# Ensure column names match for rbind
names(gbif_processed_output) <- c("scientificName", "decimalLatitude", "decimalLongitude", "year", "package")
names(bien_processed_output) <- c("scientificName", "decimalLatitude", "decimalLongitude", "year", "package")

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
inval_data <- merge_data[!inval_flags, ]
if (nrow(inval_data) > 0) {
  message(paste("Detected", nrow(inval_data), "invalid coordinates based on cc_val. Saving to invalid_coords.csv"))
  data.table::fwrite(file = file.path(dir.output, "invalid_coords.csv"), inval_data)
  # You might want to stop or filter these out completely depending on strictness
  merge_data <- merge_data[inval_flags, ] # Keep only valid ones
} else {
  message("No invalid coordinates detected by cc_val.")
}

# Remove terrestrial records that fall in the ocean
message("Filtering out marine occurrences (cc_sea)...")
pirate_flags <- CoordinateCleaner::cc_sea(
  merge_data,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged",
  verbose = FALSE
)
pirate_data <- merge_data[!pirate_flags, ]
if (nrow(pirate_data) > 0) {
  message(paste("Detected", nrow(pirate_data), "marine occurrences. Saving to pirate_coords.csv"))
  data.table::fwrite(file = file.path(dir.output, "pirate_coords.csv"), pirate_data)
  merge_data <- merge_data[pirate_flags, ] # Keep only terrestrial ones
} else {
  message("No marine occurrences detected by cc_sea.")
}

# Remove outliers per species
message("Cleaning spatial outliers (cc_outl)...")
species_with_enough_points <- merge_data %>%
  group_by(scientificName) %>%
  summarise(n_points = n()) %>%
  pull(scientificName)

# Separate data into those with enough points and those without
data_to_clean <- merge_data %>% filter(scientificName %in% species_with_enough_points)
data_to_keep_as_is <- merge_data %>% filter(!(scientificName %in% species_with_enough_points))

cleaned_outliers_data <- data.frame()
if (nrow(data_to_clean) > 0) {
  cleaned_outliers_data <- tryCatch(
    {
      CoordinateCleaner::cc_outl(
        x = data_to_clean,
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
      return(data_to_clean) # Return original data if cc_outl errors
    }
  )
}

# Combine cleaned data with data that didn't go through outlier cleaning
merge_data_clean <- rbind(cleaned_outliers_data, data_to_keep_as_is)

message(paste("Final number of cleaned records:", nrow(merge_data_clean)))

# Ensure scientificName is correctly formatted (Genus species)
# This step is good practice but might be redundant if your query_term is always "Genus species"
merge_data_clean$scientificName <- sapply(strsplit(as.character(merge_data_clean$scientificName), " "), function(x) {
  paste(x[1], x[2])
})


# --- 5. Output ---
message("Saving final cleaned occurrence data...")

# Ensure correct column order and content for final output
final_output_cols <- c("scientificName", "decimalLatitude", "decimalLongitude", "year", "package")
merge_data_final <- merge_data_clean %>% select(all_of(final_output_cols))

merge_data_final$LUI <- source_term # Add the source term for tracking

# Save the data
# Using data_name in the filename ensures it's specific to the queried species
data.table::fwrite(file = paste0(dir.output, '/', data_name, '_occurrences_clean.csv'), merge_data_final)

message(paste("Processing complete for", query_term, ". Output saved to:", dir.output))