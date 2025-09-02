library(dplyr)
library(CoordinateCleaner)
library(data.table)

setDTthreads(8)
cols <- c("species","decimalLongitude", "decimalLatitude")
df <- fread("data/filtered_gbif_data.csv",
            select = cols)

df <- df %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

df <- df %>%
  mutate(sum_coords = decimalLatitude + decimalLongitude,
         diff_coords = decimalLatitude - decimalLongitude) %>%
  filter(!(sum_coords == 0 & diff_coords == 0)) %>%  # Exclude points at (0,0)
  select(-sum_coords, -diff_coords)

df_unique <- unique(df, by = c("decimalLatitude", "decimalLongitude", "species"))

val_flags <- cc_val(
  df_unique,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged",
  verbose = FALSE
)
df_unique <- df_unique[val_flags, ]

landlubber_flags <- cc_sea(
  df_unique,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "flagged",
  verbose = FALSE
)
df_unique <- df_unique[landlubber_flags, ]
cat(dim(df_unique))
saveRDS(df_unique, "data/bulkCleanOccs.rds")

gbif_data <- readRDS("data/bulkCleanOccs.rds")

# Get unique species for processing
unique_species <- unique(gbif_data$species)

#tax_data <- fread("data/taxonomy_info.csv")
#length(unique(tax_data$Organism))
message(paste("Processing outliers for", length(unique_species), "species"))

# Process species in batches to manage memory
cleaned_data_list <- list()

for (i in 1:length(unique_species)) {
  cat(i)
  spec <- unique_species[i]

  # Get data for this batch of species
  batch_data <- gbif_data[species %in% spec]
  
  # Apply outlier detection to the batch
  cleaned_batch <- tryCatch(
    {
      cc_outl(
        x = batch_data,
        lon = "decimalLongitude",
        lat = "decimalLatitude",
        species = "species",
        method = 'mad',  # Median Absolute Deviation
        mltpl = 10,      # Multiplier for MAD
        verbose = FALSE
      )
    },
    error = function(e) {
      message(paste("Error during cc_outl for batch:", e$message, 
                    "Returning original data for this batch"))
      return(batch_data)
    }
  )
  
  cleaned_data_list[[length(cleaned_data_list) + 1]] <- cleaned_batch
}
saveRDS(cleaned_data_list, "data/cleaned_data_list.Rds")
cleaned_data_list <- readRDS("data/cleaned_data_list.Rds")
# Combine all cleaned batches
message("Combining cleaned data...")
final_cleaned_data <- do.call(rbind, cleaned_data_list)

message(paste("Final number of cleaned records:", nrow(final_cleaned_data)))
message(paste("Removed", nrow(gbif_data) - nrow(final_cleaned_data), 
              "outlier records"))

# Save the cleaned data
message("Saving cleaned data...")
fwrite(final_cleaned_data, output_file)

# Summary statistics
message("\n=== CLEANING SUMMARY ===")
message(paste("Input records:", initial_count))
message(paste("Final records:", nrow(final_cleaned_data)))
message(paste("Total removed:", initial_count - nrow(final_cleaned_data)))
message(paste("Retention rate:", 
              round(100 * nrow(final_cleaned_data) / initial_count, 2), "%"))
message(paste("Unique species in final dataset:", 
              length(unique(final_cleaned_data$scientificName))))
message(paste("Output saved to:", output_file))

