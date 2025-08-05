require(rtry)
library(data.table)
library(stringr)

N_data <- rtry::rtry_import("data/43299.txt")
summary(N_data)

N_data$TraitName <- as.factor(N_data$TraitName)
unique(N_data$TraitName)

fixation_status <- N_data[N_data$TraitName == "Plant nitrogen(N) fixation capacity", ]
fixation_status <- fixation_status[,c("AccSpeciesName","OrigValueStr")]
table(fixation_status$OrigValueStr)

tax_data <- read.csv('data/taxonomy_info.csv')

head(tax_data)
head(fixation_status)

setDT(tax_data)
setDT(fixation_status)

clean_species <- function(name) {
  #specifiers <- c("x", "sp", "var", "subsp", "f", "spp", "cf", "aff", "ex")
  words <- str_split(name, "\\s+")[[1]]
  #cleaned <- words[!tolower(str_remove(words, "\\.")) %in% specifiers]
  paste(head(words, 2), collapse = " ")
}

# Match the Python pattern: alphanumeric ID + alphanumeric species name (no underscore)
tax_data[, Species := sapply(Organism, clean_species)]
tax_data[, LUI := paste0(gsub("[^A-Za-z0-9]", "", ID), gsub("[^A-Za-z0-9]", "", gsub(" ", "", Species)))]

fixation_status[, AccSpeciesName := tolower(AccSpeciesName)]
tax_data[, Organism := tolower(Organism)]
tax_data_unique <- unique(tax_data, by = "Organism")

merged <- fixation_status[tax_data_unique, on = .(AccSpeciesName = Organism)]
length(unique(merged$AccSpeciesName))

sum(is.na(merged$OrigValueStr))

merged <- na.omit(merged)
table(merged$OrigValueStr)

merged[, FixationLabel := fcase(
  tolower(OrigValueStr) %in% c("n-fixer", "n fixer", "n2 fixing", "n2 fixing?", "yes", "y", "true", "1", "yes, an n fixer", "high"), "Yes",
  tolower(OrigValueStr) %in% c("no", "no-n-fixer", "non fixer", "no, not an n fixer", "none", "0", "n"), "No",
  default = "Uncertain"
)]
table(merged$FixationLabel)

merged[, .(n_labels = uniqueN(FixationLabel)), by = AccSpeciesName][n_labels > 1]

conflicts <- merged[AccSpeciesName %in% merged[, .N, by = .(AccSpeciesName, FixationLabel)][, unique(AccSpeciesName[duplicated(AccSpeciesName)])]]

cleaned <- merged[FixationLabel != "Uncertain"]
resolved <- cleaned[, .N, by = .(AccSpeciesName, FixationLabel)][
  order(-N), .SD[1], by = AccSpeciesName
][, .(AccSpeciesName, FixationLabel)]

table(resolved$FixationLabel)

resolved <- merge(resolved, tax_data[, .(Organism, LUI, Species)], 
                  by.x = "AccSpeciesName", by.y = "Organism", all.x = TRUE)

resolved <- resolved[!duplicated(LUI), ]

resolved <- resolved[, .N, by = .(LUI, AccSpeciesName, FixationLabel)][
   order(-N), .SD[1], by = LUI
][, .(LUI, AccSpeciesName, FixationLabel)]

resolved <- merge(resolved, tax_data[, .(LUI, Taxonomy)], by = "LUI", all.x = TRUE)

# Filter to include only Magnoliopsida
resolved <- resolved[grepl("Magnoliopsida", Taxonomy, ignore.case = TRUE), ]

sum(is.na(resolved$FixationLabel))
length(unique(resolved$LUI))
plot(as.factor(resolved$FixationLabel),main="nitrogen fixation status")

output_file <- "data/n_status_luis.txt"

# Write only the LUI column to the file
fwrite(resolved[, .(LUI)], file = output_file, col.names = FALSE)

mers <- fread("jellyfish_output/k31/AB0422403Triticumaestivum_counts.tsv")
head(mers)
hist(mers$V2)

library(parallel)

# Function to extract LUI from filename
extract_lui_from_filename <- function(filename) {
  # Remove the _counts.tsv extension and extract the LUI part
  base_name <- gsub("_counts\\.tsv$", "", filename)
  # The LUI is the part before the species name (alphanumeric ID + alphanumeric species)
  # This matches the pattern from your Python script
  return(base_name)
}


# Function to process a single kmer file and return presence/absence
process_kmer_file <- function(file_path) {
  tryCatch({
    # Read the kmer file
    kmers <- fread(file_path, select = 1)  # Only read the kmer column (V1)
    
    # Create a data.table with LUI and presence indicator
    lui <- extract_lui_from_filename(basename(file_path))
    
    cat(lui)
    
    # Return a data.table with LUI and all kmers present
    result <- data.table(
      LUI = lui,
      kmer = kmers$V1,
      present = 1
    )
    
    return(lui)
  }, error = function(e) {
    warning(paste("Error processing file:", file_path, "-", e$message))
    return(NULL)
  })
}

# Get all kmer count files
kmer_dir <- "jellyfish_output/k31"
kmer_files <- list.files(kmer_dir, pattern = "_counts\\.tsv$", full.names = TRUE)


results_list <- mclapply(kmer_files, process_kmer_file, mc.cores = 20)
results_list <- results_list[!sapply(results_list, is.null)]

all_kmers <- rbindlist(results_list, fill = TRUE)

presence_absence <- dcast(all_kmers, LUI ~ kmer, value.var = "present", fill = 0)

cat("Sparsity (percentage of zeros):", 
    round((1 - sum(presence_absence[, -1]) / (nrow(presence_absence) * (ncol(presence_absence) - 1))) * 100, 2), "%\n")
print(head(presence_absence[, 1:10]))

cat("\nMemory usage of presence/absence matrix:", 
    round(object.size(presence_absence) / 1024^2, 2), "MB\n")

