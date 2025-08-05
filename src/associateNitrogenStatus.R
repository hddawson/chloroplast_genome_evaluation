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


clean_species <- function(name) {
  specifiers <- c("x", "sp", "var", "subsp", "f", "spp", "cf", "aff", "ex")
  words <- str_split(name, "\\s+")[[1]]
  cleaned <- words[!tolower(str_remove(words, "\\.")) %in% specifiers]
  paste(head(cleaned, 2), collapse = " ")
}

# Match the original Python pattern: accession_id + species_name (no spaces)
tax_data[, Species := sapply(Organism, clean_species)]
tax_data[, LUI := paste0(ID, "_", gsub(" ", "", Species))]
tax_data[, LUI := gsub("[^A-Za-z0-9]", "", LUI)]
tax_data[, Species := gsub("[^A-Za-z0-9]", "", Species)]

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
plot(as.factor(resolved$FixationLabel),main="nitrogen fixation status")

output_file <- "data/n_status_luis.txt"

# Write only the LUI column to the file
fwrite(resolved[, .(LUI)], file = output_file, col.names = FALSE)

