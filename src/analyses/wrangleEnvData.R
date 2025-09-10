library(data.table)

#GBIF occurrences, bulk download, filter, 

gbif_data <- fread("data/combinedOccurrences.csv")
length(unique(gbif_data$queryTerm))

gbif_occs <- readRDS("data/cleaned_data_list.Rds")
gbifs <- rbindlist(gbif_occs)
length(unique(gbifs$species))

#bien data, query by species 
bien_data <- fread("data/combinedOccurrences_BIEN.csv")
table(bien_data$queryTerm==bien_data$scientificName)
length(unique(bien_data$queryTerm))
table(bien_data$package) #yay! 
table(bien_data$year < 2000) #I trust BIEN more so it is okay

gbif_hits <- unique(gbifs$species)
bien_hits <- unique(bien_data$queryTerm)

total_hits <- unique(c(gbif_hits, bien_hits)) #14239

data <- fread("data/selected_genomes.csv")
genome_hits <- unique(data$Organism)
sum(genome_hits %in% total_hits)

bien_data_clean <- bien_data[,c("queryTerm","decimalLatitude","decimalLongitude")]

names(gbifs)[1] <- "queryTerm"

occs <- rbind(bien_data_clean, gbifs)

fwrite(occs, "data/combinedCleanOccurrences.csv")
