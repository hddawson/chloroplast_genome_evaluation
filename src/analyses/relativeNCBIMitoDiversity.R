library(stringr)

mammal_mitos <- read.csv("results/mitoDiversity/mammal_taxonomy_details.csv")
plant_mitos <- read.csv("results/mitoDiversity/plant_taxonomy_details.csv")

length(unique(mammal_mitos$TaxId))
length(unique(plant_mitos$TaxId))

plant_mitos$Genus <- str_split_i(plant_mitos$Organism, " ", 1)
mammal_mitos$Genus <- str_split_i(mammal_mitos$Organism, " ", 1)

length(unique(mammal_mitos$Genus))
length(unique(plant_mitos$Genus))
