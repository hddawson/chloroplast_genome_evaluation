require(rtry)
library(data.table)
library(stringr)
setwd("/workdir/hdd29/chloroplast_genome_evaluation/")
N_data <- rtry::rtry_import("data/TRY_data/43299_04082025234316/43299.txt")
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

#how many of fixation_status$AccSpeciesName are in tax_data$Organism



fixation_status[, FixationLabel := fcase(
  tolower(OrigValueStr) %in% c("n-fixer", "n fixer", "n2 fixing", "yes", "y", "true", "1", "yes, an n fixer", "high"), "Yes",
  tolower(OrigValueStr) %in% c("no", "no-n-fixer", "non fixer", "no, not an n fixer", "none", "0", "n"), "No",
  default = "Uncertain"
)]

fixation_status[, .(n_labels = uniqueN(FixationLabel)), by = AccSpeciesName][n_labels > 1]

conflicts <- fixation_status[AccSpeciesName %in% fixation_status[, .N, by = .(AccSpeciesName, FixationLabel)][, unique(AccSpeciesName[duplicated(AccSpeciesName)])]]

cleaned <- fixation_status[FixationLabel != "Uncertain"]
resolved <- cleaned[, .N, by = .(AccSpeciesName, FixationLabel)][
  order(-N), .SD[1], by = AccSpeciesName
][, .(AccSpeciesName, FixationLabel)]

table(resolved$FixationLabel)

intx <- intersect(resolved$AccSpeciesName, tax_data$Organism)
length(unique(intx))

write.csv(resolved, "data/fixationStatus.csv")

phData <- rtry::rtry_import("data/TRY_data/42869_15072025181342/42869.txt")

dim(phData)
str(phData)

###cleaning the dataset
phData <- phData[which(phData$DataName=="Plant photosynthetic pathway"),]
table(phData$OrigValueStr)

phData$OrigValueStr[which(phData$OrigValueStr=="c3")] <- "C3"
phData$OrigValueStr[which(phData$OrigValueStr=="C3.")] <- "C3"
phData$OrigValueStr[which(phData$OrigValueStr=="c4")] <- "C4"

clean_phData <- phData[which((phData$OrigValueStr) %in% c("C3", "C4", "CAM")),
                       c("AccSpeciesName","OrigValueStr")] 
table(clean_phData$OrigValueSt)

tax_data <- read.csv('data/taxonomy_info.csv')

head(tax_data)
head(clean_phData)

setDT(tax_data)
setDT(clean_phData)

#how many of fixation_status$AccSpeciesName are in tax_data$Organism

resolved_phData <- clean_phData[, .N, by = .(AccSpeciesName, OrigValueStr)][
  order(-N), .SD[1], by = AccSpeciesName
][, .(AccSpeciesName, OrigValueStr)]

intx <- intersect(resolved_phData$AccSpeciesName, tax_data$Organism)
length(unique(intx))

write.csv(resolved_phData,"data/photosyntheticPathways.csv")
