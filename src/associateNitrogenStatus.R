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

intx <- intersect(fixation_status$AccSpeciesName, tax_data$Organism)
length(unique(intx))

fixation_status[, FixationLabel := fcase(
  tolower(OrigValueStr) %in% c("n-fixer", "n fixer", "n2 fixing", "yes", "y", "true", "1", "yes, an n fixer", "high"), "Yes",
  tolower(OrigValueStr) %in% c("no", "no-n-fixer", "non fixer", "no, not an n fixer", "none", "0", "n"), "No",
  default = "Uncertain"
)]

table(fixation_status$FixationLabel)
write.csv(fixation_status, "data/fixationStatus.csv")





mean(fixation_status$AccSpeciesName %in% tax_data$Organism)

clean_species <- function(name) {
  #specifiers <- c("x", "sp", "var", "subsp", "f", "spp", "cf", "aff", "ex")
  words <- str_split(name, "\\s+")[[1]]
  #cleaned <- words[!tolower(str_remove(words, "\\.")) %in% specifiers]
  paste(head(words, 2), collapse = " ")
}

# Match the pattern: alphanumeric ID + alphanumeric species name (no underscore)
tax_data[, Species := sapply(Organism, clean_species)]
tax_data[, LUI := paste0(gsub("[^A-Za-z0-9]", "", ID), gsub("[^A-Za-z0-9]", "", gsub(" ", "", Species)))]

#fixation_status[, AccSpeciesName := tolower(AccSpeciesName)]
#tax_data[, Organism := tolower(Organism)]
tax_data_unique <- unique(tax_data, by = "Organism")

merged <- fixation_status[tax_data_unique, on = .(AccSpeciesName = Organism)]
length(unique(merged$AccSpeciesName))

sum(is.na(merged$OrigValueStr))

merged <- na.omit(merged)
table(merged$OrigValueStr)

merged[, FixationLabel := fcase(
  tolower(OrigValueStr) %in% c("n-fixer", "n fixer", "n2 fixing", "yes", "y", "true", "1", "yes, an n fixer", "high"), "Yes",
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
#save the data
write.csv(
  resolved,"data/nitrogen_status.csv"
)
resolved <- read.csv("data/nitrogen_status.csv")

# Filter to include only Streptophyta
setDT(resolved)
resolved <- resolved[grepl("Streptophyta", Taxonomy, ignore.case = TRUE), ]

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

prabs <- readRDS("data/kmer_presence_absence_matrix.rds")

dim(prabs)
rownames(prabs[1:10])

# Extract LUIs from matrix
matrix_luis <- prabs@Dimnames[[1]]

# Filter resolved to matching LUIs
resolved_sub <- resolved[LUI %in% matrix_luis]

# Reorder matrix rows to match resolved_sub
match_idx <- match(resolved_sub$LUI, matrix_luis)
prabs_aligned <- prabs[match_idx, ]

# Confirm alignment
stopifnot(all(rownames(prabs_aligned) == resolved_sub$LUI))

library(parallel)
y <- as.integer(resolved_sub$FixationLabel == "Yes")

kmer_cols <- lapply(seq_len(ncol(prabs_aligned)), function(j) prabs_aligned[, j])

# New version of the function
run_fisher_vec <- function(col) {
  present <- y[col != 0]
  absent <- y[col == 0]
  
  if (length(present) == 0 || length(absent) == 0) return(NA_real_)
  
  tbl <- matrix(c(sum(present), length(present) - sum(present),
                  sum(absent), length(absent) - sum(absent)), nrow=2)
  
  fisher.test(tbl)$p.value
}

# Run in parallel
cl <- makeCluster(detectCores())
clusterExport(cl, c("kmer_cols", "y", "run_fisher_vec"))
pvals <- parSapply(cl, seq_len(ncol(prabs_aligned)), run_fisher)
stopCluster(cl)



data <- read.csv("data/nitro_merged_data.csv")

data$FixationLabel <- as.factor(data$FixationLabel)
data$FixationLabel <- relevel(data$FixationLabel, ref = "No")


mod <- glm(FixationLabel ~ GC_Content + Largest_IR_Length, data = data, family = binomial)
summary(mod)

fabs <- data[grepl("Fabales", data$Taxonomy_x),]

plot(data$Largest_IR_Length, data$GC_Content)
library(ggplot2)
ggplot(data, aes(x = Genome_Length, y = GC_Content, color = FixationLabel)) +
  geom_point(alpha = 0.6) +
  labs(title = "GC Content vs Genome Length by Nitrogen Fixation Status",
       x = "Genome Length",
       y = "GC Content (%)",
       color = "Fixation Status") +
  theme_minimal()

ggplot(data, aes(x = Largest_IR_Length, y = GC_Content, color = FixationLabel)) +
  geom_point(alpha = 0.6) +
  labs(title = "GC Content vs Largest_IR_Length by Nitrogen Fixation Status",
       x = "Largest_IR_Length",
       y = "GC Content (%)",
       color = "Fixation Status") +
  theme_minimal()

ggplot(data, aes(x = Genome_Length, y = GC_Content, color = FixationLabel)) +
  geom_point(alpha = 0.6) +
  geom_point(data = subset(data, is.na(Largest_IR_Length)),
             aes(x = Genome_Length, y = GC_Content),
             color = "black", shape = 1, size = 2, stroke = 1) +
  labs(title = "GC Content vs Genome Length by Nitrogen Fixation Status",
       x = "Genome Length",
       y = "GC Content (%)",
       color = "Fixation Status") +
  theme_minimal()



mod2 <- lm(GC_Content ~ Largest_IR_Length, data=data)
summary(mod2)
