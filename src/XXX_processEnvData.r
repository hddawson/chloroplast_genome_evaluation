library(data.table)
library(stringr)

env_dir <- "data/geoDataOut"
#recursively find all csv files in the directory and subdirectories
csv_files <- list.files(env_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)[1:10]

dt <- fread(csv_files[1])

#define a function to process each file,
# read the file, group by queryTerm, and calculate the 10th 50th and 90th percentiles of the value column
# the basename of the file is the value column
# also return as n_occurrences, the number of the subgroups

process_file <- function(file) {
  dt <- fread(file)
  dt <- na.omit(dt)
  
  # Extract the column name from the filename
  colname <- tools::file_path_sans_ext(basename(file))
  colname <- substr(colname, 0, nchar(colname)-8)
  if (!colname %in% names(dt)) {
    stop(paste("Column", colname, "not found in", file))
  }
  
  # Compute percentiles
  result <- dt[, .(
    p10 = quantile(get(colname), 0.1, na.rm = TRUE),
    p50 = quantile(get(colname), 0.5, na.rm = TRUE),
    p90 = quantile(get(colname), 0.9, na.rm = TRUE)
    ), by = queryTerm]
  
  # Rename percentile columns dynamically
  setnames(result, old = c("p10", "p50", "p90"), 
           new = paste0(colname, c("_p10", "_p50", "_p90")))
  
  return(result)
}

# Process files
results <- lapply(csv_files, process_file)
result <- Reduce(function(x, y) merge(x, y, by = "queryTerm", all = TRUE, sort = FALSE), results)

#get n_occurrences from a representative file 
#not T_Opt_site - NA at noisy sites

dt <- fread(csv_files[2])
dt <- na.omit(dt)

n_occs <- dt[, .(
  n_occurrences = .N
), by = queryTerm]

result <- na.omit(result)

data <- merge(result,n_occs,by = "queryTerm")

sum(data$n_occurrences < 4)
par(mfrow=c(1,1))
hist(data$Topt_site_p50)
hist(data$Topt_site_p10)
hist(data$Topt_site_p90)
hist(log10(data$n_occurrences))
hist(log10(data$n_occurrences))

plot(data$wc2.1_2.5m_bio_10_p90,data$Topt_site_p90*0.01)
abline(a=0,b=1,col="red")


filtered_data <- data[which(data$n_occurrences > 4),]

#check thge merge
#g_data <- fread("data/selected_genomes.csv")

#g_hits <- unique(g_data$Organism)
#o_hits <- unique(filtered_data$queryTerm)

#hits <- intersect(g_hits, o_hits)

fwrite(filtered_data, "data/pheno.csv")

cat(dim(filtered_data))

