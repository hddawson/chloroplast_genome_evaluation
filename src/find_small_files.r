#checking annotation resultrs
library(data.table)

data <- fread("data/jobResults.tsv")

hist(data$diskUsage)
# Find "small" files (e.g., less than 1000 KB = 1 MB)
small_files <- data[diskUsage < 10000]

#write the filenames out to a .txt

writeLines(small_files$directoryName, "data/small_files.txt")
