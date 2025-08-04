library(data.table)
data <- fread("data/combined_occurrences_clean.csv")
data$LUI <- as.factor(data$LUI)
count <- data[, .N, by = LUI]
hist(count$N)

sum(count$N > 5)
