library(dplyr)
data <- read.csv("data/aaindex1.csv")

pca_result <- prcomp(data[,-1], scale. = TRUE)
summary(pca_result)
str(summary(pca_result))

pcastats <- data.frame(t(summary(pca_result)$importance))
loadings <- data.frame(pca_result$rotation)

aaPCs <- data.frame(pca_result$x)
aaPCs$AminoAcid <- data$X

write.csv(aaPCs, "data/aaPCs.csv")