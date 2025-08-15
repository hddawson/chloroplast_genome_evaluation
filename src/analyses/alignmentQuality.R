library(stringr)
gapcounts <- read.csv("data/mergedGapcounts.csv")

hist(gapcounts$gaps)

summary(gapcounts$gaps)


genes <- unique(str_split_i(gapcounts$X, "_", 1))

for (gene in genes) {
  sset <- gapcounts[grep(gene, gapcounts$X),]
  hist(sset$gaps, main=gene)
  #n_notIncluded <- min(sset$gaps)
  #n_included <- max(sset$gaps) - n_notIncluded
  #cutoff <- nmax - nrow()
  abline(v= cutoff, col="red")
  cat()
}


# basically, I want to take the total number of samples aligned per gene
# and keep only the sites with >99% of the samples present