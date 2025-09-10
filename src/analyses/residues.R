library(data.table)
library(parallel)
library(arrow)

df        <- as.data.table(read_parquet("data/processed_data.parquet"))
hist(df$Total_Amino_Acids, main="Total AA per plastome")
hist(df$geno_genomeLength)

aln <- as.data.table(read_parquet("data/tmp/aa_supermatrix.parquet"))
aln[1:10,1:10]
unique_counts <- sapply(aln[, -1], function(x) length(unique(x[!is.na(x) & x != ""])))
par(mfrow=c(1,1))
hist(unique_counts)

maj_counts <- sapply(aln[, -1], function(x) {
  clean_x <- x[!is.na(x) & x != ""]  # Remove missing data
  if(length(clean_x) == 0) return(0)
  
  # Find modal (most common) residue
  modal_residue <- names(sort(table(clean_x), decreasing=TRUE))[1]
  
  # Count non-modal residues
  sum(clean_x == modal_residue)
})

# View results
sum(maj_counts==10961)

sum(alternate_counts > 4) / length(aln)
boxplot(alternate_counts)
summary(maj_counts)



ggplot(df, aes(x = Order, y = Total_Amino_Acids)) +
  geom_boxplot() +
  xlab("Order") +
  title("Total Amino Acids per plastome") +
  ylab("Total AA") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

plot(df$geno_genomeLength, df$Total_Amino_Acids)
points(
  df[grep("Fabales", df$Order), ][["geno_genomeLength"]],
  df[grep("Fabales", df$Order), ][["Total_Amino_Acids"]],col="yellow"
)
points(
  df[grep("Poales", df$Order), ][["geno_genomeLength"]],
  df[grep("Poales", df$Order), ][["Total_Amino_Acids"]],col="green"
)
points(
  df[grep("Pinales", df$Order), ][["geno_genomeLength"]],
  df[grep("Pinales", df$Order), ][["Total_Amino_Acids"]],col="blue"
)
points(
  df[grep("Asparagales", df$Order), ][["geno_genomeLength"]],
  df[grep("Asparagales", df$Order), ][["Total_Amino_Acids"]],col="purple"
)
