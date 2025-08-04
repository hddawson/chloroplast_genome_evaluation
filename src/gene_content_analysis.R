gene_content_analysis <- "yippee!"
library(data.table)
library(ggplot2)

df <- fread("/workdir/hdd29/chloroplast_genome_evaluation/data/chloroplast_genes_gc_content.csv")
df$LUI <- str_split_i(gc_data$sample,"_",1)
str(df)
#convert this to a binary
df_binary <- unique(df[, .(LUI, gene)])

# Add a column to mark presence
df_binary[, present := 1L]

# Reshape to wide format: rows = LUI, columns = gene names
df_wide <- dcast(df_binary, LUI ~ gene, value.var = "present", fill = 0)

# View the result
print(df_wide)


# Assuming df_wide from before: rows = LUI, columns = gene names
# Exclude 'LUI' column and count present genes per sample
df_completeness <- df_wide[, .(n_present = rowSums(.SD)), .SDcols = !'LUI']
df_completeness[, LUI := df_wide$LUI]

# Sort by number of present genes
setorder(df_completeness, n_present)
df_completeness[, rank := .I]  # assign a rank for plotting x-axis

# Plot raw gene counts
ggplot(df_completeness, aes(x = rank, y = n_present)) +
  geom_line(color = "darkgreen") +
  labs(
    x = "Samples (sorted by number of genes present)",
    y = "Number of Genes Present",
    title = "Gene Recovery Across Samples"
  ) +
  theme_minimal()
