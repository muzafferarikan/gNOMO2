# Title: Preparing Peptide Abundance Tables
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script generates a peptide abundance table using msgf results and 
#   prepares a peptide list for unipept search.

# Load packages (installed using conda)
library(tidyverse)
library(phyloseq)

# Get peptide files
temp <- list.files(path = "results/intermediate_files/msgf_plus", pattern = "peptides_*", full.names = TRUE)

# Import files as vectors
for (i in 1:length(temp)){
  print(temp[i])
  assign(temp[i], readLines(temp[i]))
}

# Determine unique peptides and their counts and convert to df for each sample
for(i in 1:length(temp)){
  assign(temp[i],as.data.frame(table(get(temp[i]))))
}

# Create a list of dataframes
df_list <- mget(ls(pattern = "results/intermediate_files/msgf_plus/peptides_*"))

# Merge all dataframes in the list by first column (peptide seq)
merged_df <- df_list %>% reduce(full_join, by='Var1')

# Remove unnecessary strings and get sample names
pats = c("results/intermediate_files/msgf_plus/peptides_|.txt")
temp_samplenames <- str_replace_all(temp, pats, "")

# Rename columns with sample names
colnames(merged_df) <- append("peptide",temp_samplenames)

# Replace NAs to 0
merged_df[is.na(merged_df)] <- 0

# Export total peptide list for unipept search
write(as.vector(merged_df$peptide), "results/final/MP/unipept_list.txt")

# Save final dataframe for differential abundance analysis
write.table(merged_df, "results/final/MP/pep_abundance_table.txt", quote = FALSE, 
            col.names=TRUE,row.names = FALSE, sep="\t")






