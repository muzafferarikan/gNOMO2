# Title: Differential Abundance Analysis of Taxa Identified by Metatranscriptomics
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a differential abundance analysis of taxa identified 
#   by metatranscriptomics using Maaslin2.

# Load packages (installed using conda)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(optparse)

commandArgs(trailing=TRUE)

option_list <- list( 
  make_option(c("-g", "--group"), type="character", default="Group", 
              help="enter group name"))

opt <- parse_args(OptionParser(option_list=option_list))


# Get metadata
metadata <- import_qiime_sample_data("resources/metadata.txt")

# Get files 
temp <- list.files(path = "results/intermediate_files/kaiju/kaiju_output/MT", 
                   pattern = "*kaiju_summary.tsv", 
                   full.names = TRUE
                   )

# import files as dataframes
for (i in 1:length(temp)){
  print(temp[i])
  assign(temp[i], as.data.frame(read.table(temp[i],
                                           sep="\t", 
                                           header=TRUE)[,c(5,3)]))
  }

# create a list of dataframes
df_list <- mget(ls(pattern = "*kaiju_summary.tsv"))

# merge all dataframes in the list by first column (taxon_id)
merged_df <- df_list %>% reduce(full_join, by='taxon_name') %>% 
                          column_to_rownames(var = "taxon_name")

# remove unnecessary strings and get sample names
pats = c("results/intermediate_files/kaiju/kaiju_output/MT/MT_|_kaiju_summary.tsv")
temp_samplenames <- str_replace_all(temp, pats, "")

# rename columns with sample names
names(merged_df) <- c(temp_samplenames)

# remove "unclassified" and viral sequences
merged_df <- merged_df[!(rownames(merged_df) == "unclassified"), ]
merged_df <- merged_df[!(rownames(merged_df) == "Viruses"), ]

maaslin_results = Maaslin2::Maaslin2(input_data = merged_df,
                                     input_metadata = metadata,
                                     output = "results/final/diff_abun/taxa-maaslin2-MT",
                                     fixed_effects = opt$group,
                                     standardize = FALSE
                                     )

merged_df2 <- rownames_to_column(merged_df, var="taxa")
merged_df2[is.na(merged_df2)] <- 0
merged_df2$taxa <- paste(merged_df2$taxa, "_MT", sep = "")

write.table(merged_df2,
            file = "results/final/MT/mt_taxa_abundance.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
            )






