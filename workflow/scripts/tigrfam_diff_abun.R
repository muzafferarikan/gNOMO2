# Title: Differential Abundance Analysis of Functional Roles
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a differential abundance analysis of functional roles 
#   (or subroles) assigned by TIGRFAM.

# Load packages (installed using conda env R.yaml)
library(phyloseq)
library(tidyverse)
library(Maaslin2)
library(optparse)

# Define command line arguments
commandArgs(trailing=TRUE)

option_list <- list( 
  make_option(c("-g", "--group"), type="character", default="Group", help="enter group name")
  make_option(c("-c", "--covariates"), type = "character", default = "", help = "Enter covariates"),
  make_option(c("-f", "--transformation"), type = "character", default = "LOG", help = "Enter transformation type"),
  make_option(c("-m", "--normalization"), type = "character", default = "TSS", help = "Enter normalization type")
  )
opt <- parse_args(OptionParser(option_list=option_list))

### Import and format metadata, TIGRFAM, and tigrfam result files ###
metadata <- read.delim("resources/metadata.txt", header=TRUE, sep="\t")
rownames(metadata) <- metadata$SampleID

Table_TIGRFAM <- read.delim("resources/Table_TIGRFAM.tsv", header=FALSE, sep="\t")
colnames(Table_TIGRFAM) <- c("TIGR", "ID", "role", "subrole")
Table_TIGRFAM$role_subrole <- paste(Table_TIGRFAM$role, Table_TIGRFAM$subrole, sep = "_")

tigrfam_files <- list.files(path = "results/intermediate_files/tigrfam", pattern = "*TIGRFAM_v2.tsv", full.names = TRUE)
tigrfam_data <- lapply(tigrfam_files, function(file_path) {
  read.table(file_path, sep = "\t", header = FALSE)
})
df_list <- setNames(tigrfam_data, tigrfam_files)

# Read the Snakemake config file
config <- yaml::yaml.load_file("config/config.yaml")

# Get the analysis options from the config file
omics <- config$omics_combination

# Extract coverage values and merge identical TIGRFAM IDs
df_list_modified <- lapply(df_list, function(dataframe) {
  dataframe %>%
    mutate(V1 = as.numeric(str_remove_all(V1, ".*cov_|\\.g|_.*"))) %>%
    group_by(V2) %>%
    summarise(summed_value = sum(V1)) %>%
    replace_na(list(summed_value = 0))
})

# Custom function to perform differential abundance analysis steps
perform_diff_abundance_analysis <- function(prefix, output_suffix) {
  # Filter the list of dataframes based on the pattern
  list_of_df <- df_list_modified[grepl(prefix, names(df_list_modified))]
  
  # Merge all dataframes in the list by first column (V1)
  merged_df <- list_of_df %>%
    reduce(full_join, by = "V2") %>%
    column_to_rownames(var = "V2")
  
  # Remove unnecessary strings and get sample names
  pats <- paste0("results/intermediate_files/tigrfam/", prefix, "_|_TIGRFAM_v2.tsv")
  temp <- list.files(path = paste0("results/intermediate_files/tigrfam"), pattern = paste0(prefix), full.names = TRUE)
  sample_names <- str_replace_all(temp, pats, "")
  
  # Rename columns with sample names
  colnames(merged_df) <- sample_names
  
  merged_df2 <- rownames_to_column(merged_df, var = "TIGR")
  merged_df2[is.na(merged_df2)] <- 0
  
  sample_names_role <- c("role_subrole", sample_names)
  
  merged_annotated <- inner_join(merged_df2, Table_TIGRFAM, by = "TIGR") %>%
    distinct(TIGR, .keep_all = TRUE) %>%
    select(all_of(sample_names_role))
  
  # Define function to sum values within each role
  sum_duplicates <- function(dataframe) {
    dataframe %>%
      group_by(role_subrole) %>%
      summarise(across(where(is.numeric), sum))
  }
  
  # Apply the sum function to the merged annotated dataframe
  merged_annotated_final <- sum_duplicates(merged_annotated) %>%
    column_to_rownames(var = "role_subrole")
  
# Set fixed effects for differential abundance analysis
fixed_effects <- c(opt$group, opt$covariates)

  # Perform Maaslin2 analysis
  maaslin_results <- Maaslin2::Maaslin2(
    input_data = merged_annotated_final,
    input_metadata = metadata,
    output = paste0("results/final/diff_abun/", output_suffix),
    fixed_effects = fixed_effects,
    transform = opt$transformation,
    normalization = opt$normalization,
    standardize = FALSE
  )
}

# Determine which analysis options to run based on the number of options in the config
if (omics >= 3) {
  perform_diff_abundance_analysis("MG", "tigrfam-maaslin2-MG")
}

if (omics >= 4) {
  perform_diff_abundance_analysis("MT", "tigrfam-maaslin2-MT")
}

