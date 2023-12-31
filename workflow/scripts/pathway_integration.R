# Title: Pathway Level Integrated Analysis for Multi-Omics Data
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a pathway level integrated analysis of multi-omics
#   data. It covers data loading (emapper results), preprocessing, differential
#   abundance analysis for each omics level and pathway based visualization.

# Load necessary packages
library(phyloseq)
library(tidyverse)
library(optparse)
library(pathview)
library(KEGGREST)
library(gtools)

# Define command line arguments
commandArgs(trailing=TRUE)

option_list <- list( 
  make_option(c("-g", "--group"), type="character", default="Group", 
              help="enter group name"))

opt <- parse_args(OptionParser(option_list=option_list))

################################ FUNCTIONS ####################################

# Function for filling empty cells in result files with "unclassified"
fill_empty_cells <- function(df) {
  df[is.na(df) | df == ""] <- "unclassified"
  return(df)
}

# Function for removing the redundant parts in names in coverage stat results
remove_first_redundant_part <- function(df) {
  df[, 1] <- sub("\\s*#.*$", "", df[, 1])
  df[, 1] <- gsub("cov_[0-9]+\\.[0-9]+\\.NODE_[0-9]+_length_[0-9]+_", "", df[, 1], perl = TRUE)
  return(df)
}

# Function for processing emapper results and output tables for pathview
process_seq <- function(prefix) {
  
  ## Read emapper result files
  eggnog_files <- list.files(
    path = "results/intermediate_files/eggnog",
    pattern = ".*_eggnog.emapper.annotations",
    full.names = TRUE
  )
  
  ## Import emapper result files
  eggnog_data <- lapply(eggnog_files, function(file_path) {
    data <- read.table(
      file_path,
      sep = "\t",
      header = FALSE,
      skip = 4,
      quote = ""
    )
    data[, c(1, 7)]
  })
  
  ## Extract the sample names from the full file names
  sample_name <- gsub("_eggnog.emapper.annotations", "", basename(eggnog_files))
  
  ## Create a dataframe list of emapper results
  emapper_df_list <- setNames(eggnog_data, sample_name)
  
  ## Create a new list of data frames without empty cells
  emapper_df_list_filled <- lapply(emapper_df_list, fill_empty_cells)
  
  ## Read coverage stat files
  aug_prod_file <- file.path("results/intermediate_files/aug_prod", sample_name, "stats.txt")
  
  ## Import coverage result files
  coverage_data <- lapply(aug_prod_file, function(file_path) {
    data <- read.table(
      file_path,
      sep = "\t",
      header = FALSE,
      quote = "",
      skip = 1,
      comment.char = "@"
    )
    data[, 1:2]
  })
  
  ## Create a dataframe list of coverage results
  coverage_df_list <- setNames(coverage_data, sample_name)
  
  ## Apply the redundancy removal function to each data frame in the list
  coverage_df_list_modified <- lapply(coverage_df_list, remove_first_redundant_part)
  
  ## Initialize an empty list to store the merged data frames
  merged_list <- list()
  
  ## Loop through the data frame names
  for (name in names(emapper_df_list_filled)) {
    ## Merge data frames with the same name based on the first column (ID)
    merged_df <- merge(emapper_df_list_filled[[name]], coverage_df_list_modified[[name]], by.x = 1, by.y = 1, all = FALSE)
    
    ## Add the merged data frame to the new list
    merged_list[[name]] <- merged_df
  }
  
  ## Extract coverage values and merge identical KEGG IDs
  emapper_df_list_clean <- lapply(merged_list, function(dataframe) {
    dataframe %>%
      group_by(V7) %>%
      summarise(summed_value = sum(V2)) %>%
      replace_na(list(summed_value = 0))
  })
  
  ## Subset specified omics dataframes based on prefix 
  emapper_df_final_list <- emapper_df_list_clean[grepl(prefix, names(emapper_df_list_clean))] #subset specified omics dataframes
  
  ## Merge dataframes and format
  merged_df <- emapper_df_final_list %>%
    reduce(full_join, by = "V7") %>%
    column_to_rownames(var = "V7")
  
  ## Prepare column names for the dataframe
  pats <- paste0("results/intermediate_files/eggnog/", prefix, "_|_eggnog.emapper.annotations")
  temp <- list.files(path = "results/intermediate_files/eggnog", pattern = paste0(prefix, ".*eggnog\\.emapper\\.annotations"), full.names = TRUE)
  sample_names <- str_replace_all(temp, pats, "")
  colnames(merged_df) <- sample_names
  
  ## Format and prepare for the analysis
  merged_df2 <- rownames_to_column(merged_df, var = "KEGG")
  merged_df2[is.na(merged_df2)] <- 0
  merged_df_uniq <- merged_df2[!grepl(",", merged_df2$KEGG), ]
  rownames(merged_df_uniq) <- NULL
  
  kegg_df <- column_to_rownames(merged_df_uniq, var = "KEGG")
  
  ## Perform differential abundance analysis (eggnog based)
  maaslin_results <- Maaslin2::Maaslin2(
    input_data = kegg_df,
    input_metadata = metadata,
    output = paste0("results/final/diff_abun/eggnog-maaslin2-", prefix),
    fixed_effects = opt$group,
    standardize = FALSE
  )
  
  ## Process results for omics comparisons
  diff_abun <- maaslin_results$results %>% data.frame()
  diff_abun_kegg <- diff_abun$feature
  kegg_df2 <- kegg_df[diff_abun_kegg, ] 
  kegg_df3 <- kegg_df2 %>% mutate_all(~ (. / sum(.)) * 100)
  result <- data.frame(ID = rownames(kegg_df3))
  rownames(result) <- result$ID
  
  unique_groups <- unique(metadata[[opt$group]])
  
  ## Calculate fold changes and log2 ratios
  for (group in unique_groups) {
    group_samples <- metadata$SampleID[metadata[[opt$group]] == group]
    group_values <- kegg_df3[, group_samples, drop = FALSE]
    group_means <- rowMeans(group_values)
    result[group] <- group_means
  }
  
  foldchange <- result[[unique_groups[[2]]]] / result[[unique_groups[[1]]]]
  log2_ratios <- log2(foldchange)
  
  ## Add the prefix to the "log2_ratios" columnnames
  log2_ratios_column_name <- paste0("log2_ratios_", prefix)
  result[[log2_ratios_column_name]] <- log2_ratios
  
  ## Clean up and return the results
  result$ID <- NULL
  
  return(result)
}

# Function for processing unipept based EC abundance table for pathview 
process_ec_abundance <- function(ec_abundance_file, metadata) {
  ## Read MP based EC abundance table
  ec_abundance <- read.delim(ec_abundance_file, header=TRUE, sep="\t")

  ## Prepare KEGG id - EC table
  ec2ko <- as.matrix(keggLink("enzyme", "ko"))
  ec2ko <- as.data.frame(ec2ko) %>%
    rownames_to_column(var = "KEGG") %>%
    mutate(KEGG = gsub("ko.", "", KEGG),
           V1 = gsub("ec:", "", V1))

  ## Keep only EC with a KEGG id match
  ec_matched <- subset(ec_abundance, rownames(ec_abundance) %in% ec2ko$V1)

  ## Replace EC row names with KEGG id row names
  row_indices <- match(rownames(ec_matched), ec2ko$V1)
  new_row_names <- ec2ko$KEGG[row_indices]
  rownames(ec_matched) <- new_row_names

  ## Perform differential abundance analysis (ec based)
  maaslin_results_mp <- Maaslin2::Maaslin2(
    input_data = ec_matched,
    input_metadata = metadata,
    output = paste0("results/final/diff_abun/eggnog-maaslin2-", "MP"),
    fixed_effects = opt$group,
    standardize = FALSE
  )
  
  ## Filter differentially abundant features
  diff_abun_mp <- maaslin_results_mp$results %>% data.frame()
  diff_abun_mp_kegg <- diff_abun_mp$feature
  ec_matched2 <- ec_matched[diff_abun_mp_kegg, ]
  ec_matched3 <- ec_matched2 %>% mutate_all(~ (. / sum(.)) * 100)
  result_mp <- data.frame(ID = rownames(ec_matched3))
  rownames(result_mp) <- result_mp$ID

  unique_groups <- unique(metadata[[opt$group]])

  for (group in unique_groups) {
    group_samples <- metadata$SampleID[metadata[[opt$group]] == group]
    group_values_mp <- ec_matched3[, group_samples, drop = FALSE]
    group_means_mp <- rowMeans(group_values_mp)
    result_mp[group] <- group_means_mp
  }
  
  ## Calculate fold changes and log2 ratios
  foldchange_mp <- result_mp[[unique_groups[[2]]]] / result_mp[[unique_groups[[1]]]]
  log2_ratios_mp <- log2(foldchange_mp)
  
  ## Add the prefix to the "log2_ratios" column name
  log2_ratios_column_name <- paste0("log2_ratios_", "MP")
  result_mp[[log2_ratios_column_name]] <- log2_ratios_mp 
  
  ## Clean up and return the results
  result_mp$ID <- NULL
  
  return(result_mp)

}

# Function for replacing Inf values in log2ratio tables
replace_inf <- function(column) {
  column[column == Inf] <- 1
  column[column == -Inf] <- -1
  return(column)
}

# Function for performing integrated pathway analysis  
perform_pathway_analysis <- function(omics, metadata) {
  
  if (omics == 3) {
    result_mg <- process_seq("MG")
    result_mp <- process_ec_abundance("results/final/MP/ec_abundance.txt", metadata)
    common_kegg <- intersect(row.names(result_mg), row.names(result_mp))
    result_to_merge_mg <- result_mg[row.names(result_mg) %in% common_kegg, ]
    result_to_merge_mp <- result_mp[row.names(result_mp) %in% common_kegg, ]
    merged_results_common <- merge(result_to_merge_mg, 
                                   result_to_merge_mp,
                                   by = "row.names",
                                   all = FALSE)
  } else if (omics == 4) {
    result_mg <- process_seq("MG")
    result_mt <- process_seq("MT")
    common_kegg <- Reduce(intersect, list(rownames(result_mg), rownames(result_mt)))
    result_to_merge_mg <- result_mg[row.names(result_mg) %in% common_kegg, ]
    result_to_merge_mt <- result_mt[row.names(result_mt) %in% common_kegg, ]
    merged_results_common <- merge(result_to_merge_mg, 
                                   result_to_merge_mt,
                                   by = "row.names",
                                   all = FALSE)
    
  } else if (omics >= 5) {
    result_mg <- process_seq("MG")
    result_mt <- process_seq("MT")
    result_mp <- process_ec_abundance("results/final/MP/ec_abundance.txt", metadata)
    common_kegg <- Reduce(intersect, list(rownames(result_mg), rownames(result_mt), rownames(result_mp)))
    result_to_merge_mg <- result_mg[row.names(result_mg) %in% common_kegg, ]
    result_to_merge_mt <- result_mt[row.names(result_mt) %in% common_kegg, ]
    result_to_merge_mp <- result_mp[row.names(result_mp) %in% common_kegg, ]
    merged_results_common_temp1 <- merge(result_to_merge_mg,
                                   result_to_merge_mt,
                                   by = "row.names",
                                   all = FALSE)
    merged_results_common_temp2 <- column_to_rownames(merged_results_common_temp1, var = "Row.names")
    merged_results_common <- merge(merged_results_common_temp2,
                                   result_to_merge_mp,
                                   by = "row.names",
                                   all = FALSE)
  } 
  
  log2_columns <- grep("log2", colnames(merged_results_common), value = TRUE)
  subset_merged_results_common <- merged_results_common[, c("Row.names", log2_columns)] %>%
    column_to_rownames(var = "Row.names") %>%
    data.frame()
  modified_subset <- as.data.frame(lapply(subset_merged_results_common, replace_inf))
  rownames(modified_subset) <- rownames(subset_merged_results_common)
              
  pathways <- ko2path[ko2path$KEGG %in% common_kegg, "V1"] %>% unique()
  setwd("results/final/integrated/pathview")
  pv.out <- pathview(
    gene.data = modified_subset,
    pathway.id = pathways,
    species = "ko",
    out.suffix = "gnomo2",
    kegg.native = TRUE
  )
  
  return(pv.out)
}
###############################################################################

# Create output directory
dir.create("results/final/integrated/pathview")

# Read metadata table, config file and unipept and emapper result files
metadata <- read.delim("resources/metadata.txt", header=TRUE, sep="\t")
rownames(metadata) <- metadata$SampleID

# Read snakemake config file 
config <- yaml::yaml.load_file("config/config.yaml")
omics <- config$omics_combination

# Prepare KEGG id - pathway table
ko2path <- as.matrix(keggLink("pathway", "ko"))
ko2path <- as.data.frame(ko2path) %>%
  rownames_to_column(var = "KEGG")
ko2path$KEGG <- gsub("ko.", "", ko2path$KEGG)
ko2path$V1 <- gsub("path:map", "", ko2path$V1)

# Perform integrated pathway analysis
pathway_analysis_result <- perform_pathway_analysis(omics, metadata)

# List all files in the directory
files <- list.files(full.names = TRUE)

# Filter out .txt and .multi.png files
files_to_keep <- files[grep("\\.multi\\.png$", files)]

# Remove unwanted files (.xml, .png)
files_to_delete <- setdiff(files, files_to_keep)
file.remove(files_to_delete)

# Extract unique_groups[[2]] and unique_groups[[1]] values as strings
unique_groups <- unique(metadata[[opt$group]])
group_2 <- as.character(unique_groups[[2]])
group_1 <- as.character(unique_groups[[1]])

# Open a text file in write mode (this will create a new file)
file_path <- "pathview_analysis_summary.txt"
file_conn <- file(file_path, "w")

# Write explanations or summary of your analysis
writeLines("Analysis Summary:", file_conn)
writeLines(paste("1. This analysis demonstrates results for:",  group_2, "/", group_1), file_conn)
writeLines("2. The findings from each omics type are depicted separately, represented as split nodes.", file_conn)
writeLines("3. Nodes are color-coded in the sequence: MG -> MT -> MP.", file_conn)
writeLines("   In cases like MG+MP, the left node represents MG, and the right node represents MP.", file_conn)
writeLines("   Similarly, for MG+MT, MG appears on the left while MT appears on the right.", file_conn)
writeLines("   In instances of MG+MT+MP, MG is positioned on the left, MT in the middle, and MP on the right.", file_conn)

# Close the file connection
close(file_conn)