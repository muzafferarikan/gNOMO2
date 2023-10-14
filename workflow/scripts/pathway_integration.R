# Title: Pathway Level Integrated Analysis
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a pathway level integrated analysis of multi-omics
#   data. It covers data loading (emapper results), preprocessing, differential
#   abundance analysis and pathway based data integration and visualization.


# Load necessary packages
library(phyloseq)
library(tidyverse)
library(optparse)
library(pathview)
library(KEGGREST)
library(gtools)

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

# Function for processing emapper results and output tables for pathview  
process_seq <- function(prefix) {
  
  ## Read emapper result files
  eggnog_files <- list.files(
    path = "results/intermediate_files/eggnog",
    pattern = "*eggnog.emapper.annotations",
    full.names = TRUE
  )
  
  ## Import emapper result files
  eggnog_data <- lapply(eggnog_files, function(file_path) {
    read.table(
      file_path,
      sep = "\t",
      header = FALSE,
      skip = 4,
      quote = ""
    )
  })
  
  ## Create a df list of emapper results
  emapper_df_list <- setNames(eggnog_data, eggnog_files)
  
  ## Create a new list of data frames without empty cells
  emapper_df_list_filled <- lapply(emapper_df_list, fill_empty_cells)
  
  ## Extract coverage values and merge identical KEGG IDs
  emapper_df_list_clean <- lapply(emapper_df_list_filled, function(dataframe) {
    dataframe %>%
      mutate(V1 = as.numeric(str_remove_all(V1, ".*cov_|\\.g|_.*"))) %>%
      group_by(V7) %>%
      summarise(summed_value = sum(V1)) %>%
      replace_na(list(summed_value = 0))
  })
  
  emapper_df_final_list <- emapper_df_list_clean[grepl(prefix, names(emapper_df_list_clean))]
  
  merged_df <- emapper_df_final_list %>%
    reduce(full_join, by = "V7") %>%
    column_to_rownames(var = "V7")
  
  pats <- paste0("results/intermediate_files/eggnog/", prefix, "_|_eggnog.emapper.annotations")
  temp <- list.files(path = "results/intermediate_files/eggnog", pattern = paste0(prefix, ".*eggnog\\.emapper\\.annotations"), full.names = TRUE)
  sample_names <- str_replace_all(temp, pats, "")
  
  colnames(merged_df) <- sample_names
  
  merged_df2 <- rownames_to_column(merged_df, var = "KEGG")
  merged_df2[is.na(merged_df2)] <- 0
  merged_df_uniq <- merged_df2[!grepl(",", merged_df2$KEGG), ]
  rownames(merged_df_uniq) <- NULL
  
  kegg_df <- column_to_rownames(merged_df_uniq, var = "KEGG")
  
  maaslin_results <- Maaslin2::Maaslin2(
    input_data = kegg_df,
    input_metadata = metadata,
    output = paste0("results/final/diff_abun/eggnog-maaslin2-", prefix),
    fixed_effects = opt$group,
    standardize = FALSE
  )
  
  diff_abun <- maaslin_results$results %>% data.frame() %>% filter(qval < 1)
  diff_abun_kegg <- diff_abun$feature
  kegg_df2 <- kegg_df[diff_abun_kegg, ]
  result <- data.frame(ID = rownames(kegg_df2))
  rownames(result) <- result$ID
  result$ID <- NULL
  
  unique_groups <- unique(metadata[[opt$group]])
  
  for (group in unique_groups) {
    group_samples <- metadata$SampleID[metadata[[opt$group]] == group]
    group_values <- kegg_df2[, group_samples, drop = FALSE]
    group_means <- rowMeans(group_values)
    result[group] <- group_means
  }
  
  foldchange <- result[[unique_groups[[2]]]] / result[[unique_groups[[1]]]]
  log2_ratios <- log2(foldchange)
  
  ## Add the prefix to the "log2_ratios" column name
  log2_ratios_column_name <- paste0("log2_ratios_", prefix)
  result[[log2_ratios_column_name]] <- log2_ratios  
  
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

  ## Replacing EC row names with KEGG id row names
  row_indices <- match(rownames(ec_matched), ec2ko$V1)
  new_row_names <- ec2ko$KEGG[row_indices]
  rownames(ec_matched) <- new_row_names

  ## Maaslin analysis for MP
  maaslin_results_mp <- Maaslin2::Maaslin2(
    input_data = ec_matched,
    input_metadata = metadata,
    output = paste0("results/final/diff_abun/eggnog-maaslin2-", "MP"),
    fixed_effects = opt$group,
    standardize = FALSE
  )
  
  ## Filtering differentially abundant features
  diff_abun_mp <- maaslin_results_mp$results %>% data.frame() %>% filter(qval < 1)
  diff_abun_mp_kegg <- diff_abun_mp$feature
  ec_matched2 <- ec_matched[diff_abun_mp_kegg, ]
  result_mp <- data.frame(ID = rownames(ec_matched2))
  rownames(result_mp) <- result_mp$ID
  result_mp$ID <- NULL

  unique_groups <- unique(metadata[[opt$group]])

  for (group in unique_groups) {
    group_samples <- metadata$SampleID[metadata[[opt$group]] == group]
    group_values_mp <- ec_matched2[, group_samples, drop = FALSE]
    group_means_mp <- rowMeans(group_values_mp)
    result_mp[group] <- group_means_mp
  }
  
  foldchange_mp <- result_mp[[unique_groups[[2]]]] / result_mp[[unique_groups[[1]]]]
  log2_ratios_mp <- log2(foldchange_mp)
  
  ## Add the prefix to the "log2_ratios" column name
  log2_ratios_column_name <- paste0("log2_ratios_", "MP")
  result_mp[[log2_ratios_column_name]] <- log2_ratios_mp 

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
    common_kegg <- intersect(row.names(result_mg), row.names(result_mt))
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
    merged_results_common <- merge(result_to_merge_mg,
                                   result_to_merge_mt,
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

sink("results/final/integrated/pathview/log.txt", append = TRUE)


# Perform integrated pathway analysis
pathway_analysis_result <- perform_pathway_analysis(omics, metadata)

sink()