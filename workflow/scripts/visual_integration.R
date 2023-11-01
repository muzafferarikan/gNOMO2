# Title: Visual Integration of Omics Results
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a visual integration of multi-omics data. It covers
#   data loading (omics abundance tables), preprocessing and visual integration.

# Load necessary packages
library(phyloseq)
library(combi)
library(tidyverse)
library(microbiome)
library(optparse)

commandArgs(trailing=TRUE)

option_list <- list( 
  make_option(c("-g", "--group"), type="character", default="Group", 
              help="enter group name"))

opt <- parse_args(OptionParser(option_list=option_list))

dir.create("results/final/integrated")

# Read metadata table, config file and unipept and emapper result files
metadata <- read.table("resources/metadata.txt", header=TRUE, sep = "\t", fileEncoding = "UTF-8")
metadata2 <- column_to_rownames(metadata, var = "SampleID")
#metadata2$DiseaseStage <- NULL
metadata2$Season <- as.factor(metadata2$Season)
#metadata2$Stage <- NULL
metadata3 <- sample_data(metadata2)

# Read snakemake config file
config <- yaml::yaml.load_file("config/config.yaml")
omics <- config$omics_combination

if (omics == 2) {
  # Read MG taxa abundance table
  as_physeq <- readRDS("results/final/AS/as_physeq.rds")
  as_physeq_genus <- aggregate_taxa(as_physeq, level = "Genus")
  
  # Read MP taxa abundance table
  ec_abundance <- read.delim("results/final/MP/pep_abundance_table.txt", header=TRUE, sep="\t")
  ec_abundance <- column_to_rownames(ec_abundance, var = "peptide")
  ec_abundance2 <- otu_table(ec_abundance, taxa_are_rows = TRUE)

  # merge mp phyloseq object
  mp_physeq <- merge_phyloseq(ec_abundance2, metadata3)
  
  # Perform visual integration by combi analysis
  microMetaboIntConstr <- combi(
    list(microbiome = as_physeq_genus, metaproteomics = mp_physeq),
    distributions = c("quasi", "gaussian"),
    compositional = c(TRUE, FALSE),
    logTransformGaussian = FALSE,
    covariates = metadata2,
    verbose = TRUE,
    prevCutOff = 0.01,
    allowMissingness = TRUE
  )
  
} else if (omics == 3) {
  # Read MG taxa abundance table
  mg_taxa_abundance <- read.delim("results/final/MG/mg_taxa_abundance.txt", header=TRUE, sep="\t")
  mg_taxa_abundance <- column_to_rownames(mg_taxa_abundance, var = "taxa")
  mg_taxa_abundance[is.na(mg_taxa_abundance)] <- 0
  mg_taxa_abundance2 <- otu_table(mg_taxa_abundance, taxa_are_rows = TRUE)
  
  # Read MP taxa abundance table
  ec_abundance <- read.delim("results/final/MP/pep_abundance_table.txt", header=TRUE, sep="\t")
  ec_abundance <- column_to_rownames(ec_abundance, var = "peptide")
  ec_abundance2 <- otu_table(ec_abundance, taxa_are_rows = TRUE)
  
  # merge phyloseq objects
  mg_physeq <- merge_phyloseq(mg_taxa_abundance2, metadata3)
  mp_physeq <- merge_phyloseq(ec_abundance2, metadata3)
  
  # Perform visual integration by combi analysis
  microMetaboIntConstr <- combi(
    list(microbiome = mg_physeq, metaproteomics = mp_physeq),
    distributions = c("quasi", "gaussian"),
    compositional = c(TRUE, FALSE),
    logTransformGaussian = FALSE,
    covariates = metadata2,
    verbose = TRUE,
    prevCutOff = 0.1,
    allowMissingness = TRUE
  )
  
} else if (omics == 4) {
  # Read MG taxa abundance table
  mg_taxa_abundance <- read.delim("results/final/MG/mg_taxa_abundance.txt", header=TRUE, sep="\t")
  mg_taxa_abundance <- column_to_rownames(mg_taxa_abundance, var = "taxa")
  mg_taxa_abundance[is.na(mg_taxa_abundance)] <- 0
  mg_taxa_abundance2 <- otu_table(mg_taxa_abundance, taxa_are_rows = TRUE)
  
  # Read MT taxa abundance table
  mt_taxa_abundance <- read.delim("results/final/MT/mt_taxa_abundance.txt", header=TRUE, sep="\t")
  mt_taxa_abundance <- column_to_rownames(mt_taxa_abundance, var = "taxa")
  mt_taxa_abundance[is.na(mt_taxa_abundance)] <- 0
  mt_taxa_abundance2 <- otu_table(mt_taxa_abundance, taxa_are_rows = TRUE)
  
  # merge phyloseq objects
  mg_physeq <- merge_phyloseq(mg_taxa_abundance2, metadata3)
  mt_physeq <- merge_phyloseq(mt_taxa_abundance2, metadata3)
  
  # Perform visual integration by combi analysis
  microMetaboIntConstr <- combi(
    list(microbiome = mg_physeq, microbiome = mt_physeq),
    distributions = c("quasi", "quasi"),
    compositional = c(TRUE, TRUE),
    covariates = metadata2,
    verbose = TRUE,
    prevCutOff = 0.1,
    allowMissingness = TRUE
  )
  
} else if (omics >= 5) {
  # Read MG taxa abundance table
  mg_taxa_abundance <- read.delim("results/final/MG/mg_taxa_abundance.txt", header=TRUE, sep="\t")
  mg_taxa_abundance <- column_to_rownames(mg_taxa_abundance, var = "taxa")
  mg_taxa_abundance[is.na(mg_taxa_abundance)] <- 0
  mg_taxa_abundance2 <- otu_table(mg_taxa_abundance, taxa_are_rows = TRUE)
  
  # Read MT taxa abundance table
  mt_taxa_abundance <- read.delim("results/final/MT/mt_taxa_abundance.txt", header=TRUE, sep="\t")
  mt_taxa_abundance <- column_to_rownames(mt_taxa_abundance, var = "taxa")
  mt_taxa_abundance[is.na(mt_taxa_abundance)] <- 0
  mt_taxa_abundance2 <- otu_table(mt_taxa_abundance, taxa_are_rows = TRUE)
  
  # Read MP taxa abundance table
  ec_abundance <- read.delim("results/final/MP/pep_abundance_table.txt", header=TRUE, sep="\t")
  ec_abundance <- column_to_rownames(ec_abundance, var = "peptide")
  ec_abundance2 <- otu_table(ec_abundance, taxa_are_rows = TRUE)
  
  # merge phyloseq objects
  mg_physeq <- merge_phyloseq(mg_taxa_abundance2, metadata3)
  mt_physeq <- merge_phyloseq(mt_taxa_abundance2, metadata3)
  mp_physeq <- merge_phyloseq(ec_abundance2, metadata3)
  
  # Perform visual integration by combi analysis
  microMetaboIntConstr <- combi(
    list(microbiome = mg_physeq, microbiome = mt_physeq, metaproteomics = mp_physeq),
    distributions = c("quasi", "quasi", "gaussian"),
    compositional = c(TRUE, TRUE, FALSE),
    logTransformGaussian = FALSE,
    covariates = metadata2,
    verbose = TRUE,
    prevCutOff = 0.1,
    allowMissingness = TRUE
  )
  
}

# Set up graphics device for PNG output
svg("results/final/integrated/combi_plot.svg",  width=10, height=10)
plot(microMetaboIntConstr)
dev.off()
