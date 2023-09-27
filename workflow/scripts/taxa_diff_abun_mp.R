# Title: Differential Abundance Analysis of Taxa Identified by Metaproteomics
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a differential abundance analysis of taxa identified 
#   by metaproteomics using Maaslin2.

# Load packages (installed using conda)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(optparse)

# define arguments
commandArgs(trailing=TRUE)

option_list <- list(
  make_option(c("-g", "--group"), type="character", default="Group", help="enter group name"),
  make_option(c("-t", "--taxa_rank"), type = "character", default="Genus", help = "enter taxa rank")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Get metadata, abundance and taxonomy tables
metadata <- import_qiime_sample_data("resources/metadata.txt")

abundance <- read_delim("results/final/diff_abun/MP/pep_abundance_table.txt", 
                        delim = "\t", 
                        escape_double = FALSE, 
                        trim_ws = TRUE
                        )

taxonomy <- read_csv("results/final/diff_abun/MP/unipept_results.csv",
                     col_names = FALSE, 
                     trim_ws = TRUE
                     )

# keep only peptides that have taxonomic assignment by unipept and vice versa
taxonomy$X9 <- gsub("^;", "", taxonomy$X9)
taxonomy_cleaned <- subset(taxonomy, !grepl("[;-]", X9))
taxonomy_uniq = taxonomy_cleaned[!duplicated(taxonomy_cleaned$X1),]
abundance_matched <- subset(abundance, peptide %in% taxonomy_uniq$X1)
taxonomy_matched <- subset(taxonomy_uniq, X1 %in% abundance_matched$peptide)

colnames(taxonomy_matched) <- c("peptide", "Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus", "Species", "EC"
                             )

# prepare abundance table with EC ids
taxonomy_ec <- taxonomy_matched[c("peptide","EC")]
ec_abundance <- merge(abundance_matched, taxonomy_ec, by="peptide", all=FALSE) %>% filter(EC != "")
ec_abundance$peptide <- NULL


# Define function to sum values within each role
sum_duplicates <- function(dataframe) {
  dataframe %>%
    group_by(EC) %>%
    summarise(across(where(is.numeric), sum))
}

# Apply the sum function to the merged annotated dataframe
ec_abundance2 <- sum_duplicates(ec_abundance) %>%
  column_to_rownames(var = "EC")

# save final dataframe for pathway level integrated analysis
write.table(ec_abundance2, "results/final/diff_abun/MP/ec_abundance.txt", quote = FALSE, 
            col.names=TRUE,row.names = TRUE, sep="\t")

# prepare taxonomy tables for phyloseq merge
taxonomy_noec <- taxonomy_matched[-c(9)]

#subset the table by user defined taxonomic rank
taxonomy_noec_subsetted <- taxonomy_noec[c("peptide","Kingdom")]

# prepare abundance table for phyloseq
abundance_matched2 <- column_to_rownames(abundance_matched, var = "peptide")
abundance_matched2[is.na(abundance_matched2)] <- "0"
abundance_matched2[] <- lapply(abundance_matched2, as.numeric)

# prepare phyloseq components
abundance_matched_physeq <- otu_table(abundance_matched2, 
                                      taxa_are_rows = TRUE
                                      )

taxonomy_matched_physeq <- column_to_rownames(taxonomy_noec, var = "peptide") %>%
                           as.matrix() %>%
                           tax_table()

# create a phyloseq object
pep_physeq <- merge_phyloseq(abundance_matched_physeq, metadata, taxonomy_matched_physeq)     

# remove peptides that don't have a taxonomic assignment at kingdom level
pep_physeq_2 <- pep_physeq %>% subset_taxa(Kingdom != "")

# modify taxonomy table - change "NA" assignments 
tax.clean <- as(tax_table(pep_physeq_2),"matrix") %>% as.data.frame()
tax.clean[is.na(tax.clean)] <- ""


for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = "_")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = "_")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = "_")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = "_")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = "_")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = "_")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

TAX = tax_table(as.matrix(tax.clean))

# create a new phyloseq object with modified taxonomy table
pep_physeq_3 <- phyloseq(otu_table(pep_physeq_2), TAX, sample_data(pep_physeq_2))

# write peptide physeq for further analyses 
saveRDS(pep_physeq_3, file = "results/final/diff_abun/MP/peptide_phyloseq.Rdata")

pep_physeq_4 <- tax_glom(pep_physeq_3, taxrank=opt$taxa_rank, NArm=TRUE)
taxa_names(pep_physeq_4) <- tax_table(pep_physeq_4)[,opt$taxa_rank]

otu <- as(otu_table(pep_physeq_4),"matrix") %>% as.data.frame()

# perform differential abundance analysis
maaslin_results = Maaslin2::Maaslin2(input_data = otu,
                                     input_metadata = metadata,
                                     output = "results/final/diff_abun/taxa-maaslin2-MP",
                                     fixed_effects = opt$group,
                                     standardize = FALSE
                                     )

