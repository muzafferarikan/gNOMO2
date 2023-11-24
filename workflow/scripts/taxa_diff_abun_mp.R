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
library(microbiome)
library(optparse)

# Define command line arguments
commandArgs(trailing=TRUE)

option_list <- list(
  make_option(c("-g", "--group"), type="character", default="Group", help="enter group name"),
  make_option(c("-t", "--taxa_rank"), type = "character", default="Genus", help = "enter taxa rank"),
  make_option(c("-n", "--top_taxa"), type = "character", default = "5", help = "Enter most abundant taxa number")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Get metadata, abundance and taxonomy tables
metadata <- import_qiime_sample_data("resources/metadata.txt")

abundance <- read_delim("results/final/MP/pep_abundance_table.txt", 
                        delim = "\t", 
                        escape_double = FALSE, 
                        trim_ws = TRUE
                        )

taxonomy <- read_csv("results/final/MP/unipept_results.csv",
                     col_names = FALSE, 
                     trim_ws = TRUE
                     )

# Keep only peptides that have taxonomic assignment by unipept and vice versa
taxonomy$X9 <- gsub("^;", "", taxonomy$X9)
taxonomy_cleaned <- subset(taxonomy, !grepl("[;-]", X9))
taxonomy_uniq = taxonomy_cleaned[!duplicated(taxonomy_cleaned$X1),]
abundance_matched <- subset(abundance, peptide %in% taxonomy_uniq$X1)
taxonomy_matched <- subset(taxonomy_uniq, X1 %in% abundance_matched$peptide)

colnames(taxonomy_matched) <- c("peptide", "Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus", "Species", "EC"
                             )

# Prepare abundance table with EC ids
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

# Save final dataframe for pathway level integrated analysis
write.table(ec_abundance2, "results/final/MP/ec_abundance.txt", quote = FALSE, 
            col.names=TRUE,row.names = TRUE, sep="\t")

# Prepare taxonomy table for phyloseq merge
taxonomy_noec <- taxonomy_matched[-c(9)] 

# Prepare abundance table for phyloseq
abundance_matched2 <- column_to_rownames(abundance_matched, var = "peptide")
abundance_matched2[is.na(abundance_matched2)] <- "0"
abundance_matched2[] <- lapply(abundance_matched2, as.numeric)

# Prepare phyloseq components
abundance_matched_physeq <- otu_table(abundance_matched2, 
                                      taxa_are_rows = TRUE
                                      )

taxonomy_matched_physeq <- column_to_rownames(taxonomy_noec, var = "peptide") %>%
                           as.matrix() %>%
                           tax_table()

# Create a phyloseq object
pep_physeq <- merge_phyloseq(abundance_matched_physeq, metadata, taxonomy_matched_physeq)     
pep_physeq_2 <- pep_physeq %>% subset_taxa(opt$taxa_rank != "")

# Generate barplot based on relative abundances
ps_rel <- microbiome::transform(pep_physeq_2, "compositional")
ps_rel_taxa <- aggregate_taxa(ps_rel, level = opt$taxa_rank)
topn_taxa <- top_taxa(ps_rel_taxa, n = opt$top_taxa)
ps_rel_taxa_topn <- prune_taxa(topn_taxa, ps_rel_taxa)
abun_plot <- plot_bar(ps_rel_taxa_topn, x = "SampleID", fill = opt$taxa_rank) +
                  facet_grid(as.formula(paste("~", opt$group)), scales = "free", space = "free")

# Save abundance barplot
svg("results/final/MP/mp_abundance_plot.svg")
abun_plot
dev.off()

# Write peptide physeq for further analyses 
saveRDS(pep_physeq_2, file = "results/final/MP/peptide_phyloseq.Rdata")

# Prepare abundance table for maaslin2 analysis
pep_physeq_3 <- aggregate_taxa(pep_physeq_2, level=opt$taxa_rank)
otu <- as(otu_table(pep_physeq_3),"matrix") %>% as.data.frame()

# Perform differential abundance analysis
maaslin_results = Maaslin2::Maaslin2(input_data = otu,
                                     input_metadata = metadata,
                                     output = "results/final/diff_abun/taxa-maaslin2-MP",
                                     fixed_effects = opt$group,
                                     standardize = FALSE
                                     )

