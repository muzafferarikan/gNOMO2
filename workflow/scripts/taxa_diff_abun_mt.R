# Title: Differential Abundance Analysis of Taxa Identified by Metagenomics
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script performs a differential abundance analysis of taxa identified 
#   by metatranscriptomics using Maaslin2.

# Load packages (installed using conda)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(microbiome)
library(optparse)

commandArgs(trailing=TRUE)

option_list <- list(
  make_option(c("-g", "--group"), type="character", default="Group", help="enter group name"),
  make_option(c("-t", "--taxa_rank"), type = "character", default="Genus", help = "enter taxa rank"),
  make_option(c("-n", "--top_taxa"), type = "character", default = "5", help = "Enter most abundant taxa number")
)

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

# sum abundances for same taxons
sum_abun_by_taxon <- function(df) {
  df %>%
    group_by(taxon_name) %>%
    summarize(total_abun = sum(reads))
}

result_df_list <- lapply(df_list, sum_abun_by_taxon)

# merge all dataframes in the list by first column (taxon_name)
merged_df <- result_df_list %>% reduce(full_join, by='taxon_name') %>% 
                                column_to_rownames(var = "taxon_name")

# remove unnecessary strings and get sample names
pats = c("results/intermediate_files/kaiju/kaiju_output/MT/MT_|_kaiju_summary.tsv")
temp_samplenames <- str_replace_all(temp, pats, "")

# rename columns with sample names
names(merged_df) <- c(temp_samplenames)

# remove "unclassified" and viral sequences
merged_df <- merged_df[!(rownames(merged_df) == "unclassified"), ]
merged_df <- merged_df[!(rownames(merged_df) == "Viruses"), ]
#merged_df <- merged_df[!(rownames(merged_df) == "cannot be assigned to a (non-viral) genus"), ]
rownames(merged_df)[rownames(merged_df) == "cannot be assigned to a (non-viral) genus"] <- "Unassigned_genus"

#perform maaslin2 analysis
maaslin_results = Maaslin2::Maaslin2(input_data = merged_df,
                                     input_metadata = metadata,
                                     output = "results/final/diff_abun/taxa-maaslin2-MT",
                                     fixed_effects = opt$group,
                                     standardize = FALSE
                                     )

# prepare phyloseq components
merged_df[is.na(merged_df)] <- "0"
merged_df[] <- lapply(merged_df, as.numeric)
abundance_physeq <- otu_table(merged_df, taxa_are_rows = TRUE)

# prepare a taxonomy table for phyloseq merge
taxonomy <- data.frame(taxa = rownames(merged_df))
rownames(taxonomy) <- taxonomy$taxa
taxonomy_physeq <- taxonomy %>% as.matrix() %>% tax_table()

# create a phyloseq object
pep_physeq <- merge_phyloseq(abundance_physeq, metadata, taxonomy_physeq)     

# Generate barplot based on relative abundances
ps_rel <- microbiome::transform(pep_physeq, "compositional")

topn_taxa <- top_taxa(ps_rel, n = opt$top_taxa)
ps_rel_topn <- prune_taxa(topn_taxa, ps_rel)
abun_plot <- plot_bar(ps_rel_topn, x = "SampleID", fill="taxa") +
                  facet_grid(as.formula(paste("~", opt$group)), scales = "free", space = "free")

# Write output files (abundance barplot, phyloseq object, and most abundant taxa names)
svg("results/final/MT/mt_abundance_plot.svg")
abun_plot
dev.off()

# prepare abundance table for dowstream analysis
merged_df2 <- rownames_to_column(merged_df, var="taxa")
merged_df2[is.na(merged_df2)] <- 0
merged_df2$taxa <- paste(merged_df2$taxa, "_MT", sep = "")

# save abundance table
write.table(merged_df2,
            file = "results/final/MT/mt_taxa_abundance.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
            )