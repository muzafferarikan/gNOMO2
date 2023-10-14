# Title: Processing Amplicon Sequencing Data
# Author: Muzaffer Arikan
# Date: July 2023
# Description:
#   This script processes amplicon sequencing data to determine most abundant 
#   taxa present in the samples. It generates a phyloseq object and a taxonomic
#   composition plot.

# Load packages (installed using conda)
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(microbiome)
library(Maaslin2)
library(tidyverse)
library(optparse)

# Define command line arguments
commandArgs(trailing = TRUE)

option_list <- list(
  make_option(c("-g", "--group"), type = "character", default = "Group", help = "Enter group name"),
  make_option(c("-t", "--taxa_rank"), type = "character", default = "Genus", help = "Enter taxa rank"),
  make_option(c("-n", "--top_taxa"), type = "character", default = "5", help = "Enter most abundant taxa number")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Import metadata file 
metadata <- read.delim("resources/metadata.txt", header=TRUE, sep="\t")
rownames(metadata) <- metadata$SampleID


# Get merged amplicon sequencing data files 
merged_as <- list.files(
  path = "results/intermediate_files/merged/AS",
  pattern = "*extendedFrags.fastq",
  full.names = TRUE
)

# Extract sample names 
sample_names <- sub("^results/intermediate_files/merged/AS/(.*?)\\.extendedFrags\\.fastq$", "\\1", merged_as)

# Assign the filenames for the filtered files.
filts <- file.path("results/intermediate_files/dada_filtered", paste0(sample_names, "_filt.fastq.gz"))
names(filts) <- sample_names

# Apply filtering and trimming
out <- filterAndTrim(
  merged_as,
  filts,
  truncLen = 0,
  maxEE = 10,
  truncQ = 4,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# Learn error rates
errF <- learnErrors(filts, multithread = TRUE)

# Infer sequence variants
dds <- vector("list", length(sample_names))
names(dds) <- sample_names
for (sam in sample_names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err = errF, multithread = TRUE)
}

# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab_nochim, "results/intermediate_files/silva_db/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

# Create a phyloseq object
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(metadata))

# Rename ASVs
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Remove junk taxa
ps_clean <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "Mitochondria" &
      Order != "Chloroplast" &
      Phylum != "none"
  )

# Generate barplot based on relative abundances
ps_rel <- microbiome::transform(ps_clean, "compositional")
ps_rel_taxa <- aggregate_taxa(ps_rel, level = opt$taxa_rank)
topn_taxa <- top_taxa(ps_rel_taxa, n = opt$top_taxa)
ps_rel_taxa_topn <- prune_taxa(topn_taxa, ps_rel_taxa)
abun_plot <- plot_bar(ps_rel_taxa_topn, x = "SampleID", fill = opt$taxa_rank) +
                  facet_grid(as.formula(paste("~", opt$group)), scales = "free", space = "free")

# Prepare the list of most abundant taxa
topn_taxa_2 <- gsub("Other", "", gsub("Unknown", "", topn_taxa))
topn_taxa_3 <- as.character(topn_taxa_2[topn_taxa_2 != ""])

# perform differential abundance analysis
ps_clean_rank <- aggregate_taxa(ps_clean, level = opt$taxa_rank)
otu <- as(otu_table(ps_clean_rank),"matrix") %>% as.data.frame()

maaslin_results = Maaslin2::Maaslin2(input_data = otu,
                                     input_metadata = metadata,
                                     output = "results/final/diff_abun/taxa-maaslin2-AS",
                                     fixed_effects = opt$group,
                                     standardize = FALSE
                                     )

# Write output files (abundance barplot, phyloseq object, and most abundant taxa names)
svg("results/final/AS/as_abundance_plot.svg")
abun_plot
dev.off()

write.table(topn_taxa_3, "results/intermediate_files/top_taxa.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

saveRDS(ps_clean, "results/final/AS/as_physeq.rds")

print("Completed successfully.")
