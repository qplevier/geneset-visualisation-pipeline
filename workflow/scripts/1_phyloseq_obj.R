# Title     : Phyloseq object
# Objective : Use the output of the HUMAnN pipeline to create a phyloseq object with the Counts per Million
# Created by: Quinten Plevier
# Created on: 22-2-2024

# Load libraries
library(phyloseq)
library(tidyverse, quietly = TRUE)

# Load count and meta data
countData <- read.delim(snakemake@input[["countdata"]])
metaData <- read.delim(snakemake@input[["metadata"]])

# Prepare metadata for a phyloseq object
samples <- metaData %>%
  column_to_rownames("sample_name") %>%
  sample_data()

# Prepare countdata for the phyloseq object including RXN identifiers
data <- countData %>%
  rename(enzyme = X..Gene.Family) %>%
  column_to_rownames("enzyme") %>%
  rename_with(~gsub("_.*", "", .x) %>%
    gsub("\\.", "_", .)) %>%
  otu_table(taxa_are_rows = T)

# Create phyloseq object with RXN identifiers and write to a .RDS file
physeq <- phyloseq(data, samples)
write_rds(physeq, snakemake@output[["physeq_out_rxn"]])

# Prepare countdata for the phyloseq object
data <- countData %>%
  rename(enzyme = X..Gene.Family) %>%
  mutate(enzyme = str_remove(enzyme, ".*: ")) %>%
  distinct() %>%
  group_by(enzyme) %>% # Remove duplicate enzymes if available by grouping and summing their values
  summarise(across(everything(), sum)) %>%
  ungroup() %>%
  column_to_rownames("enzyme") %>%
  rename_with(~gsub("_.*", "", .x) %>%
    gsub("\\.", "_", .)) %>%
  otu_table(taxa_are_rows = T)

# Create phyloseq object and write to a .RDS file
physeq <- phyloseq(data, samples)
write_rds(physeq, snakemake@output[["physeq_out"]])
