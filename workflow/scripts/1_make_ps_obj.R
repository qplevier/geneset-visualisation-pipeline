 # ---
# title: "Create the phyloseq object for the other scripts"
# author: "Quinten Plevier"
# date: "Last compiled on `r format(Sys.time(), '%d-%m-%Y')`"
# output:
#   html_document: default
# editor_options:
#   chunk_output_type: console
# ---
### Load libraries
library(phyloseq)
library(tidyverse, quietly = TRUE)
# library(tidyr)

### Load count and meta data
countData <- read.delim(snakemake@input[["countdata"]])
metaData <- read.delim(snakemake@input[["metadata"]])

### Prepare metadata for a phyloseq object
samples <- metaData %>% 
  column_to_rownames("sample_name") %>%
  sample_data()

### Prepare countdata for a phyloseq object
data <- countData %>% 
  rename(enzyme = X..Gene.Family) %>% 
  column_to_rownames("enzyme") %>% 
  rename_with(~ gsub("_.*", "", .x) %>% 
                gsub("\\.", "_", .)) %>% 
  otu_table(taxa_are_rows = T)

### Create phyloseq object and write to a .RDS file
physeq <- phyloseq(data, samples)
write_rds(physeq, snakemake@output[["physeq_out_rxn"]])

data <- countData %>%
  rename(enzyme = X..Gene.Family) %>%
  mutate(enzyme = str_remove(enzyme, ".*: ")) %>%
  distinct() %>%
  group_by(enzyme) %>%
  summarise(across(everything(), sum)) %>%
  ungroup() %>%
  column_to_rownames("enzyme") %>%
  rename_with(~ gsub("_.*", "", .x) %>%
                gsub("\\.", "_", .)) %>%
  otu_table(taxa_are_rows = T)

### Create phyloseq object and write to a .RDS file
physeq <- phyloseq(data, samples)
write_rds(physeq, snakemake@output[["physeq_out"]])


### Session info
# sessioninfo::package_info("attached")