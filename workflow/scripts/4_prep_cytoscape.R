# Title     : Prepare Cytoscape
# Objective : Prepare for Cytoscape EnrichmentMap by creating a expression file of a subset of samples
# Created by: Quinten Plevier
# Created on: 27-2-2024

# Load libraries
library(tidyverse, quietly = TRUE)
library(readxl)
library(phyloseq)

# Preprocess count data and write out to expression file
read_rds(snakemake@input[["physeq"]]) %>%
  subset_samples(micro == snakemake@params[["micro"]] &
                   (time == snakemake@params[["baseline"]] |
                     time == snakemake@params[["timepoint"]])) %>%
  otu_table() %>%
  as.data.frame() %>%
  {num_columns <<- ncol(.); .} %>%
  rownames_to_column("gene") %>%
  filter(grepl("\\|", gene)) %>%
  mutate(gene = str_extract(gene, ".*?(?=:)")) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>%
  mutate(description = gene) %>%
  select(gene, description, everything()) %>%
  write.table(snakemake@output[["expression"]],
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

data.frame()

glue::glue(
  "{num_columns} 2 1
  # {snakemake@params[['baseline']]} {snakemake@params[['timepoint']]}
  {toString(rep(c(snakemake@params[['baseline']], snakemake@params[['timepoint']]), num_columns / 2)) %>% str_remove_all(., ',')}"
) %>%
  writeLines(snakemake@output[["classes"]])