# Title     : MaAsLin2
# Objective : Differential gene expression analysis with MaAsLin2
# Created by: Quinten Plevier
# Created on: 26-2-2024

# Load libraries
library(Maaslin2)
library(phyloseq)
library(tidyverse, quietly = TRUE)

# Load phyloseq object
physeq <- read_rds(snakemake@input[["physeq"]])

# Get meta data
metadata <- sample_data(physeq) %>%
  data.frame() %>%
  mutate(micro_time = sprintf("%s_%s", micro, time),
         id_sample = str_trim(id_sample))

# Get count data
data <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  otu_table() %>%
  as.data.frame() %>%
  filter(rowSums(.) != 0) %>%
  rownames_to_column("enzymes") %>%
  filter(!grepl("\\|", enzymes)) %>%
  filter(enzymes != "UNMAPPED" | enzymes != "UNGROUPED") %>%
  mutate(enzymes = str_remove(enzymes, ".*\\((metacyc|expasy)\\) ")) %>%
  group_by(enzymes) %>%
  summarise(across(where(is.numeric), mean)) %>%
  data.table::transpose(keep.names = "enzymes") %>%
  column_to_rownames("enzymes") %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  mutate_all(function(x) as.numeric(as.character(x)))

# Execute MaAsLin2
# Filter on minumum prevelance 0.1 (default)
fit_data <- data %>%
  Maaslin2(
    input_data = .,
    input_metadata = metadata,
    output = snakemake@output[["out_dir_maaslin"]],
    # min_prevalence = 0,
    normalization = "NONE", # TSS, CLR, NONE
    transform = snakemake@params[["transformation"]],
    analysis_method = snakemake@params[["analysis_method"]],
    max_significance = 0.05,
    random_effects = "id_sample",
    fixed_effects = "time",
    reference = glue::glue("time,{snakemake@params[['baseline']]}"),
    cores = snakemake@threads
  )