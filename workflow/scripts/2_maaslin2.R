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
data <- otu_table(physeq) %>%
  as.data.frame() %>%
  rownames_to_column("enzymes") %>%
  filter(!grepl("\\|", enzymes)) %>%
  filter(enzymes != "UNMAPPED" & enzymes != "UNGROUPED") %>%
  mutate(enzymes = str_remove(enzymes, ".*: "),
         enzymes = str_remove(enzymes, "\\((metacyc|expasy)\\) ")) %>%
  separate(enzymes, c("enzymes", "ec"), sep = "(?=\\[\\d\\..*\\])", fill = "right") %>%
  distinct(across(-enzymes), .keep_all = TRUE) %>%
  unite("enzymes", enzymes, ec, sep = "") %>%
  mutate(enzymes = make.unique(enzymes)) %>%
  data.table::transpose(keep.names = "enzymes") %>%
  column_to_rownames("enzymes") %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  filter(grepl(ifelse(snakemake@params[["micro"]] == "healthy", "Pala", "NS"), rownames(.)))

# Execute MaAsLin2
# Filter on minumum prevelance 0.1 (default)
fit <- data %>%
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
    reference = "time,T0",
    cores = 1
  )