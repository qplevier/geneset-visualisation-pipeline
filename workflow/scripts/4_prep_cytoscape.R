library(tidyverse)
library(readxl)
library(phyloseq)


countData <- read_rds(snakemake@input[["physeq"]]) %>%
  subset_samples(micro == snakemake@params[["micro"]] &
                   (time == snakemake@params[["baseline"]] |
                     time == snakemake@params[["timepoint"]])) %>%
  otu_table() %>%
  as.data.frame() %>%
  rownames_to_column("enzyme") %>%
  filter(grepl("\\|", enzyme)) %>%
  mutate(enzyme = str_extract(enzyme, ".*?(?=:)")) %>%
  na.omit() %>%
  group_by(enzyme) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("enzyme") %>%
  write.table(snakemake@output[["expression"]],
              sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
