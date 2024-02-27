# Title     : GSEA input
# Objective : Create input files for GSEA
# Created by: Quinten Plevier
# Created on: 27-2-2024

# Load libraries
library(tidyverse, quietly = TRUE)
library(readxl)

# Read differential gene expression output, preprocess and write to .rnk file
read_excel(snakemake@input[["dge"]], sheet = "allresults") %>%
  data.frame() %>%
  rename("time" = "term", "RXN" = "EC") %>%
  select(time, RXN, logFC) %>%
  mutate(time = str_remove(time, "time"),
         RXN = str_extract(RXN, ".*?(?=:)")) %>%
  filter(time == snakemake@params[["timepoint"]]) %>%
  select(RXN, logFC) %>%
  distinct(RXN, logFC, .keep_all = T) %>%
  na.omit() %>%
  arrange(desc(logFC)) %>%
  write.table(
    file = snakemake@output[["out_rank"]],
    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
  )
