library(tidyverse, quietly = TRUE)
library(readxl)

dge <- read_excel(snakemake@input[["dge"]], sheet = "allresults") %>%
  data.frame() %>%
  rename("time" = "term", "RXN" = "EC") %>%
  select("time", "RXN", "logFC") %>%
  mutate(time = str_remove(time, "time"),
         RXN = str_extract(RXN, ".*?(?=:)")) %>%
  distinct(time, RXN, logFC, .keep_all = T) %>%
  na.omit()

dge[dge$time == snakemake@params[["timepoint"]], ] %>%
  select(RXN, logFC) %>%
  arrange(desc(logFC)) %>%
  write.table(
    file = snakemake@output[["out_rank"]],
    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
  )
