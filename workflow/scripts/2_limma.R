# Title     : TODO
# Objective : TODO
# Created by: Quinten Plevier
# Created on: 22-2-2024

# Load libraries
library(tidyverse, quietly = TRUE)
library(phyloseq)
library(variancePartition)
library(BiocParallel)
library(edgeR)
library(openxlsx)
register(SnowParam(1, "SOCK", progressbar = T))

physeq <- read_rds(snakemake@input[["physeq"]])

countData <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  otu_table %>%
  as.data.frame() %>%
  rownames_to_column("enzyme") %>%
  filter(!str_detect(enzyme, '\\|')) %>% # remove stratified abundances (microbiota)
  column_to_rownames("enzyme")

sampleData <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  sample_data %>%
  as("data.frame") %>%
  mutate(id_sample = as.factor(id_sample),
         time = as.factor(time))

modelMeta <- sampleData

modelCounts <- countData %>%
  .[rowSums(.) != 0,] %>%
  DGEList() %>% # function reads count data and constructs an object that holds raw count data along with relevant information about the samples
  calcNormFactors() # These factors are used to adjust for differences in library sizes or sequencing depths among samples. Common normalization methods include the TMM (trimmed mean of M values) and the upper-quartile normalization

modelFormula <- ~time + (1 | id_sample)

vobjDream <- modelCounts %>%
  voomWithDreamWeights(modelFormula, modelMeta)

fittedModels <- vobjDream %>%
  dream(modelFormula, modelMeta)

design_matrix <- fittedModels$design

terms <- colnames(design_matrix)[2:ncol(design_matrix)]

modelSummary <- tibble(term = terms) %>%
  mutate(result = map(term,
                      ~topTable(fittedModels,
                                coef = .x,
                                number = "all")),
         result = map(result,
                      as_tibble,
                      rownames = "EC"))

modelResult <- modelSummary %>%
  unnest(result) %>%
  arrange(adj.P.Val)

SigmodelResult <- modelSummary %>%
  unnest(result) %>%
  arrange(adj.P.Val) %>%
  filter(adj.P.Val < 0.05)

sheets <- list("alldata" = countData %>% rownames_to_column("enzyme"),
               "allresults" = modelResult,
               "significant-results" = SigmodelResult)
write.xlsx(sheets, snakemake@output[["dge_out"]])