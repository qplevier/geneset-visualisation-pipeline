# Title     : limma
# Objective : Differential gene expression analysis with limma
# Created by: Quinten Plevier
# Created on: 22-2-2024

# Load libraries
library(tidyverse, quietly = TRUE)
library(phyloseq)
library(variancePartition)
library(BiocParallel)
library(edgeR)
library(openxlsx)
library(ggnewscale)
register(SnowParam(snakemake@threads[[1]], "SOCK", progressbar = T))

# Load phyloseq object
physeq <- read_rds(snakemake@input[["physeq"]])

# Get count data
countData <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  otu_table() %>%
  as.data.frame() %>%
  filter(!str_detect(rownames(.), '\\|')) %>% # remove stratified abundances (microbiota)
  filter(rownames(.) != "UNMAPPED" | rownames(.) != "UNGROUPED")

# Get meta data
modelMeta <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  sample_data() %>%
  as("data.frame") %>%
  mutate(id_sample = as.factor(id_sample),
         time = as.factor(time))

# Create the model data
modelCounts <- countData %>%
  filter(rowSums(.) != 0) %>%
  DGEList() %>% # Function reads count data and constructs an object that holds raw count data along with relevant information about the samples
  calcNormFactors() # These factors are used to adjust for differences in library sizes or sequencing depths among samples. Common normalization methods include the TMM (trimmed mean of M values) and the upper-quartile normalization

# Formulate the model formula
modelFormula <- ~time + (1 | id_sample)

# Using limma voom to calculate log2foldchanges
fittedModels <- modelCounts %>%
  voomWithDreamWeights(modelFormula, modelMeta) %>%
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

# Get all results
modelResult <- modelSummary %>%
  unnest(result) %>%
  arrange(adj.P.Val)

# Get signficant results
SigmodelResult <- modelResult %>%
  filter(adj.P.Val < 0.05)

# Write the results to a file
sheets <- list("alldata" = countData %>% rownames_to_column("enzyme"),
               "allresults" = modelResult,
               "significant-results" = SigmodelResult)
write.xlsx(sheets, snakemake@output[["dge_out"]])

# Visualise the upregulated differentials
filterData <- SigmodelResult %>%
  select(-t, -z.std) %>%
  rename(enzyme = EC) %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val) %>%
  mutate(enzyme = gsub("\\((metacyc|expasy)\\) ", "", enzyme))

selecttop <- filterData %>%
  select(enzyme, logFC) %>%
  distinct(enzyme) %>%
  slice(1:50)

all <- filterData %>%
  inner_join(selecttop)

averageEx <- all %>%
  mutate(term = gsub("time", "", term)) %>%
  select(term, enzyme, AveExpr) %>%
  mutate(typename = "AveExp", .before = term) %>%
  select(-term) %>%
  mutate(plotexpression = AveExpr)

plotDataheat <- all %>%
  mutate(term = gsub("time", "", term)) %>%
  select(term, enzyme, logFC, AveExpr) %>%
  mutate(type = term, .before = term) %>%
  mutate(typename = "logFC", .before = term) %>%
  bind_rows(averageEx) %>%
  mutate_at(c('term'), ~replace_na(., "")) %>%
  mutate(enzyme = gsub("\\(expasy)", "", enzyme))

plotDataheat %>%
  arrange(AveExpr) %>% # arrange to abundance level
  mutate(enzyme = fct_inorder(factor(enzyme, ordered = TRUE))) %>%  # Order genus to abundance level
  ggplot(aes(term, enzyme)) +
  geom_tile(aes(fill = logFC)) +
  labs(fill = "log2fc") +
  scale_fill_gradient2(low = "royalblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  new_scale("fill") +
  geom_tile(data = plotDataheat %>%
    filter(typename == "AveExp"),
            aes(fill = AveExpr)) +
  labs(fill = "Average\nexpression\n(log2)") +
  scale_fill_gradient(low = "antiquewhite", high = "aquamarine4") +
  theme_bw() +
  facet_grid(~typename,
             scales = "free_x",
             space = "free_x"
  ) +
  labs(y = "Enzymes",
       x = "Time",
       title = "Differentially expressed genes",
       subtitle = glue::glue("Upregulated in {snakemake@params[['micro']]}"),
       caption = "adj.P.Val < 0.05") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10)) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_discrete(
    label = function(x)
      stringr::str_trunc(x, 60, "center")
  )

ggsave(snakemake@output[["png_up"]],
       height = length(unique(plotDataheat$enzyme)) * 20 + 500,
       width = 1500,
       units = "px",
       dpi = 150,
       device = "png")

# Visualise the downregulated differentials
filterData <- SigmodelResult %>%
  select(-t, -z.std) %>%
  rename(enzyme = EC) %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val) %>%
  mutate(enzyme = gsub("\\((metacyc|expasy)\\) ", "", enzyme))

selecttop <- filterData %>%
  select(enzyme, logFC) %>%
  distinct(enzyme) %>%
  slice(1:50)

all <- filterData %>%
  inner_join(selecttop)

averageEx <- all %>%
  mutate(term = gsub("time", "", term)) %>%
  select(term, enzyme, AveExpr) %>%
  mutate(typename = "AveExp", .before = term) %>%
  select(-term) %>%
  mutate(plotexpression = AveExpr)

plotDataheat <- all %>%
  mutate(term = gsub("time", "", term)) %>%
  select(term, enzyme, logFC, AveExpr) %>%
  mutate(type = term, .before = term) %>%
  mutate(typename = "logFC", .before = term) %>%
  bind_rows(averageEx) %>%
  mutate_at(c('term'), ~replace_na(., "")) %>%
  mutate(enzyme = gsub("\\(expasy)", "", enzyme))

plotDataheat %>%
  arrange(AveExpr) %>% # arrange to abundance level
  mutate(enzyme = fct_inorder(factor(enzyme, ordered = TRUE))) %>%  # Order genus to abundance level
  ggplot(aes(term, enzyme)) +
  geom_tile(aes(fill = logFC)) +
  labs(fill = "log2fc") +
  scale_fill_gradient2(low = "royalblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  new_scale("fill") +
  geom_tile(data = plotDataheat %>%
    filter(typename == "AveExp"),
            aes(fill = AveExpr)) +
  labs(fill = "Average\nexpression\n(log2)") +
  scale_fill_gradient(low = "antiquewhite", high = "aquamarine4") +
  theme_bw() +
  facet_grid(~typename,
             scales = "free_x",
             space = "free_x"
  ) +
  labs(y = "Enzymes",
       x = "Time",
       title = "Differentially expressed genes",
       subtitle = glue::glue("Downregulated in {snakemake@params[['micro']]}"),
       caption = "adj.P.Val < 0.05") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10)) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_discrete(
    label = function(x)
      stringr::str_trunc(x, 60, "center")
  )

ggsave(snakemake@output[["png_down"]],
       height = length(unique(plotDataheat$enzyme)) * 20 + 500,
       width = 1500,
       units = "px",
       dpi = 150,
       device = "png")
