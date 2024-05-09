# Title     : MaAsLin2 figures
# Objective : Create figures of the MaAsLin2 output
# Created by: Quinten Plevier
# Created on: 12-4-2024 10:00

# Load libraries
library(phyloseq)
library(tidyverse, quietly = TRUE)
library(ggnewscale)

# Load phyloseq object
physeq <- read_rds(snakemake@input[["physeq"]])

countData <- physeq %>%
  subset_samples(micro == snakemake@params[["micro"]]) %>%
  otu_table() %>%
  as.data.frame() %>%
  filter(!grepl("\\|", rownames(.)) &
           rowSums(.) != 0) %>%
  rownames_to_column("enzyme") %>%
  mutate(enzyme = str_remove(enzyme, ".*\\((expasy|metacyc)\\) ")) %>%
  pivot_longer(cols = -enzyme) %>%
  group_by(enzyme) %>%
  summarise(value = mean(log2(value + 1))) %>%
  distinct() %>%
  rename(AveExpr = value) %>%
  mutate(enzyme = str_replace_all(enzyme, "[^[:alnum:]]", "."))

table <- read_delim(snakemake@input[["results"]], delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  as.data.frame() %>%
  rename(enzyme = feature,
         term = value,
         logFC = coef) %>%
  select(enzyme, term, logFC, pval, qval)

fit_table <- table %>%
  filter(qval < 0.05,
         logFC > 0) %>%
  arrange(qval) %>%
  inner_join(
    .,
    table %>%
      filter(qval < 0.05,
             logFC > 0) %>%
      arrange(qval) %>%
      select(enzyme) %>%
      distinct(enzyme) %>%
      slice_head(n = 50)
  ) %>%
  select(-pval, -qval) %>%
  mutate(enzyme = str_remove(enzyme, "X(?=([^[:alnum:]]|\\d))")) %>%
  left_join(countData, by = "enzyme") %>%
  pivot_longer(cols=c("logFC", "AveExpr")) %>%
  mutate(term = ifelse(name == "AveExpr", "", term)) %>%
  left_join(countData, by = "enzyme") %>%
  arrange(AveExpr) %>%
  mutate(enzyme = str_replace(enzyme, "\\.(\\d+\\.\\d+\\.\\d+\\.\\d+)\\.", "\\[\\1\\]"),
         enzyme = str_replace(enzyme, "\\[(\\d+\\.\\d+\\.\\d+\\.\\d+)\\]\\.(\\d+\\.\\d+\\.\\d+\\.\\d+)\\.", "\\[\\1; \\2\\]"),
         enzyme = str_replace_all(enzyme, "\\.(?=[^\\[\\]]*\\[)", " "),
         enzyme = str_trim(enzyme, side = "both") %>%
           str_replace_all(., " +", " "),
         enzyme = fct_inorder(factor(enzyme, ordered = TRUE)))

fit_table %>%
  ggplot(aes(term, enzyme)) +
  geom_tile(data = fit_table %>% filter(name == "logFC"), aes(fill = value)) +
  labs(fill = "log2fc") +
  scale_fill_gradient2(low = "royalblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits = c(0, max(fit_table %>% filter(name == "logFC") %>% .$value))) +
  new_scale("fill") +
  geom_tile(data = fit_table %>%
    filter(name == "AveExpr"),
            aes(fill = value)) +
  labs(fill = "Average\nexpression\n(log2)") +
  scale_fill_gradient(low = "antiquewhite", high = "aquamarine4") +
  theme_bw() +
  facet_grid(~name,
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
       height = length(unique(fit_table$enzyme)) * 20 + 600,
       width = 1500,
       units = "px",
       dpi = 150,
       device = "png")

fit_table <- table %>%
  filter(qval < 0.05,
         logFC < 0) %>%
  arrange(qval) %>%
  inner_join(
    .,
    table %>%
      filter(qval < 0.05,
             logFC < 0) %>%
      arrange(qval) %>%
      select(enzyme) %>%
      distinct(enzyme) %>%
      slice_head(n = 50)
  ) %>%
  select(-pval, -qval) %>%
  mutate(enzyme = str_remove(enzyme, "X(?=([^[:alnum:]]|\\d))")) %>%
  left_join(countData, by = "enzyme") %>%
  pivot_longer(cols=c("logFC", "AveExpr")) %>%
  mutate(term = ifelse(name == "AveExpr", "", term)) %>%
  left_join(countData, by = "enzyme") %>%
  arrange(AveExpr) %>%
  mutate(enzyme = str_replace(enzyme, "\\.(\\d+\\.\\d+\\.\\d+\\.\\d+)\\.", "\\[\\1\\]"),
         enzyme = str_replace(enzyme, "\\[(\\d+\\.\\d+\\.\\d+\\.\\d+)\\]\\.(\\d+\\.\\d+\\.\\d+\\.\\d+)\\.", "\\[\\1;\\2\\]"),
         enzyme = str_replace_all(enzyme, "\\.(?=[^\\[\\]]*\\[)", " "),
         enzyme = str_trim(enzyme, side = "both") %>%
           str_replace_all(., " +", " "),
         enzyme = fct_inorder(factor(enzyme, ordered = TRUE)))

fit_table %>%
  ggplot(aes(term, enzyme)) +
  geom_tile(data = fit_table %>%
    filter(name == "logFC"),
            aes(fill = value)) +
  labs(fill = "log2fc") +
  scale_fill_gradient2(low = "royalblue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits = c(min(fit_table %>% filter(name == "logFC") %>% .$value), 0)) +
  new_scale("fill") +
  geom_tile(data = fit_table %>%
    filter(name == "AveExpr"),
            aes(fill = value)) +
  labs(fill = "Average\nexpression\n(log2)") +
  scale_fill_gradient(low = "antiquewhite", high = "aquamarine4") +
  theme_bw() +
  facet_grid(~name,
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
       height = length(unique(fit_table$enzyme)) * 18 + 600,
       width = 1500,
       units = "px",
       dpi = 150,
       device = "png")
