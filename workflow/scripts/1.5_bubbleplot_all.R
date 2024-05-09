# Title     : Heatmap of top N species over all data and differences
# Objective : To visualise the abundance of species in the samples over all micro and timepoint
# Created by: Quinten Plevier
# Created on: 9-4-2024 13:36

# Load libraries
library(tidyverse, quietly = TRUE)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)

# Set seed for reproductive color creation
set.seed(123)

# Create table from phyloseq object
table <- read_rds(snakemake@input[["physeq"]]) %>%
  psmelt() %>%
  mutate(id_sample = as.factor(id_sample)) %>%
  filter(Abundance > 0) %>%
  separate(OTU, c("OTU", "Species"), "\\|") %>%
  na.omit() %>%
  mutate(Species = str_remove(Species, ".*s__")) %>%
  group_by(Species, id_sample, micro, time) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  group_by(id_sample, micro, time, .drop = FALSE) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  group_by(Species, .drop = FALSE) %>%
  mutate(Sum_abundance = sum(Abundance)) %>%
  ungroup() %>%
  filter(Sum_abundance %in% c(sort(unique(.$Sum_abundance), decreasing=TRUE)[0:snakemake@params[["top"]]]))

healthy <- table %>%
  filter(micro == "healthy") %>%
  select(id_sample) %>%
  unique() %>%
  nrow()
paro <- table %>%
  filter(micro == "parodontitis") %>%
  select(id_sample) %>%
  unique() %>%
  nrow()

diff_species <- table %>%
  group_by(Species, micro) %>%
  summarise(observations = n(),
            ave = n() / ifelse(micro=="healthy", healthy, paro)) %>%
  distinct() %>%
  group_by(Species) %>%
  mutate(diff = abs(lag(ave) / ave)) %>%
  group_by(Species) %>%
  fill(diff, .direction = "up") %>%
  mutate(diff = ifelse(diff < 1, 1/diff, diff)) %>%
  filter(diff >= 2) %>%
  distinct()

table2 <- table %>%
  filter(.$Species %in% diff_species$Species) %>%
  filter(Sum_abundance %in% c(sort(unique(.$Sum_abundance), decreasing=TRUE)[0:snakemake@params[["top"]]]))

meta <- read_rds(snakemake@input[["physeq"]]) %>%
  sample_data() %>%
  as("data.frame")

table2 <- merge(table2, meta, by=c("id_sample", "micro", "time"), all = TRUE) %>%
  fill(Species, .direction = "down")

# Create first heatmap
ggplot(table,
       aes(x = id_sample,
           y = reorder(Species, desc(Species)),
           fill = Relative_abundance)) +
  geom_tile(colour = NA, aes(width=1)) +
  facet_grid(~ micro + time,
             scales = "free",
             space = "free"
             ) +
  labs(x = "",
       y = "",
       size = "Relative\nabundance(%)",
       fill = "Relative\nabundance(%)") +
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 8,
      face = "plain",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      colour = "black",
      face = "plain",
      size = 8
    ),
    legend.text = element_text(
      size = 10,
      face = "plain",
      colour = "black"
    ),
    legend.title = element_text(size = 12, face = "plain"),
    panel.background = element_rect(fill="white", colour="white", linewidth = 0),
    panel.border = element_rect(
      colour = "grey",
      fill = NA,
      linewidth = 1.2
    )) +
  scale_fill_gradient(low="lightblue",
                      high="red",
                      limits = c(0, max(table$Relative_abundance)),
                      na.value = "white") +
  labs(x = "Individual",
       y = "Species",
       title = glue::glue("Top {snakemake@params[['top']]} most abundant species over all samples")) +
  theme(strip.text.x = element_text(size = 12))

# Save the plot
ggsave(snakemake@output[["png"]],
       height = (length(unique(table$Species)) * 40 + 450),
       width = length(unique(meta$id_sample)) * length(unique(meta$micro)) * length(unique(meta$time)) * 20 + 1500,
       units = "px",
       device = "png")

# Create second heatmap with differences
ggplot(table2,
       aes(x = id_sample,
           y = reorder(Species, desc(Species)),
           fill = Relative_abundance)) +
  geom_tile(colour = NA, aes(width=1)) +
  facet_grid(~ micro + time,
             scales = "free",
             space = "free"
             ) +
  labs(x = "",
       y = "",
       size = "Relative\nabundance(%)",
       fill = "Relative\nabundance(%)") +
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 8,
      face = "plain",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      colour = "black",
      face = "plain",
      size = 8
    ),
    legend.text = element_text(
      size = 10,
      face = "plain",
      colour = "black"
    ),
    legend.title = element_text(size = 12, face = "plain"),
    panel.background = element_rect(fill="white", colour="white", linewidth = 0),
    panel.border = element_rect(
      colour = "grey",
      fill = NA,
      linewidth = 1.2
    )) +
  scale_fill_gradient(low="lightblue",
                      high="red",
                      limits = c(0, max(table2$Relative_abundance)),
                      na.value = "white") +
  labs(x = "Individual",
       y = "Species",
       title = "Species with a >2 times difference of occurence between the groups",
       subtitle = "Difference over the total occurcence of species per group divided by amount of samples") +
  theme(strip.text.x = element_text(size = 12))

# Save the plot
ggsave(snakemake@output[["png2"]],
       height = (length(unique(table2$Species)) * 40 + 600),
       width = length(unique(meta$id_sample)) * length(unique(meta$micro)) * length(unique(meta$time)) * 20 + 1500,
       units = "px",
       device = "png")