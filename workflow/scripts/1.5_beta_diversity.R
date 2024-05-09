# Title     : Beta diversity
# Objective : Visualise the beta diversity in the samples through PCA
# Created by: Quinten Plevier
# Created on: 26-2-2024

# Load packages
library(compositions)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggnewscale)

# Load data
physeq <- read_rds(snakemake@input[["physeq"]])

# Get meta data
sampleData <- physeq %>%
  sample_data() %>%
  as("data.frame") %>%
  rownames_to_column("id") %>%
  mutate(id_sample = as.factor(id_sample))

# Create PCA data with prcomp() and centered log ratio
pcaFit <- physeq %>%
  otu_table() %>%
  clr() %>%
  scale() %>%
  prcomp()

# Get the points of the PCAs
pcaPoints <- pcaFit %>%
  unclass() %>%
  enframe() %>%
  filter(name == "rotation") %>%
  pluck("value", 1) %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  inner_join(sampleData)

# Variance explained by principal component
varExpl <- pcaFit %>%
  summary() %>%
  unclass() %>%
  enframe() %>%
  filter(name == "importance") %>%
  pluck("value", 1) %>%
  as.data.frame() %>%
  rownames_to_column("type") %>%
  pivot_longer(-type) %>%
  filter(type == "Proportion of Variance",
         name %in% c("PC1", "PC2")) %>%
  mutate(label = paste(name, scales::percent(value), sep = " - ")) %>%
  pluck("label")

# Calculate centroids of the PCA points at the four times and two micros
pcaCentroids <- pcaPoints %>%
  group_by(time, micro) %>%
  summarise_at(vars(PC1, PC2),
               mean)

# Calculate the path which the centroids take
pcaPath <- pcaCentroids %>%
  filter(time != "T0") %>%
  rename(PCA1_end = PC1,
         PCA2_end = PC2) %>%
  bind_cols(pcaCentroids %>%
              filter(time != "T30") %>%
              ungroup() %>%
              select(PCA1_start = PC1,
                     PCA2_start = PC2))

# Prepare data for the plot
plotPoints1 <- pcaPoints %>%
  mutate(newTime = time) %>%
  mutate(newTime = gsub("T", "", newTime)) %>%
  mutate(newTime = as.numeric(newTime)) %>%
  filter(micro == snakemake@params[["micros"]][1])

plotCentroids1 <- pcaCentroids %>%
  mutate(newTime = time) %>%
  mutate(newTime = gsub("T", "", newTime)) %>%
  mutate(newTime = as.numeric(newTime)) %>%
  filter(micro == snakemake@params[["micros"]][1])

plotPoints2 <- pcaPoints %>%
  mutate(newTime = time) %>%
  mutate(newTime = gsub("T", "", newTime)) %>%
  mutate(newTime = as.numeric(newTime)) %>%
  filter(micro == snakemake@params[["micros"]][2])

plotCentroids2 <- pcaCentroids %>%
  mutate(newTime = time) %>%
  mutate(newTime = gsub("T", "", newTime)) %>%
  mutate(newTime = as.numeric(newTime)) %>%
  filter(micro == snakemake@params[["micros"]][2])

# Plot the PCA
ggplot() +
  geom_point(
    data = plotPoints1,
    aes(PC1,
        PC2,
        fill = newTime),
    colour = "blue",
    shape = 21,
    size = 3
  ) +
  geom_point(
    data = plotCentroids1,
    aes(PC1,
        PC2,
        fill = newTime),
    colour = "blue",
    shape = 21,
    size = 7
  ) +
  labs(fill = sprintf("Time\n%s", snakemake@params[["micros"]][1])) +
  scale_fill_gradient(low = "white", high = "blue") +
  guides(fill = guide_legend(ncol = 1)) +
  new_scale("fill") +
  geom_point(
    data = plotPoints2,
    aes(PC1,
        PC2,
        fill = newTime),
    colour = "red",
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = plotCentroids2,
    aes(PC1,
        PC2,
        fill = newTime),
    colour = "red",
    shape = 24,
    size = 7
  ) +
  labs(fill = sprintf("Time\n%s", snakemake@params[["micros"]][2])) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_segment(
    data = pcaPath,
    aes(
      PCA1_start,
      PCA2_start,
      xend = PCA1_end,
      yend = PCA2_end,
      lty = "black"
    ),
    arrow = arrow(length = unit(.25, "cm")),
    show.legend = F
  ) +
  labs(title = "Principal Component Analysis on the functional data",
       subtitle = "With 95% confidence ellipse",
       x = varExpl[1],
       y = varExpl[2]) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  stat_ellipse(data = plotPoints1,
               mapping = aes(PC1, PC2,
                             group = snakemake@params[["micros"]][1]), color = "blue") +
  stat_ellipse(data = plotPoints2,
               mapping = aes(PC1, PC2,
                             group = snakemake@params[["micros"]][2]), color = "red")

# Save the PCA plot
ggsave(snakemake@output[["png"]],
       width = 3000,
       height = 2000,
       units = "px",
       device = "png")


## TOO LITTLE RAM FOR MY MACHINE (16GB)

# # plot distances
# d <- dist(pcaFit$x)
# df <- melt(as.matrix(d), varnames = c("row", "col"))
#
# distance_pc1 <- dist(pcaFit$x[,1:2])
# distance_pc2 <- dist(pcaFit$x[,2:3])
#
# dfpc1 <- melt(as.matrix(distance_pc1), varnames = c("row", "col"))
# dfpc2 <- melt(as.matrix(distance_pc2), varnames = c("row", "col"))
#
# df <- dfpc1
#
# distancePlot <- df %>%
#   mutate(TrueFalse = case_when(
#     # str_detect(row, "Inulin_1_") & str_detect(col, "Inulin_1_T0") ~ "True",
#     # str_detect(row, "Inulin_2_") & str_detect(col, "Inulin_2_T0") ~ "True",
#     # str_detect(row, "Inulin_3_") & str_detect(col, "Inulin_3_T0") ~ "True",
#     # str_detect(row, "Inulin_4_") & str_detect(col, "Inulin_4_T0") ~ "True",
#     # str_detect(row, "Inulin_5_") & str_detect(col, "Inulin_5_T0") ~ "True",
#     # str_detect(row, "Inulin_6_") & str_detect(col, "Inulin_6_T0") ~ "True",
#     # str_detect(row, "Inulin_7_") & str_detect(col, "Inulin_7_T0") ~ "True",
#     # str_detect(row, "Inulin_8_") & str_detect(col, "Inulin_8_T0") ~ "True",
#     # str_detect(row, "Inulin_9_") & str_detect(col, "Inulin_9_T0") ~ "True",
#     # str_detect(row, "Inulin_10_") & str_detect(col, "Inulin_10_T0") ~ "True",
#      str_detect(row, "Palatinose_11_") & str_detect(col, "Palatinose_11_T0") ~ "True",
#     str_detect(row, "Palatinose_12_") & str_detect(col, "Palatinose_12_T0") ~ "True",
#     str_detect(row, "Palatinose_13_") & str_detect(col, "Palatinose_13_T0") ~ "True",
#     str_detect(row, "Palatinose_14_") & str_detect(col, "Palatinose_14_T0") ~ "True",
#     str_detect(row, "Palatinose_15_") & str_detect(col, "Palatinose_15_T0") ~ "True",
#     str_detect(row, "Palatinose_16_") & str_detect(col, "Palatinose_16_T0") ~ "True",
#     str_detect(row, "Palatinose_17_") & str_detect(col, "Palatinose_17_T0") ~ "True",
#     str_detect(row, "Palatinose_18_") & str_detect(col, "Palatinose_18_T0") ~ "True",
#     str_detect(row, "Palatinose_19_") & str_detect(col, "Palatinose_19_T0") ~ "True",
#     str_detect(row, "Palatinose_20_") & str_detect(col, "Palatinose_20_T0") ~ "True",
#   .default = as.character("False")
# )) %>%
#   filter(TrueFalse == "True") %>%
#   separate(row, c("str","person", "timepoint"), sep="_")
#
# distancePlot %>%
#   ggplot(aes(timepoint, value,group = person, color=person))+
#   # geom_point(size=2)+
#   geom_point(aes(fill = person), size = 3, shape = 21) +
#   geom_line()+
#   theme_bw()+
#   labs(title = "Distances PCA",
#        subtitle = "Isomaltulose - periodontitis individuals",
#        x="Exposure\ntime (min)",
#        y="Euclidean distance")+
#   guides(color=guide_legend(ncol =2)) +
#   labs(fill="Individual") + guides(color="none")
#   theme_bw()
#
# # hierarchical clustering
#
# # Aitchison distance most appropiate for CPM normalized data -> not existing anymore. CLR -> euclidean as alternative
#
# countData = physeq %>%
#   otu_table %>%
#   t() %>% # samples as columns and components as rows
#   as.data.frame() %>%
#   rownames_to_column("names") %>%
#   filter(names != "UNMAPPED") %>%
#   filter(names != "UNGROUPED") %>%
#   column_to_rownames("names") %>%
#   clr %>%
#   t() %>%
#   as.data.frame
#
# dist_mat <- dist(countData, method = 'euclidean')
# hclust_avg <- hclust(dist_mat, method = 'average')
# plot(hclust_avg)





