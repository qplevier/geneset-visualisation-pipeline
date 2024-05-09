# Title     : Differential genes in GSEA output
# Objective : Visualise DEGs in the output of GSEA
# Created by: Quinten Plevier
# Created on: 22-4-2024 09:59

### Load libraries
library(tidyverse)
library(reshape2)
library(phyloseq)
library(readxl)
library(ggplot2)

# Set variables
if (!dir.exists(file.path(snakemake@output[["outdir"]]))){
  dir.create(file.path(snakemake@output[["outdir"]]), showWarnings = FALSE)
}
filter <- "adj.P.Val"
filterValueDiffs <- 0.05
filterValueGenesets <- 0.25
plot_on <- "species"

dir_gsea <- "data/gsea_output"

### Load phyloseq object of countdata and metadata
physeq <- read_rds(snakemake@input[["physeq"]])

# counts <- otu_table(physeq) %>%
#   as.data.frame() %>%
#   filter(grepl("\\|", rownames(.))) %>%
#   rownames_to_column("rxn") %>%
#   separate("rxn", c("rxn", "EC"), sep = ": ", fill = "left") %>%
#   separate("EC", c("EC", "species"), sep = "\\|") %>%
#   separate("species", c("genus", "species"), sep = "\\.") %>%
#   mutate(genus = str_remove(genus, "g__"),
#          species = str_remove(species, "s__")) %>%
#   replace_na(list(species = "unclassified"))

# meta <- sample_data(physeq) %>%
#   data.frame() %>%
#   rownames_to_column("variable")

### Load differential gene expression data
for (diff in snakemake@input[["diffs"]]) {
  assign(paste0("dge_", str_remove(diff, ".*_")), read_excel(diff, sheet = "allresults") %>%
    rename("time" = "term") %>%
    mutate(time = str_remove(time, "time"),
           EC = str_remove(EC, "\\((expasy|metacyc)\\) ")) %>%
    separate(EC, c("rxn", "EC"), sep = ": "))
}

### Loops for creating the figures
for (Micro in snakemake@params[["micro"]]) {
  counts <- physeq %>%
    subset_samples(micro == Micro) %>%
    otu_table() %>%
    as.data.frame() %>%
    filter(grepl("\\|", rownames(.))) %>%
    rownames_to_column("rxn") %>%
    separate("rxn", c("rxn", "EC"), sep = ": ", fill = "left") %>%
    separate("EC", c("EC", "species"), sep = "\\|") %>%
    separate("species", c("genus", "species"), sep = "\\.") %>%
    mutate(genus = str_remove(genus, "g__"),
           species = str_remove(species, "s__")) %>%
    replace_na(list(species = "unclassified"))

  meta <- physeq %>%
    subset_samples(micro == Micro) %>%
    sample_data() %>%
    data.frame() %>%
    rownames_to_column("variable")

  for (timepoint in snakemake@params[["time"]]) {
    # Input gene set and positive and negative GSEA output data
    dir <-
      sprintf("%s/%s_%s", dir_gsea, Micro, timepoint)
    inputGenesets <- file.path(dir, "edb/gene_sets.gmt")
    inputNegSets <-
      list.files(
        pattern = "gsea_report_for_na_neg.*\\.tsv",
        recursive = TRUE,
        full.names = TRUE,
        path = dir
      )
    inputPosSets <-
      list.files(
        pattern = "gsea_report_for_na_pos.*\\.tsv",
        recursive = TRUE,
        full.names = TRUE,
        path = dir
      )

    # Load gene sets
    gene_sets <- inputGenesets %>%
      read.table(
        sep = "\t",
        fill = TRUE,
        quote = "\"",
        col.names = paste0("V", seq_len(max(
          count.fields(inputGenesets, sep = "\t"), na.rm = TRUE
        )))
      ) %>%
      select(-V2) %>%
      rename("geneset" = "V1") %>%
      melt("geneset", value.name = "rxn") %>%
      select(-variable) %>%
      filter(rxn != "")

    # Input gene sets
    posSets <- inputPosSets %>%
      read.table(sep = "\t",
                 fill = TRUE,
                 header = TRUE) %>%
      filter(FDR.q.val < filterValueGenesets) %>%
      select(NAME, SIZE, ES, NES, NOM.p.val, FDR.q.val)

    negSets <- inputNegSets %>%
      read.table(sep = "\t",
                 fill = TRUE,
                 header = TRUE) %>%
      filter(FDR.q.val < filterValueGenesets) %>%
      select(NAME, SIZE, ES, NES, NOM.p.val, FDR.q.val)

    signSets <- rbind(posSets, negSets) %>%
      rename("geneset" = "NAME")

    ### Filter differentials on time and join with significant gene sets
    df_data <- get(paste0("dge_", Micro, ".xlsx")) %>%
      inner_join(gene_sets, by = "rxn", relationship = "many-to-many") %>%
      inner_join(signSets, by = "geneset") %>%
      distinct() %>%
      arrange(desc(NES))

    # Figure for the significant differentials in significant gene sets
    df_data2 <- df_data %>%
      filter(eval(parse(text = filter)) < filterValueDiffs) %>%
      mutate(geneset = tolower(geneset))

    # If statement implemented to avoid error when the dataframe is empty
    if (nrow(df_data2) != 0) {
      ggplot(df_data2, aes(time, reorder(EC, logFC), fill = logFC)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%0.2f", round(logFC, digits = 2))),
                  color = "black",
                  size = 2) +
        scale_fill_gradient2(
          midpoint = 0,
          low = "blue",
          mid = "white",
          high = "red"
        ) +
        labs(
          title = sprintf(
            "%s, %s",
            str_to_sentence(Micro),
            timepoint
          ),
          subtitle = "Up and down in the differential genes if available",
          x = "Gene set & Time",
          y = "Enzyme",
          caption = sprintf(
            "%s < %s for differentials, FDR q-value < %s for gene sets",
            filter,
            filterValueDiffs,
            filterValueGenesets
          )
        ) +
        theme(axis.text.y = element_text(size = 0.8),
              strip.text.x = element_text(size = 3)) +
        theme_bw() +
        facet_grid(
          . ~ factor(geneset, levels = unique(geneset)),
          scales = "free_x",
          space = "free_x",
          switch = "x",
          labeller = label_wrap_gen(20)
        )

      # Save figure
      ggsave(
        sprintf(
          "%s/%s_%s_sign_differentials_in_sign_genesets.png",
          snakemake@output[["outdir"]],
          Micro,
          timepoint
        ),
        width = length(unique(df_data2$geneset)) * 200 + 800,
        height = nrow(df_data2) * 15 + 300,
        dpi = 150,
        units = "px", limitsize = FALSE
      )

      # Preprocess data for figure
      top3 <- signSets %>%
        arrange(FDR.q.val) %>%
        .[1:3,]

      for (Geneset in top3$geneset) {
        df_data2 <- df_data %>%
          filter(geneset == Geneset)

        plot_data <- counts %>%
          mutate(EC = str_remove(EC, "\\((expasy|metacyc)\\) ")) %>%
          melt() %>%
          inner_join(df_data2, by = c("rxn", "EC"), relationship = "many-to-many") %>%
          filter(value != 0) %>%
          left_join(meta, by = "variable") %>%
          filter(time.x == timepoint &
                   time.y == timepoint &
                   micro == Micro) %>%
          mutate(id_sample = fct_relevel(id_sample, sort(unique(id_sample))),
                 EC = sprintf("%s\nLogFC: %s\n%s",
                              str_replace(EC, " (?=\\[\\d\\..*\\])", "\n"),
                              round(as.numeric(logFC), 3),
                              ifelse(adj.P.Val < 0.05, "*", "ns"))) %>%
          select(-variable, -time.x, -time.y, -t, -z.std, -micro)

        if (plot_on == "genus") {
          plot_data <- plot_data %>%
            group_by(EC, genus, id_sample, logFC, FDR.q.val, NES) %>%
            summarise(value = sum(value), .groups = "keep") %>%
            rename("species" = "genus")
        }

        # Create figure of enzymes and species
        ggplot(plot_data, aes(reorder(EC, desc(logFC)), species,
                              color = id_sample,
                              size = value)) +
          geom_point(position = position_dodge(width = .9), aes(group = id_sample)) +
          labs(
            title = sprintf(
              "%s, %s, %s",
              str_to_sentence(Micro),
              timepoint,
              tolower(Geneset)
            ),
            subtitle = sprintf(
              "Significance gene set: %s\nNormalized Enrichment Score: %s\nCPM of species and differential enzymes",
              plot_data$FDR.q.val, round(plot_data$NES, 3)
            ),
            caption = "* = adj. p-value < 0.05\nNot significant = ns = adj. p-value > 0.05",
            y = "Species",
            x = "Enzymes",
            fill = "CPM",
            size = "CPM",
            color = "Sample"
          ) +
          scale_y_discrete(limits = rev) +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            panel.spacing = unit(.2, "lines"),
          ) +
          geom_vline(xintercept = 2:length(unique(plot_data$EC)) - 1 + 0.5,
                     colour = "grey50") +
          geom_hline(yintercept = 2:length(unique(plot_data$species)) - 1 + 0.5,
                     colour = "grey50") +
          scale_colour_brewer(palette = "Paired") +
          scale_size(range = c(2, 8))

        # Save figure
        ggsave(
          sprintf(
            "%s/%s_%s_%s_%s_CPM_diffs_in_genesets.png",
            snakemake@output[["outdir"]],
            plot_on,
            Micro,
            timepoint,
            tolower(Geneset)
          ),
          width = length(unique(plot_data$EC)) * 200 + 400,
          height = length(unique(plot_data$species)) * 30 + 600,
          dpi = 150,
          units = "px"
        )
      }
    }
  }
}
