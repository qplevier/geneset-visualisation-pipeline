# Title     : Weighted Gene Co-expression Network Analysis
# Objective : To identify clusters of highly correlated genes, knows as modules, and relate them to a phenotypic trait
# Created by: Quinten Plevier
# Created on: 22-3-2024

# Load libraries
library(tidyverse)
library(phyloseq)
library(WGCNA)
library(ggnewscale)
library(readxl)
options(stringsAsFactors = FALSE)

# Extract parameters and output directory from snakemake object
micro0 <- snakemake@params[["micro"]]
outdir <- snakemake@output[["outdir"]]

# Create subdirectory for output if it does not exist
if (!file.exists(outdir)) { dir.create(file.path(outdir)) }

# Data input and cleaning

# Load Phyloseq and expression data
physeq <- read_rds(snakemake@input[["physeq"]]) %>%
  subset_samples(micro == micro0)

datExpr <- otu_table(physeq) %>%
  as.data.frame() %>%
  rownames_to_column("enzyme") %>%
  filter(!grepl("\\|", enzyme) &
           enzyme != "UNGROUPED" &
           enzyme != "UNMAPPED") %>%
  mutate(enzyme = str_remove(enzyme, ".*(?=\\[\\d\\.)")) %>%
  distinct() %>%
  group_by(enzyme) %>%
  summarise(across(everything(), sum), .groups = 'drop') %>%
  distinct() %>%
  column_to_rownames("enzyme") %>%
  t() %>%
  as.data.frame()

# Cluster samples to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")

png(
  file.path(outdir, "1_dendrogram_outliers.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "Samples",
  cex.lab = 1.5,
  cex.main = 2
)
abline(h = snakemake@params[["cutoff"]], col = "red") # Add a line to indicate the cut-off for outlier detection
dev.off()

# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = snakemake@params[["cutoff"]], minSize = 10)
datExpr <- datExpr[(clust == 1),] # Clust 1 contains the samples we want to keep.

# Re-cluster samples
sampleTree <- hclust(dist(datExpr), method = "average")

# Filter samples and genes with a lot of zero values
gsg <- goodSamplesGenes(datExpr, verbose = 3)

# Remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))

  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))

  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Load meta data
datTraits <- sample_data(physeq) %>%
  as("data.frame") %>%
  mutate(time = as.numeric(str_remove(time, "T"))) %>%
  .[intersect(rownames(.), rownames(datExpr)),] %>%
  binarizeCategoricalColumns(
    dropFirstLevelVsAll = FALSE,
    includePrefix = FALSE,
    levelSep = "",
    nameForAll = "")

# Convert continuous data traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
png(
  file.path(outdir, "1_dendrogram_samples.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap",
  xlab = "Samples",
  cex.lab = 1.5,
  cex.main = 2
)
dev.off()

# Automatic network construction and module detection

# Construct Network Automatic
enableWGCNAThreads(nThreads = snakemake@threads[[1]])

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
png(
  file.path(outdir, "2_topology_soft_power.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
par(mfrow = c(1, 2))

# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  font = ifelse(powers %in% sft$powerEstimate, 4, 1),
  col = "red"
)

# This line corresponds to using an R^2 cut-off of h
abline(h = 0.85, col = "red")

# Plot the mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  font = ifelse(powers %in% sft$powerEstimate, 4, 1),
  col = "red"
)
dev.off()

# Determine power and create modules

# Create modules
net <- blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(outdir, paste0("wgcna_", micro0)),
  verbose = 3,
  nThreads = snakemake@threads
)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
png(
  file.path(outdir, "3_dendro_samples_modules.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

# Store module Labels, Colors and Module Eigengene
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

save(
  MEs,
  moduleLabels,
  moduleColors,
  geneTree,
  file = file.path(outdir, "networkConstruction-auto.RData")
)

# Relating modules to external information and identifying important genes

# Recalculate MEs with color labels
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes %>%
  orderMEs()

# Compute module-trait correlations and p-values
nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create text matrix for displaying in the plot
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
png(
  file.path(outdir, "4_corr_traits_modules.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)
dev.off()

moduleTraitCorNA <- moduleTraitCor
moduleTraitPvalueNA <- moduleTraitPvalue
is.na(moduleTraitCorNA) <- moduleTraitPvalue > 0.05
is.na(moduleTraitPvalueNA) <- moduleTraitPvalue > 0.05

# Create text matrix for displaying in the plot
textMatrix <- ifelse(!is.na(moduleTraitCorNA), paste0(signif(moduleTraitCorNA, 2), "\n(", signif(moduleTraitPvalueNA, 2), ")"), "")
dim(textMatrix) <- dim(moduleTraitCorNA)

# Display the correlation values within a heatmap plot
png(
  file.path(outdir, "4_corr_traits_modules_signif.png"),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCorNA,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)
dev.off()

# Define variable time containing the time column of datTrait
time <- as.data.frame(datTraits$time) %>%
  rename("time" = "datTraits$time")

# Define variable modNames containing the names (colors) of the modules
modNames <- substring(names(MEs), 3)

# Correlation (membership) between genes and the modules and calculate corresponding p-values
geneModuleMembership <- datExpr %>%
  cor(MEs, use = "p") %>%
  as.data.frame()
names(geneModuleMembership) <- paste0("MM", modNames)

MMPvalue <- geneModuleMembership %>%
  as.matrix() %>%
  corPvalueStudent(nSamples) %>%
  as.data.frame()
names(MMPvalue) <- paste0("p.MM", modNames)

# Correlation (membership) between genes and the trait(s) and calculate corresponding p-values
geneTraitSignificance <- datExpr %>%
  cor(time, use = "p") %>%
  as.data.frame()
names(geneTraitSignificance) <- paste0("GS.", names(time))

GSPvalue <- geneTraitSignificance %>%
  as.matrix() %>%
  corPvalueStudent(nSamples) %>%
  as.data.frame()
names(GSPvalue) <- paste0("p.GS.", names(time))

module <- moduleTraitPvalue %>%
  as.data.frame() %>%
  arrange(time) %>%
  slice(1) %>%
  rownames() %>%
  substring(3)
column <- match(module, modNames)
moduleGenes <- moduleColors == module

# Plot module membership vs. gene significance
png(
  file.path(outdir, "5_module_gene_signif.png"),
  width = 7,
  height = 7,
  units = "in",
  res = 300
)
par(mfrow = c(1, 1))
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for time",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module
)
dev.off()

# Create initial gene information dataframe
geneInfo0 <- data.frame(
  enzymes = names(datExpr),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)

# Order modules by their significance for time
modOrder <- MEs %>%
  cor(time, use = "p") %>%
  abs() %>%
  order(decreasing = TRUE)

# Add module membership information in the chosen order
for (mod in seq_len(ncol(geneModuleMembership)))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <-
    data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
               MMPvalue[, modOrder[mod]])
  names(geneInfo0) <-
    c(oldNames,
      paste0("MM.", modNames[modOrder[mod]]),
      paste0("p.MM.", modNames[modOrder[mod]]))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneInfo <-
  geneInfo0[order(geneInfo0$moduleColor, -abs(geneInfo0$GS.time)),]

# Write geneInfo to CSV file
write.csv(geneInfo, file = sprintf("%s/geneInfo.csv", outdir))
collectGarbage()

# Network visualization using WGCNA functions

# Calculate topological overlap anew.
# Transform TOM with a power to make moderately strong connections more visible in the heatmap
plotTOM <-
  (1 - TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate))^7

# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA

# Call the plot function
# How lighter the olor in the heatmap, the more overlap
png(
  file.path(outdir, "6_topology_overlap_matrix.png"),
  width = 12,
  height = 12,
  units = "in",
  res = 300
)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

# Add the time to existing module eigengenes
MET <- orderMEs(cbind(MEs, time))

# Plot the relationships among the eigengenes and the trait
png(
  file.path(outdir, "7_eigengene_network.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 300
)
par(cex = 0.9)
plotEigengeneNetworks(
  MET,
  "",
  marDendro = c(0, 4, 1, 5),
  marHeatmap = c(3, 4, 1, 2),
  cex.lab = 0.8,
  xLabelsAngle = 90
)
dev.off()

# Exporting a gene network to external visualization software

# Recalculate topological overlap if needed
TOM <- TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)

# Select modules
modules <- moduleTraitPvalue %>%
  as.data.frame() %>%
  arrange(time) %>%
  filter(time < 0.25) %>%
  rownames(.) %>%
  substring(3)

# Link back EC to enzyme names
fullGeneName <- physeq %>%
  otu_table() %>%
  as.data.frame() %>%
  rownames_to_column("enzyme") %>%
  filter(!grepl("\\|", enzyme) &
           enzyme != "UNGROUPED" &
           enzyme != "UNMAPPED") %>%
  select(enzyme) %>%
  mutate(enzyme = str_remove(enzyme, ".*\\((expasy|metacyc)\\) ")) %>%
  separate(enzyme, into = c("gene", "enzyme"), sep = "(?=\\[\\d\\.)", fill = "left") %>%
  mutate(gene = ifelse(is.na(gene), "", gene))

datExpr2 <- datExpr %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("enzyme") %>%
  left_join(fullGeneName, by = "enzyme", multiple = "any") %>%
  unite("enzyme", gene, enzyme, sep = "") %>%
  column_to_rownames("enzyme") %>%
  t() %>%
  as.data.frame()

# Select module probes
inModule <- is.finite(match(moduleColors, modules))

modProbes <- names(datExpr2)[inModule]

modGenes <- modProbes %>%
  str_remove(".*(?=\\[\\d\\.)")

# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste0(outdir, "/", "CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt"),
  nodeFile = paste0(outdir, "/", "CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt"),
  weighted = TRUE,
  threshold = 0.2,
  nodeNames = modProbes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]
)

ecs <- geneInfo %>%
  filter(moduleColor == module) %>%
  select(enzymes) %>%
  mutate(enzymes = str_remove_all(enzymes, "(\\[|\\])")) %>%
  pull(enzymes) %>%
  str_replace_all("\\.", "\\\\.") %>%
  paste0(., "\\D")

for (pval in c("P.Value", "adj.P.Val")) {
  dge <- read_excel(snakemake@input[["diff"]], sheet = "allresults") %>%
    filter(grepl(paste(ecs, collapse = "|"), EC)) %>%
    filter(eval(parse(text = pval)) < 0.05) %>%
    rename("time" = "term") %>%
    mutate(time = str_remove(time, "time"),
           EC = str_remove(EC, "\\((expasy|metacyc)\\) "))

  if (nrow(dge) != 0) {
    averageEx <- dge %>%
      select(time, EC, AveExpr) %>%
      mutate(typename = "AveExp", .before = time) %>%
      select(-time)

    dge <- dge %>%
      select(time, EC, logFC, AveExpr) %>%
      mutate(typename = "logFC", .before = time) %>%
      bind_rows(averageEx) %>%
      mutate_at("time", ~replace_na(., ""))

    # Determine the module to plot

    # Create plot
    p <- dge %>%
      arrange(AveExpr) %>% # arrange to abundance level
      mutate(EC = fct_inorder(factor(EC, ordered = TRUE))) %>%
      ggplot(aes(time, reorder(EC, logFC))) +
      geom_tile(aes(fill = logFC)) +
      labs(fill = "Log Fold Change") +
      scale_fill_gradient2(
        midpoint = 0,
        low = "royalblue",
        mid = "white",
        high = "red"
      ) +
      geom_text(aes(label = sprintf("%0.2f", round(logFC, digits = 2))),
                color = "black",
                size = 2) +
      new_scale("fill") +
      geom_tile(data = dge %>%
        filter(typename == "AveExp"),
                aes(fill = AveExpr)) +
      labs(fill = "Average\nexpression\n(log2)") +
      scale_fill_gradient(low = "antiquewhite", high = "aquamarine4") +
      facet_grid(~typename,
                 scales = "free_x",
                 space = "free_x") +
      labs(title = paste(micro0, "- Module:", module),
           subtitle = "Up and down in the differential genes",
           x = "Time",
           y = "Enzyme",
           caption = paste0(pval, " < 0.05")) +
      theme(axis.text.y = element_text(size = 0.8)) +
      theme_bw() +
      scale_y_discrete(
        label = function(x)
          stringr::str_trunc(x, 60, "center")
      )
    plot(p)

    # Save plot
    ggsave(
      file.path(outdir, paste0("diffs_", micro0, "_", pval, ".png")),
      width = 2000,
      height = 500 + length(unique(dge$EC)) * 40,
      dpi = 200,
      units = "px"
    )
  }
}