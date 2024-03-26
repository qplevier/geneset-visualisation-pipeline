# Title     : Weighted Gene Co-expression Network Analysis
# Objective : To identify clusters of highly correlated genes, knows as modules, and relate them to a consensus of the phenotypic traits
# Created by: Quinten Plevier
# Created on: 22-3-2024

# Load libraries
library(tidyverse)
library(phyloseq)
library(WGCNA)
library(ggnewscale)
library(readxl)
options(stringsAsFactors = FALSE)

# Set variables
micros <- snakemake@params[["micro"]]
outdir <- snakemake@output[["outdir"]]

# Create subdirectory for output if it does not exist
if (!file.exists(outdir)) { dir.create(file.path(outdir)) }

# Data input and cleaning

# Load Phyloseq and expression data
physeq <- read_rds(snakemake@input[["physeq"]])

# Create dataframes for consensus
for (micro0 in micros) {
  assign(paste0("datExpr", micro0), physeq %>%
    subset_samples(micro == micro0) %>%
    otu_table() %>%
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
    as.data.frame())
}

# Data cleaning and pre-processing

# We work with two sets:
nSets <- length(grep("datExpr", ls()))
# For easier labeling of plots, create a vector holding descriptive names of the two sets.

setLabels <- str_to_sentence(paste(micros, "group"))
shortLabels <- str_to_sentence(micros)

# Form multi-set expression data:
multiExpr <- vector(mode = "list", length = nSets)

# Create list of all the dataframes
dataList <- lapply(grep("datExpr", ls(), value = TRUE), get)

# Get the common column names across all data frames
commonNames <- ls() %>%
  grep("datExpr", ., value = TRUE) %>%
  lapply(get) %>%
  lapply(colnames) %>%
  reduce(intersect)

# Create a new list 'multiExpr' to store the intersected data from each data frame
multiExpr <- lapply(dataList, function(df) list(data = df[, commonNames]))
exprSize <- checkSets(multiExpr)

# Cluster trees
sampleTrees <- list()
for (set in 1:nSets)
{
  sampleTrees[[set]] <- hclust(dist(multiExpr[[set]]$data), method = "average")
}
sampleTrees <- lapply(multiExpr, function(expr) hclust(dist(expr$data), method = "average"))

# Plot and save cluster tree
png(
  file = file.path(outdir, "1_SampleClustering.png"),
  width = 12,
  height = 6 * nSets,
  units = "in",
  res = 300
)
par(mfrow = c(nSets, 1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets) {
  plot(
    sampleTrees[[set]],
    main = paste("Sample clustering on all genes in", setLabels[set]),
    xlab = "",
    sub = "",
    cex = 0.7
  )
}
dev.off()

# Adjust the cut height for the data sets for the number of samples
cutHeights <-
  c(snakemake@params[["cutoff"]] * (exprSize$nSamples[1] / exprSize$nSamples[2]),
    snakemake@params[["cutoff"]])
# Re-plot the dendrograms including the cut lines
png(
  file = file.path(outdir, "1_SampleClustering_cut.png"),
  width = 12,
  height = 6 * nSets,
  units = "in",
  res = 300
)
par(mfrow = c(nSets, 1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(
    sampleTrees[[set]],
    main = paste("Sample clustering on all genes in", setLabels[set]),
    xlab = "",
    sub = "",
    cex = 0.7
  )
  abline(h = cutHeights[set], col = "red")
}
dev.off()

# Remove outliers above threshold
for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels <- cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set], minSize = 10)
  # Keep the largest one (labeled by the number 1)
  multiExpr[[set]]$data <- multiExpr[[set]]$data[(labels == 1),]
}

# Check the size of the leftover data
exprSize <- checkSets(multiExpr)

# Remove genes and samples if needed
gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets) {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste(
        "In set",
        setLabels[set],
        "removing samples",
        paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")
      ))
    # Remove the offending genes and samples
    multiExpr[[set]]$data <- multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
  }
  # Update exprSize
  exprSize <- checkSets(multiExpr)
}

# Load meta data
datTraits <- sample_data(physeq) %>%
  as("data.frame") %>%
  mutate(time = as.numeric(str_remove(time, "T"))) %>%
  .[intersect(rownames(.), c(rownames(multiExpr[[1]]$data), rownames(multiExpr[[2]]$data))),] %>%
  rownames_to_column("samples") %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = id_sample,
    values_from = value,
    values_fill = list(value = 0)
  ) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = micro,
    values_from = value,
    values_fill = list(value = 0)
  ) %>%
  column_to_rownames("samples")

# Form a multi-set structure that will hold the clinical traits.
Traits <- lapply(1:nSets, function(set) {
  traitRows <- match(rownames(multiExpr[[set]]$data), rownames(datTraits))
  list(data = datTraits[traitRows,])
})

# Define data set dimensions
nGenes <- exprSize$nGenes
nSamples <- exprSize$nSamples

save(multiExpr,
     Traits,
     nGenes,
     nSamples,
     setLabels,
     shortLabels,
     exprSize,
     file = file.path(outdir, "Consensus-dataInput.RData"))

# One-step automatic network construction and module detection

# Choose a set of soft-thresholding powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

# Initialize a list to hold the results of scale-free analysis
powerTables <- vector(mode = "list", length = nSets)

# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector =
    powers, verbose = 2)[[2]])

# Plot the results:
colors <- c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols <- c(2, 5, 6, 7)
colNames <- c(
  "Scale Free Topology Model Fit",
  "Mean connectivity",
  "Median connectivity",
  "Max connectivity"
)

# Get the minima and maxima of the plotted points
ylim <- matrix(NA, nrow = 2, ncol = 4)

for (set in 1:nSets)
{
  for (col in seq_along(plotCols)) {
    ylim[1, col] <- min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] <- max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
png(
  file = file.path(outdir, "2_scaleFreeAnalysis.png"),
  width = 12,
  height = 12,
  units = "in",
  res = 300
)
par(mfcol = c(2, 2))
par(mar = c(4.2, 4.2, 2.2, 0.5))
cex1 <- 0.7

for (col in seq_along(plotCols)) {
  for (set in 1:nSets) {
    if (set == 1) {
      plot(
        powerTables[[set]]$data[, 1],
        -sign(powerTables[[set]]$data[, 3]) * powerTables[[set]]$data[, 2],
        xlab = "Soft Threshold (power)",
        ylab = colNames[col],
        type = "n",
        ylim = ylim[, col],
        main = colNames[col]
      )
      addGrid()
    }
    if (col == 1)
    {
      text(
        powerTables[[set]]$data[, 1],
        -sign(powerTables[[set]]$data[, 3]) * powerTables[[set]]$data[, 2],
        labels = powers,
        cex = cex1,
        col = colors[set]
      )
      abline(h = 0.85, col = "red")
    } else {
      text(
        powerTables[[set]]$data[, 1],
        powerTables[[set]]$data[, plotCols[col]],
        labels = powers,
        cex = cex1,
        col = colors[set]
      ) }
    if (col == 1) {
      legend("bottomright",
             legend = setLabels,
             col = colors,
             pch = 20)
    } else {
      legend("topright",
             legend = setLabels,
             col = colors,
             pch = 20)
    }
  }
}
dev.off()

# Identify the optimal power per dataset
sftPowers <- sapply(powerTables, function(x) which(x$data$SFT.R.sq > 0.85)[1])

net <- blockwiseConsensusModules(
  multiExpr,
  power = sftPowers,
  minModuleSize = 30,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE,
  verbose = 5,
  nThreads = snakemake@threads
)

consMEs <- net$multiMEs
moduleLabels <- net$colors
# Convert the numeric labels to color labels
moduleColors <- labels2colors(moduleLabels)
consTree <- net$dendrograms[[1]]

png(
  file = file.path(outdir, "3_ConsensusDendrogram-auto.png"),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plotDendroAndColors(
  consTree,
  moduleColors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Consensus gene dendrogram and module colors"
)
dev.off()

save(consMEs,
     moduleLabels,
     moduleColors,
     consTree,
     file = file.path(outdir, "NetworkConstruction-auto.RData"))

micro_to_relate <- snakemake@params[["micro_to_relate"]]
datExprRelate <- eval(parse(text = paste0("datExpr", micro_to_relate)))
datExpr <- datExprRelate[, intersect(names(datExprRelate),
                                     names(multiExpr[[which(micros == micro_to_relate)]]$data))]

# # Cluster samples to detect outliers
# sampleTree <- hclust(dist(datExpr), method = "average")
#
# # Determine cluster under the line
# clust <- cutreeStatic(sampleTree, cutHeight = 4400, minSize = 10)
# # Clust 1 contains the samples we want to keep.
# datExpr <- datExpr[(clust==1), ]
# nGenes <- ncol(datExpr)
# nSamples <- nrow(datExpr)
# # Re-cluster samples
# sampleTree <- hclust(dist(datExpr), method = "average")

# Filter samples and genes with a lot of zero values
gsg <- goodSamplesGenes(datExpr, verbose = 3)

# if (!gsg$allOK)
# {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes) > 0)
#     printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
#
#   if (sum(!gsg$goodSamples) > 0)
#     printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
#
#   # Remove the offending genes and samples from the data:
#   datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
# }

# Load meta data
datTraits <- sample_data(physeq) %>%
  as("data.frame") %>%
  mutate(time = as.numeric(str_remove(time, "T"))) %>%
  .[intersect(rownames(.), rownames(datExpr)),]

# Automatic network construction and module detection

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


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
  saveTOMFileBase = file.path(outdir, "wgcna"),
  verbose = 3,
  nThreads = snakemake@threads
)

# Store module Labels, Colors and Module Eigengene
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# Rename variables to avoid conflicts
relateLabels <- moduleLabels
relateColors <- moduleColors
relateTree <- geneTree
relateMEs <- orderMEs(MEs, greyName = "ME0")

lnames <- load(file.path(outdir, "NetworkConstruction-auto.RData"))
lnames

# Isolate the module labels in the order they appear in ordered module eigengenes
relateModuleLabels <- substring(names(relateMEs), 3)
consModuleLabels <- substring(names(consMEs[[1]]$data), 3)

# Convert the numeric module labels to color labels
relateModules <- labels2colors(as.numeric(relateModuleLabels))
consModules <- labels2colors(as.numeric(consModuleLabels))

# Numbers of relating and consensus modules
nRelateMods <- length(relateModules)
nConsMods <- length(consModules)

# Initialize tables of p-values and of the corresponding counts
pTable <- matrix(0, nrow = nRelateMods, ncol = nConsMods)
CountTbl <- matrix(0, nrow = nRelateMods, ncol = nConsMods)

# Execute all pairwaise comparisons
for (fmod in 1:nRelateMods)
  for (cmod in 1:nConsMods)
  {
    relateMembers <- (relateColors == relateModules[fmod])

    consMembers <- (moduleColors == consModules[cmod])

    pTable[fmod, cmod] <- -log10(fisher.test(relateMembers, consMembers, alternative = "greater")$p.value)

    CountTbl[fmod, cmod] <- sum(relateColors == relateModules[fmod] &
                                  moduleColors ==
                                    consModules[cmod])
  }

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] <- 1.3 * max(pTable[is.finite(pTable)])
pTable[pTable > 50] <- 50

# Marginal counts (really module sizes)
relateModTotals <- apply(CountTbl, 1, sum)
consModTotals <- apply(CountTbl, 2, sum)

# Actual plotting
png(
  file = file.path(outdir, "4_ConsensusVsrelateModules.png"),
  width = 10,
  height = 7,
  units = "in",
  res = 300
)
par(mfrow = c(1, 1))
par(cex = 1.0)
par(mar = c(8, 10.4, 2.7, 1) + 0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(
  Matrix = pTable,
  xLabels = paste(" ", consModules),
  yLabels = paste(" ", relateModules),
  colorLabels = TRUE,
  xSymbols = paste("Cons", consModules, ":", consModTotals),
  ySymbols = paste(str_to_sentence(micro_to_relate), relateModules, ":", relateModTotals),
  textMatrix = CountTbl,
  colors = blueWhiteRed(100)[50:100],
  main = paste("Correspondence of",
               micro_to_relate,
               "set-specific and",
               paste(str_to_title(micros), collapse = "-"),
               "consensus modules"),
  cex.text = 1.0,
  cex.lab = 1.0,
  setStdMargins = FALSE,
  legendLabel = expression("-log"[10] * "(fisher.test p-value)")
)

dev.off()

# Relating consensus modules to external microarray sample information

# Set up variables to contain the module-trait correlations
moduleTraitCor <- list()
moduleTraitPvalue <- list()
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] <-
    cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p")
  moduleTraitPvalue[[set]] <-
    corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
}

# Convert numerical labels to colors for labeling of modules in the plot
MEColors <-
  labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
MEColorNames <- paste0("ME", MEColors)

# Plot the module-trait relationship table for set number 1
set <- 1
textMatrix <- paste0(signif(moduleTraitCor[[set]], 2), "\n(", signif(moduleTraitPvalue[[set]], 2), ")")
dim(textMatrix) <- dim(moduleTraitCor[[set]])

# Open a suitably sized window (the user should change the window size if necessary)
png(
  file = file.path(outdir, paste("5_ModuleTraitRelationships-", micro_to_relate, ".png")),
  width = 10,
  height = 7,
  units = "in",
  res = 300
)
par(mar = c(6, 8.8, 3, 2.2))

labeledHeatmap(
  Matrix = moduleTraitCor[[set]],
  xLabels = names(Traits[[set]]$data),
  yLabels = MEColorNames,
  ySymbols = MEColorNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships in", setLabels[set])
)
dev.off()

# Plot the module-trait relationship table for set number 2
set <- 2
textMatrix <- paste0(signif(moduleTraitCor[[set]], 2), "\n(", signif(moduleTraitPvalue[[set]], 2), ")")
dim(textMatrix) <- dim(moduleTraitCor[[set]])

png(
  file = file.path(outdir, "5_ModuleTraitRelationships-parodontitis.png"),
  width = 10,
  height = 7,
  units = "in",
  res = 300
)
par(mar = c(6, 8.8, 3, 2.2))

labeledHeatmap(
  Matrix = moduleTraitCor[[set]],
  xLabels = names(Traits[[set]]$data),
  yLabels = MEColorNames,
  ySymbols = MEColorNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships in", setLabels[set])
)
dev.off()

for (set in 1:nSets) {
  moduleTraitCor[[set]] <- moduleTraitCor[[set]] %>%
    as.data.frame()
  moduleTraitPvalue[[set]] <- moduleTraitPvalue[[set]] %>%
    as.data.frame()
}

# Initialize matrices to hold the consensus correlation and p-value
consensusCor <- matrix(NA,
                       nrow = nrow(moduleTraitCor[[1]]),
                       ncol = ncol(moduleTraitCor[[1]]))
consensusPvalue <- matrix(NA,
                          nrow = nrow(moduleTraitPvalue[[1]]),
                          ncol = ncol(moduleTraitPvalue[[1]]))

# Find consensus negative correlations
negative <- Reduce(`&`, lapply(moduleTraitCor, `<`, 0))
negative[is.na(negative)] <- FALSE
consensusCor[negative] <- pmax(Reduce(pmax, lapply(moduleTraitCor, function(x) x[negative])))
consensusPvalue[negative] <- pmax(Reduce(pmax, lapply(moduleTraitPvalue, function(x) x[negative])))

# Find consensus positive correlations
positive <- Reduce(`&`, lapply(moduleTraitCor, `>`, 0))
positive[is.na(positive)] <- FALSE
consensusCor[positive] <- pmin(Reduce(pmin, lapply(moduleTraitCor, function(x) x[positive])))
consensusPvalue[positive] <- pmax(Reduce(pmax, lapply(moduleTraitPvalue, function(x) x[positive])))

# Prepare text matrix for heatmap
textMatrix <- paste0(signif(consensusCor, 2),
                     "\n(",
                     signif(consensusPvalue, 2),
                     ")")

# Change the dimensions of textMatrix to match moduleTraitCor
dim(textMatrix) <- dim(moduleTraitCor[[set]])
png(
  file = file.path(outdir, "6_ModuleTraitRelationships-consensus.png"),
  width = 10,
  height = 7,
  units = "in",
  res = 300
)
par(mar = c(6, 8.8, 3, 2.2))

# Generate heatmap
labeledHeatmap(
  Matrix = consensusCor,
  xLabels = names(Traits[[set]]$data),
  yLabels = MEColorNames,
  ySymbols = MEColorNames,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste(
    "Consensus module-trait relationships across\n",
    paste(setLabels, collapse = " and ")
  )
)
dev.off()

consMEs.unord <- multiSetMEs(multiExpr,
                             universalColors = moduleLabels,
                             excludeGrey = TRUE)

# Initialize lists to hold correlation results
GS <- list()
kME <- list()

# Calculate correlations for each set
for (set in 1:nSets)
{
  GS[[set]] <- corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data)
  kME[[set]] <- corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
}

# Combine correlation results from different sets
GS.metaZ <- (GS[[1]]$Z + GS[[2]]$Z) / sqrt(2)
kME.metaZ <- (kME[[1]]$Z + kME[[2]]$Z) / sqrt(2)
GS.metaP <- 2 * pnorm(abs(GS.metaZ), lower.tail = FALSE)
kME.metaP <- 2 * pnorm(abs(kME.metaZ), lower.tail = FALSE)

# Prepare matrices for correlation results
GSmat <- rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP)
kMEmat <- rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP)

# Retrieve amount and names of traits
nTraits <- checkSets(Traits)$nGenes
traitNames <- colnames(Traits[[1]]$data)

# Adjust dimensions and column/row names for GS matrix
dim(GSmat) <- c(nGenes, 6 * nTraits)
probes <- names(multiExpr[[1]]$data)
rownames(GSmat) <- probes

colnames(GSmat) <- spaste(
  c(
    "GS.set1.",
    "GS.set2.",
    "p.GS.set1.",
    "p.GS.set2.",
    "Z.GS.meta.",
    "p.GS.meta"
  ),
  rep(traitNames, rep(6, nTraits))
)

# Same code for kME:
MEnames <- colnames(consMEs.unord[[1]]$data)
nMEs <- checkSets(consMEs.unord)$nGenes
dim(kMEmat) <- c(nGenes, 6 * nMEs)
rownames(kMEmat) <- probes

colnames(kMEmat) <- spaste(
  c(
    "kME.set1.",
    "kME.set2.",
    "p.kME.set1.",
    "p.kME.set2.",
    "Z.kME.meta.",
    "p.kME.meta"
  ),
  rep(MEnames, rep(6, nMEs))
)

# Store all usefull data
info <- data.frame(
  Probe = probes,
  ModuleLabel = moduleLabels,
  ModuleColor = labels2colors(moduleLabels),
  GSmat,
  kMEmat
) %>%
  select_if(~!all(is.na(.)))

write.csv(info,
          file = file.path(outdir, "consensusAnalysis-CombinedNetworkResults.csv"),
          row.names = FALSE,
          quote = FALSE)

# Determine the module of interest
moduleColor <- consModules[which(consModuleLabels == moduleTraitCor[[1]] %>%
  select(time) %>%
  rownames(.) %>%
  .[which.min(consensusPvalue)] %>%
  substring(3))]

# Get ECs in this module
ecs <- info %>%
  filter(ModuleColor == moduleColor) %>%
  select(Probe) %>%
  mutate(Probe = str_remove_all(Probe, "(\\[|\\])")) %>%
  pull(Probe) %>%
  str_replace_all("\\.", "\\\\.") %>%
  paste0(., "\\D")

diffs <- snakemake@input[["diffs"]]

for (diff in diffs) {
  micro1 <- diff %>%
    str_remove(".*DGE_") %>%
    str_remove("\\.xlsx") %>%
    str_to_title()

  for (pval in c("P.Value", "adj.P.Val")) {
    dge <- read_excel(diff, sheet = "allresults") %>%
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
        arrange(logFC) %>%
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
        labs(title = paste(micro1, "- Module:", moduleColor),
             subtitle = "Up and down in the differential genes",
             x = "Time",
             y = "Enzyme",
             caption = paste0(pval, " < 0.05")) +
        theme(axis.text.y = element_text(size = 0.8)) +
        theme_bw()
      plot(p)

      # Save plot
      ggsave(
        file.path(outdir, paste0("diffs_", micro1, "_", pval, ".png")),
        width = 2000,
        height = 500 + length(unique(dge$EC)) * 40,
        dpi = 200,
        units = "px"
      )
    }
  }
}