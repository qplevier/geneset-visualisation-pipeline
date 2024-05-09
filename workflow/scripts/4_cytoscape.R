# Title     : Cytoscape
# Objective : Visualise the GSEA output with cytoscape and save images of the visualisation
# Created by: Quinten Plevier
# Created on: 27-2-2024

# Load libraries
library(RCy3)
library(tidyverse, quietly = TRUE)
library(glue)
library(tools)

# Connect to a open Cytoscape GUI
cytoscapePing()
# closeSession(save.before.closing = FALSE)

# Set folder and files paths
rootFolder <- file_path_as_absolute(snakemake@input[["gsea"]])
gmtFile <- file_path_as_absolute(snakemake@input[["gmt"]])
expressionFile <- file_path_as_absolute(snakemake@input[["expression"]])
classFile <- file_path_as_absolute(snakemake@input[["classes"]])

# Create enrichment map
glue(
  "enrichmentmap mastermap \\
   rootFolder={rootFolder} \\
   coefficients=COMBINED \\
   pvalue=1.0 \\
   qvalue=0.25 \\
   similaritycutoff=0.375 \\
   combinedConstant=0.5 \\
   filterByExpressions=false \\
   commonGMTFile={gmtFile} \\
   commonExpressionFile={expressionFile} \\
   commonClassFile={classFile}"
) %>% commandsGET
Sys.sleep(5)

# setNodeSizeMapping('EnrichmentMap::gs_size', sizes = c(20, 50))
setNodeFontSizeDefault(8, style.name = getCurrentStyle())
Sys.sleep(1)

# Annotate the clusters with AutoAnnotate
# Sys.sleep is used to give the program some time to calculate everything, because the API returns immediately
tryCatch(
{
  if (length(getAllEdges()) != "NULL") {
    commandsGET('autoannotate annotate-clusterBoosted labelColumn="EnrichmentMap::GS_DESCR" edgeWeightColumn="EnrichmentMap::overlap_size"')
    Sys.sleep(3)

    # Create layout of the nodes
    commandsGET("layout autoannotate-cose-cluster incremental=true springStrength=100 repulsionStrength=10 gravityStrength=50 compoundGravityStrength=50 gravityRange=50 compoundGravityRange=50 smartRepulsionRangeCalc=true smartEdgeLengthCalc=true useCatchallCluster=false")
    # commandsGET("autoannotate layout")
    Sys.sleep(3)
  } else {
    # If there are no edges, set nodes in a grid layout
    commandsGET("layout grid nodeHorizontalSpacing=200 nodeVerticalSpacing=50")
    Sys.sleep(3)
  } },
  error = function(e) {
    commandsGET("layout grid nodeHorizontalSpacing=200 nodeVerticalSpacing=50")
    Sys.sleep(3)
  }
)

# Fit the content to the screen
fitContent()
Sys.sleep(2)

# Create image and export
exportImage(snakemake@output[["image"]], type = "SVG")
