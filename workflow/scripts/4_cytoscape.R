library(RCy3)
library(tidyverse, quietly = TRUE)
library(glue)
library(tools)

cytoscapePing()
# closeSession(FALSE)

rootFolder <- file_path_as_absolute(snakemake@input[["gsea"]])
gmtFile <- file_path_as_absolute(snakemake@input[["gmt"]])
expressionFile <- file_path_as_absolute(snakemake@input[["expression"]])

#################################
# Create enrichment map
#################################
glue(
  "enrichmentmap mastermap \\
   rootFolder={rootFolder} \\
   coefficients=COMBINED \\
   pvalue=1.0 \\
   qvalue=0.25 \\
   similaritycutoff=0.375 \\
   combinedConstant=0.5 \\
   filterByExpressions=false \\
   commonGMTFile={gmtFile}
   commonExpressionFile={expressionFile}"
) %>% commandsGET
Sys.sleep(10)

#################################
# Filter out clusters with very few nodes
#################################
# importFilters(filename = filterFile)
# applyFilter(filter.name = 'Node filter')
# selectEdgesAdjacentToSelectedNodes()
# deleteSelectedEdges()
# deleteSelectedNodes()

#################################
# Annotate the clusters
#################################
if (typeof(getAllEdges()) != "NULL") {
  commandsGET('autoannotate annotate-clusterBoosted labelColumn="EnrichmentMap::GS_DESCR" minWordOccurrence=1')
  Sys.sleep(1)
  commandsGET("layout autoannotate-cose-cluster springStrength=40 repulsionStrength=30")# useCatchallCluster=true")#
  # This takes time, but the API returns immediately. Sleep for 5 sec and hope
  # that the annotation is complete in this time
  Sys.sleep(3)
} else {
  commandsGET("layout grid nodeHorizontalSpacing=200 nodeVerticalSpacing=50")
  Sys.sleep(3)
}

# nodes <- getAllNodes()
# setNodeLabelBypass(node.names = nodes, new.labels = '')

Sys.sleep(3)
fitContent()

Sys.sleep(2)

# system(glue('mkdir -p {dirname(snakemake@output$image)}'))
exportImage(snakemake@output[["image"]], type = "SVG")
