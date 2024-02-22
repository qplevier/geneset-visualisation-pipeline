# Geneset visualisation pipeline
### Quinten Plevier, Stage TNO
##### 22-08-2024
#
## User input
Make a folder `data/user_input/` where you put your count data with filetype `.tsv` and your metadata with filetype `.txt`. Change the folders in `config/config.yaml` with your respective file names. Further change the config file to your liking. For now, it possesses the Maaslin2 and baseline option.

## Before running
Make sure you open the interface of Cytoscape. Otherwise, this will lead to errors in the rule `rule cytoscape_enrichmentmap`. This is done with:
```shell
Cytoscape/Cytoscape
```

## Running the snakemake
Change the `-c4` with the amount of course you would like to use.
```shell
snakemake --resources cytoscape_instances=1 -c4
```

## Creating Directed Acyclic Graphs
```shell
snakemake --dag | dot -Tsvg > dag.svg

# Use the command below for a simplified DAG
snakemake --rulegraph | dot -Tsvg > rulegraph.svg
```