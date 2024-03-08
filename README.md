# Geneset visualisation pipeline
### Quinten Plevier, Stage TNO
##### 22-08-2024

## User input
Make a folder `data/user_input/` where you put your count data with filetype `.tsv` and your metadata with filetype `.txt`. Change the folders in `config/config.yaml` with your respective file names. Further change the config file to your liking. 

## Activate conda and change working directory
Activate your conda environment with Snakemake and change the working directory to your directory of choice.
```shell
conda activate [Name of your snakemake conda]
cd [Path to working directory]
```

## Before running
Make sure you open the interface of Cytoscape. Otherwise, this will lead to errors in the rule `rule cytoscape_enrichmentmap`. This is done with:
```shell
Cytoscape/Cytoscape
```

## Running the snakemake
Change the `4` with the amount of cores you would like to use.
```shell
snakemake --resources cytoscape_instances=1 --use-conda -c 4
```

## Creating Directed Acyclic Graphs
```shell
snakemake --dag -c 1 | dot -Tsvg > dag.svg

# Use the command below for a simplified DAG
snakemake --rulegraph -c 1 | dot -Tsvg > rulegraph.svg
```

![rulegraph.svg](rulegraph.svg)