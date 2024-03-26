#!/bin/bash

# Title     : GSEA
# Objective : Determine whether a predefined set of genes shows statistically significant differences between two biological states.
# Created by: Quinten Plevier
# Created on: 27-2-2024

# Execute GSEA
do_not_touch/GSEA_4.3.3/gsea-cli.sh GSEAPreranked -gmx "${snakemake_input[gmt]}" -collapse No_Collapse -norm None -nperm 1000 -rnd_seed 149 -rnk "${snakemake_input[rank]}" -scoring_scheme weighted -rpt_label "${snakemake_params[rpt_label]}" -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -set_max 1000 -set_min 5 -zip_report false -out "${snakemake_params[out_dir_gsea]}"

# Remove last part of directory name (a bunch of random numbers, looks ugly and has no use)
mv "${snakemake_output[output_dir_gsea]}".* "${snakemake_output[output_dir_gsea]}"

# Look for a directory with format "jan01" (month+day) and remove if it exists
sleep 2
if [ -d "$(date +'%b%d' | tr '[:upper:]' '[:lower:]')" ];
then
    # If "rm" throws an error, it is catched by "|| true"
    rm -rf "$(date +'%b%d' | tr '[:upper:]' '[:lower:]')" || true
fi
