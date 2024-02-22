#!/bin/bash

do_not_touch/GSEA_4.3.3/gsea-cli.sh GSEAPreranked -gmx "${snakemake_input[gmt]}" -collapse No_Collapse -norm None -nperm 1000 -rnd_seed 149 -rnk "${snakemake_input[rank]}" -scoring_scheme weighted -rpt_label "${snakemake_params[rpt_label]}" -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 5 -set_max 1000 -set_min 5 -zip_report false -out "${snakemake_params[out_dir_gsea]}"

mv "${snakemake_output[output_dir_gsea]}".* "${snakemake_output[output_dir_gsea]}"

if [ -d "$(date +'%b%d' | tr '[:upper:]' '[:lower:]')" ]; then
    rm -r "$(date +'%b%d' | tr '[:upper:]' '[:lower:]')"
fi
