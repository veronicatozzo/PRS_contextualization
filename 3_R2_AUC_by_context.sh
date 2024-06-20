#!/bin/bash

PGS_file="data/calculated_PGS.tsv"
phenotype_file="data/phenotypes.tsv"
path_to_results='results/'
output_file='results/R2_AUC_results.csv'

~/miniconda3/envs/py37/bin/python3.7 utils/PGS_analysis.py format_PGS_for_R2_analysis \
    --pgs_file $PGS_file \
    --phenotype_file $phenotype_file \
    --path_to_results $path_to_results
    

Rscript utils/r2_score.R \
    --path_to_pgs_file=$path_to_results\
    --output_file=$output_file