#!/bin/bash

PGS_file="data/calculated_PGS.tsv"
phenotype_file="data/phenotypes.tsv"
path_to_results='results/'

~/miniconda3/envs/py37/bin/python3.7 utils/PGS_analysis.py calculate_effect_sizes \
    --pgs_file $PGS_file \
    --phenotype_file $phenotype_file \
    --path_to_results $path_to_results
    
