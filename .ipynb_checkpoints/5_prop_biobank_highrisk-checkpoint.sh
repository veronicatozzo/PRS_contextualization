#!/bin/bash

pgs_file="data/calculated_PGS.tsv"
data_file="../data_updated/phenotypes.tsv"
path_to_results='results/'

module load python/3.11
#depending on your system you can either pass ~/miniconda3/envs/py37/bin/python3.7 or just python

python utils/prop_comparison.py compare_proportions \
    --pgs_file $pgs_file \
    --data_file $data_file \
    --path_to_results $path_to_results