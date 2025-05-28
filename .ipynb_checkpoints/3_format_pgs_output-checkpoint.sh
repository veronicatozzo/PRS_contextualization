#!/bin/bash

PGS_emerge_file="../results/emerge_pgs/PMBB/score/PMBB_pgs.txt.gz"
PGS_catalog_file="../results/pgs_catalog/PMBB/score/PMBB_pgs.txt.gz"
biobankname="PMBB"
path_to_results='data/'

module load python/3.11
#depending on your system you can either pass ~/miniconda3/envs/py37/bin/python3.7 or just python

python utils/format_pgs_file.py get_formatted_pgs_file \
    --pgs_emerge_file $PGS_emerge_file \
    --pgs_catalog_file $PGS_catalog_file \
    --biobankname $biobankname \
    --path_to_results $path_to_results