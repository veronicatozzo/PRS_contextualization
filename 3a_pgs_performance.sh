#!/bin/bash

PGS_file="data/calculated_PGS.tsv"
phenotype_file="data/phenotypes.tsv"
output_file='results/PGS_performance.csv'

module load R

Rscript utils/PGS_performance.R \
    --path_to_pgs_file=$PGS_file \
    --path_to_phenotype_file=$phenotype_file \
    --output_file=$output_file