#!/bin/bash

phenotype_file_path="data/biobank_demo_file.csv"
output_file_path="data/full_PMBB_demo_characteristics.tsv"
biobank_name='biobankname' #i.e PMBB

#depending on your system you can either pass ~/miniconda3/envs/py37/bin/python3.7 or just python
python full_biobank_demographics.py get_descriptive_characteristics \
    --phenotype_file $phenotype_file_path \
    --output_file $output_file_path \
    --biobank $biobank_name