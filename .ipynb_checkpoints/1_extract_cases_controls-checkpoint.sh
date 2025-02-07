#!/bin/bash


encounters_file="data/encounters_input.tsv"
bmi_file="data/bmi_input.tsv"
dem_file="data/demographics_input.tsv"
phenotype_file="data/phenotypes.tsv"
descriptive_characteristic_file='data/descriptive_characteristics.tsv'

~/miniconda3/envs/py37/bin/python3.7 utils/cases_controls.py get_cases_controls \
    --encounters_file $encounters_file \
    --bmi_file $bmi_file \
    --dem_file $dem_file \
    --output_file $phenotype_file
    
~/miniconda3/envs/py37/bin/python3.7  utils/demographics.py get_descriptive_characteristics \
    --phenotype_file $phenotype_file \
    --output_file $descriptive_characteristic_file