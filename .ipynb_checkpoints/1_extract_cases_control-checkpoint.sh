#!/bin/bash

encounters_file="data/encounters_input.tsv"
bmi_file="data/bmi_input.tsv"
dem_file="data/demographics_input.tsv"
phenotype_inc_file="data/phenotypes_incidence.tsv"
phenotype_prev_file="data/phenotypes_incidence.tsv"
descriptive_characteristic_inc_file='data/descriptive_characteristics_incidence.tsv'
descriptive_characteristic_prev_file='data/descriptive_characteristics_prevalence.tsv'

# Get Incidence 
~/miniconda3/envs/py37/bin/python3.7 utils/cases_control_incidence.py get_cases_controls \
    --encounters_file $encounters_file \
    --bmi_file $bmi_file \
    --dem_file $dem_file \
    --output_file $phenotype_inc_file

~/miniconda3/envs/py37/bin/python3.7  utils/demographics.py get_descriptive_characteristics \
    --phenotype_file $phenotype_inc_file \
    --output_file $descriptive_characteristic_inc_file

# Get Prevalence
~/miniconda3/envs/py37/bin/python3.7 utils/cases_control_prevalence.py get_cases_controls \
    --encounters_file $encounters_file \
    --bmi_file $bmi_file \
    --dem_file $dem_file \
    --output_file $phenotype_prev_file

~/miniconda3/envs/py37/bin/python3.7  utils/demographics.py get_descriptive_characteristics \
    --phenotype_file $phenotype_prev_file \
    --output_file $descriptive_characteristic_prev_file