Code for the context-based analysis of PRS in breast cancer and coronary heart disease. For information or bug reporting contact sandra.lapinska@pennmedicine.upenn.edu. 

Dependencies:
- python >3.7
    - pandas 
    - numpy
    - fire
- R/4.2 
     - rcompanion
     - pROC
     - optparse
- pgsc_calc (see https://pgsc-calc.readthedocs.io/en/latest/ and pgsc_calc_guide.md for instructions)



## Step 1
The first step of the pipeline requires the identification of cases and controls and associated demographics. It can be obtained by launching the script `./1_extract_cases_control.sh`. It might take a while to extract everything to match codes in the best way possible, so if you want to leave the script and forget about it run everything in the background as `./1_extract_cases_controls.sh > log.txt &` The user will need to specify four paths within the script: 

- the path to a tab separated file containing all the **in-person encounters** with corresponding diagnosis codes of the patients in the biobank. The definition of in-person is left to each biobank. 

- the path to a tab separated file containing all available BMI measurements.

- the path to a tab separated file containing the demographic information of the patients as well as the first 10 PCs obtained projecting onto 1000 Genomes. 

- the path where to save the processed data table

- the path where to save the descriptive characteristics to put at the spreadsheet (https://docs.google.com/spreadsheets/d/1b7pdyeMqVFnSUUYCU4NW5sv78hS2CSA4fELrU5gnnHk/edit?usp=sharing)

The headers of each of the files should be the following 
encounters_file = [‘ID’, ‘encounter_date’, ‘inpatient’, ‘code’]
bmi_file = [‘ID’, ‘encounter_date’, ‘bmi’]
dem_file =  ['ID', 'self_identified_sex', 'birth_date',  'ancestry', 'race', 'ethnicity', 'DeathDate', 'last_recorded_encounter_date', 'First_occurrence', "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"]

All dictionaries for the input and output files are available in `utils/dictionaries.xlsx`. The dictionary for the output file is in the sheet `phenotype_table`.

**If Race and Ethnicity are available at a fine-grain level, they will be processed and grouped together using the maps specified in the file `utils/constants.py`, if there are  subcategories in your biobank that are not present in the map, please add them before running the code or define your own mapping to macro-categories (Asian, Black, Hispanic, White, Other).** 

**This script will run both incidence and prevalence case definitions. If the biobank does not want to calculate incident cases, then just comment out the lines that produce the incidence files** 

## Step 2 
The second step of the pipeline is the computation of the PGS using pgsc_calc. It can be obtained by launching the script `./2_compute_PGS.sh`. As the results will depend on the specification in the samplesheet file, it is important to manually format the scores after this step so that they can be processed in step 3. 
The score will be in  file called biobanknamespecifiedinsamplesheet_pgs.txt.gz that is a long format table with ID x PGS having columns SUM (unnormalized pgs), Z1_norm (pgs adjusted to have mean 0 across ancestries), Z2_norm (pgs adjusted to have mean 0 and unit standard deviation across ancestries). You might have two files for running the scores file from eMERGE and the PGS from the pgs catalog. Please combined them in one unique file with columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'], by the selecting the unnormalized pgs column SUM. For code on how to do this, please look at `pgsc_calc_setup.md` as well as run `./3_format_pgs_output.sh`. This will take in the two scores file from eMERGE and the PGS from the pgs catalog as well as the name of the biobank specific in the samplesheet passed as a string.

Also look into the matchfile in  biobanknamespecifiedinsamplesheet_summary.csv and extract the overlap of variants with the used score. 

## Step 3
The third step of the pipeline will do the following:
`./3a_pgs_performance.sh`: will determine the performance of the PGS prior to doing any analysis to ensure that each biobank's PGS is performing as expected 

Evaluate effect sizes for PGS across and within context by launching `3b_context_effect_size.sh`. In the script the user will need to specify the path to the phenotype file generated at step 1, the PGS file with the columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'] obtained after step 2, and the folder where to store all the results. 

To determine which effect sizes are significantly different from one another when implementing the eMERGE ancestry calibration method then launch `3c_pgs_sig_diff_ORs.sh`. In the script the user will need to specify the path to the phenotype file generated at step 1, the PGS file with the columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'] obtained after step 2, and the folder where to store all the results.

Evaluate R2 score and AUC differences by contexts. This can be done by launching `3d_R2_AUC_by_context.sh` after changing the proper filename in the script. Temporary plot visualizations can be obtained with the notebook `utils/R2_score_plots.ipynb`

## Step 4
The fourth step of the pipeline will compare the performance across contexts within each ancestry group and is located at `4_var_within_anc.sh`. In the script the user will need to specify the path to the phenotype file generated at step 1, the PGS file with the columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'] obtained after step 2, and the folder where to store all the results.

## Step 5
The fifth step of the pipeline will compare the proportion of each context in those that are at high polygenic risk to the proportion in the biobank and is located at `4_var_within_anc.sh`. In the script the user will need to specify the path to the phenotype file generated at step 1, the PGS file with the columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'] obtained after step 2, and the folder where to store all the results.

