Code for the context-based analysis of PRS in breast cancer and coronary heart disease. For information or bug reporting contact vtozzo@mednet.ucla.edu. 

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
The first step of the pipeline requires the identification of cases and controls and associated demographics. It can be obtained by launching the script `./1_extract_cases_controls.sh`. It might take a while to extract everything to match codes in the best way possible, so if you want to leave the script and forget about it run everything in the backgroun as `./1_extract_cases_controls.sh > log.txt &` The user will need to specify four paths within the script: 

- the path to a tab separated file containing all the **in-person encounters** with corresponding diagnosis codes of the patients in the biobank. The definition of in-person is left to each biobank. 

- the path to a tab separated file containing all available BMI measurements.

- the path to a tab separated file containing the demographic information of the patients as well as the first 10 PCs obtained projecting onto 1000 Genomes. 

- the path where to save the processed data table

- the path where to save the descriptive characteristics to put at the spreadsheet (https://docs.google.com/spreadsheets/d/1b7pdyeMqVFnSUUYCU4NW5sv78hS2CSA4fELrU5gnnHk/edit?usp=sharing) 

All dictionaries for the input and output files are available in `utils/dictionaries.xlsx`. The dictionary for the output file is in the sheet `phenotype_table`.

**If Race and Ethnicity are available at a fine-grain level, they will be processed and grouped together using the maps specified in the file `utils/constants.py`, if there are  subcategories in your biobank that are not present in the map, please add them before running the code or define your own mapping to macro-categories (Asian, Black, Hispanic, White, Other).** 



## Step 2 
The second step of the pipeline is the computation of the PGS using pgsc_calc. As the results will depend on the specification in the samplesheet file, it is important to manually format the scores after this step so that they can be processed in step 3. 
The score will be in  file called biobanknamespecifiedinsamplesheet_pgs.txt.gz that is a long format table with ID x PGS having columns SUM (unnormalized pgs), Z1_norm (pgs adjusted to have mean 0 across ancestries), Z2_norm (pgs adjusted to have mean 0 and unit standard deviation across ancestries). You might have two files for running the scores file from eMERGE and the PGS from the pgs catalog. Please combined them in one unique file with columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'], by the selecting the unnormalized pgs column SUM. 

Also look into the matchfile in  biobanknamespecifiedinsamplesheet_summary.csv and extract the overlap of variants with the used score. 


## Step 3 
Evaluate effect sizes for PGS across and within context by launching `3_context_effect_size.sh`. In the script the user will need to specify the path to the phenotype file generated at step 1, the PGS file with the columns ['ID', 'BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725'] obtained after step 2, and the folder where to store all the results. 


Evaluate R2 score and AUC differences by contexts. This can be done by launching `3_R2_AUC_by_context.sh` after changing the proper filename in the script. Temporary plot visualizations can be obtained with the notebook `utils/R2_score_plots.ipynb`

