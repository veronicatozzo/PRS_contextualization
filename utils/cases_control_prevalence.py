import pandas as pd 
import sys
import fire 
import numpy as np
import warnings
from constants import *

def get_first_occurrence(dems, encs, codes, phenotype):
    print(f'Filtering encounters based on IDs..')
    encs = encs[encs.ID.isin(dems.ID)]
    encs.encounter_date = pd.to_datetime(encs.encounter_date)
    
    print(f'Filtering codes based on record type..')
    condition_codes = codes[codes.Record_type == 'Condition'] #Diagnostic Codes for either BC or CHD 
    other_codes = codes[codes.Record_type == 'Other'] #Other is either history codes for BC or procedure codes for CHD 

    def _aggregate(x, condition_codes, other_codes, phenotype):
        other = {}
        diag = {}
        has_other = x[x.code.isin(other_codes.Code)]
        has_diagnostic = x[x.code.isin(condition_codes.Code)]
        other['other_count'] = len(has_other.code)
        other['First_occurrence'] = has_other['encounter_date'].min()
        diag['diag_count'] = len(has_diagnostic.code)
        diag['First_occurrence'] = has_diagnostic['encounter_date'].min()
        if phenotype=='CHD':
            combined = {'other_count': other['other_count'], 'diag_count': diag['diag_count'],
                        'First_occurrence': other['First_occurrence'] if other['other_count'] >= 1 else diag['First_occurrence']}
        elif phenotype=='BC':
            combined = {'other_count': other['other_count'], 'diag_count': diag['diag_count'],
                        'First_occurrence': diag['First_occurrence'] if diag['diag_count'] >= 1 else other['First_occurrence']}
        return pd.Series(combined, index = ['other_count', 'diag_count', 'First_occurrence'])

    encs.sort_values(by='encounter_date', inplace=True)
    
    if phenotype=='BC':
        print(f'Filtering encounters based on codes for {phenotype}')
        encs_out = encs.groupby(['ID'], as_index=False).apply(_aggregate, condition_codes, other_codes, phenotype)
        encs_out = encs_out[(encs_out.other_count>=2)|(encs_out.diag_count>=1)].reset_index()
    elif phenotype=='CHD':
        print(f'Filtering encounters based on codes for {phenotype}')
        encs_out = encs.groupby(['ID'], as_index=False).apply(_aggregate, condition_codes, other_codes, phenotype)
        encs_out = encs_out[(encs_out.other_count>=1)|(encs_out.diag_count>=3)].reset_index()
    return encs_out


def get_controls(dems, encs, codes):
    #Filter out those patients that have relevant diagnosis codes
    print('Filtering controls based on IDs and codes..')
    encs = encs[encs.ID.isin(dems.ID)]
    encs_controls = encs.groupby(['ID']).filter(lambda x: not any(x.code.isin(codes.Code)))
    print(f'Number of controls after ensuring at least 2 encounters {encs_controls.nunique()}')

    print('Making sure there are at least two encounters in the system..')
    encs_controls = encs_controls.groupby(['ID'], as_index=False).apply(lambda group: pd.Series({
        'num_encounters': group['encounter_date'].nunique(),
        'last_recorded_encounter_date': group['encounter_date'].max()}))
    controls = encs_controls[encs_controls.num_encounters>1][['ID', 'last_recorded_encounter_date']]
    print(f'Number of controls after ensuring at least 2 encounters {controls.shape[0]}')
    
    return controls      


def get_bmi(df, df_bmi):
    df.Time = pd.to_datetime(df.Time)
    df_bmi.encounter_date = pd.to_datetime(df_bmi.encounter_date)
    print(df.shape[0])
    bmis = []
    for i in range(df.shape[0]):
        aux = df_bmi[(df_bmi.ID==df.ID.iloc[i])&
                         (df_bmi.encounter_date>=(df.Time.iloc[i]-pd.Timedelta(days=365*2)))&
                          (df_bmi.encounter_date<=(df.Time.iloc[i]+pd.Timedelta(days=60)))]
        bmis.append(np.nanmedian(aux.bmi))
    df['bmi'] = bmis
    return df

def get_phenotype(encs, bmis, dems, phenotype):
    # Load diagnosis codes corresponding to phenotype and filter to female for BC
    if phenotype=='BC':
        full_dems = dems.copy()
        dems = dems[dems.self_identified_sex=='Female']
        codes = pd.read_csv("data_updated/phenotyping_codes_BC_final.csv")
    elif phenotype=='CHD':
        full_dems = dems.copy()
        codes = pd.read_csv("data_updated/phenotyping_codes_CHD_final.csv")
    else:
        print(f'Unknown phenotype {phenotype}. Exiting..')
        sys.exit(0)

    print('\n-----------------------------------------------')
    print(f'Started extraction of {phenotype} cases/controls')
    print('-----------------------------------------------\n')

    # Get first occurrence of diagnosis code in population of interest
    # first_occurrence = get_first_occurrence(dems, encs, codes, phenotype)

    # print(f'Found {first_occurrence.shape[0]} individuals with {phenotype} diagnosis codes')
    # print('Checking to see if they had encounters in two years prior..')
    # # Check that they had two encounters in two years prior to first occurrence
    # cases = check_presence_in_EHR_before_first_occurrence(encs, first_occurrence)
    cases = get_first_occurrence(dems, encs, codes, phenotype)
    print(f'Found N={cases.shape[0]} cases of {phenotype} in data.')

    # Get controls
    print('Getting controls..')
    controls = get_controls(dems, encs, codes)
    print(f'Found N={controls.shape[0]} controls of {phenotype} in data.')

    # Format results 
    cases = cases[['ID', 'First_occurrence']]
    cases.columns = ['ID', 'Time']
    cases[phenotype] = 1
    controls = controls[['ID', 'last_recorded_encounter_date']]
    controls.columns =  ['ID', 'Time']
    controls[phenotype] = 0
    cases_controls = pd.concat([cases, controls], axis=0)
    
    cases_controls = get_bmi(cases_controls, bmis)
    cases_controls.columns = ['ID', f'{phenotype}_date', phenotype,f'{phenotype}_BMI']
    full_dems = full_dems[['ID', 'birth_date']]
    full_dems = full_dems.merge(cases_controls, on='ID', how='left')
    full_dems[f'{phenotype}_date'] = pd.to_datetime(full_dems[f'{phenotype}_date'])
    full_dems.birth_date = pd.to_datetime(full_dems.birth_date)
    full_dems[f'{phenotype}_age'] = (full_dems[f'{phenotype}_date']-full_dems.birth_date).dt.days/365.25
    # if phenotype=='BC':
    #     print(f'\n-----Making sure {phenotype} controls has age ≥18 years')
    #     full_dems = full_dems[(full_dems[phenotype] == 1)|((full_dems[phenotype] == 0) & (full_dems[f'{phenotype}_age'] >= 18))]
        
    print(f'\n-----Finished extraction of {phenotype} cases/controls')
    
    return full_dems[['ID',  phenotype,f'{phenotype}_date',f'{phenotype}_BMI',f'{phenotype}_age']] 

def f_SIRE(x):
    if x.ethnicity_grouped=='Hispanic':
        return 'Hispanic'
    elif str(x.race_grouped)=='nan':
        return 'Other'
    else:
        return x.race_grouped

def get_cases_controls(encounters_file, bmi_file, dem_file, output_file):
    # Load relevant tables
    encs = pd.read_csv(encounters_file, sep='\t')
    bmis = pd.read_csv(bmi_file, sep="\t")
    dems = pd.read_csv(dem_file, sep="\t")

    dems = dems.drop_duplicates(subset=['ID'])
    bmis = bmis.drop_duplicates(subset=['ID', 'encounter_date'])
    encs = encs.drop_duplicates(subset=['ID', 'encounter_date', 'code'])
    print(dems.shape)
    dems = dems[dems.self_identified_sex.isin(['Male', 'Female'])]
    print(dems.shape)

    #Extract BC cases/controls
    BC_cases_controls = get_phenotype(encs, bmis, dems, 'BC')

    #Extract CHD cases/controls
    CHD_cases_controls = get_phenotype(encs, bmis, dems, 'CHD')

    #Merge dataframes 
    dems = dems.merge(BC_cases_controls, on='ID')
    dems = dems.merge(CHD_cases_controls, on='ID')
    dems = dems.drop_duplicates(subset=['ID'])

    dems['race_grouped'] = dems.race.astype(str).map(map_race)
    dems['ethnicity_grouped'] = dems.ethnicity.astype(str).map(map_ethnicity)

    dems['SIRE'] = dems.apply(f_SIRE, axis=1)
    dems

    dems= dems[['ID', 'self_identified_sex', 'SIRE', 'ancestry_pred']+
             [f'PC{i}' for i in range(1, 11)]+
             ['BC', 'BC_date', 'BC_BMI', 'BC_age']+
             ['CHD', 'CHD_date', 'CHD_BMI', 'CHD_age']]
    #Save dataframe to output 
    dems.to_csv(output_file, index=None, sep='\t')

if __name__=='__main__':
    warnings.simplefilter('ignore')
    fire.Fire()