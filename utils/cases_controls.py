import pandas as pd 
import sys
import fire 
import numpy as np
import warnings
from constants import *



def get_first_occurrence(dems, encs, codes):
    # filter encounters that have codes
    print(f'Filtering encounters based on IDs and diagnosis code..')
    encs = encs[encs.ID.isin(dems.ID)]
    encs = encs[encs.code.apply(lambda x: any(str(x).startswith(str(s)) for s in codes.Code))]
    encs.encounter_date = pd.to_datetime(encs.encounter_date)

    print('Processing outpatient encounters..')
    def _aggregate(x):
        d = {}
        distances =x.encounter_date.diff()
        d['at_least_one_30_days_distance'] = np.any(np.abs(distances.dt.days)>=30)
        d['Count'] = len(x)
        d['encounter_date'] = x['encounter_date'].min()
        return pd.Series(d, index=['at_least_one_30_days_distance', 'Count', 'encounter_date'])

    encs_out = encs[encs.inpatient==0].drop_duplicates(subset=['ID', 'encounter_date'])
    encs_out.sort_values(by='encounter_date', inplace=True)
    encs_out = encs_out.groupby(['ID'], as_index=False).apply(_aggregate)
    encs_out = encs_out[(encs_out.Count>=2)&(encs_out.at_least_one_30_days_distance)]
    encs_out = encs_out[['ID', 'encounter_date']]
    
    print('Processing inpatient encounters..')
    encs_in = encs[encs.inpatient==1]
    encs_in = encs_in.groupby(['ID']).agg({'encounter_date':'min'}).reset_index()
    data = pd.concat([encs_out, encs_in])
    data = data.groupby(['ID'])['encounter_date'].agg(First_occurrence=min).reset_index()
    return data

def check_presence_in_EHR_before_first_occurrence(encs, data):
    encs.encounter_date = pd.to_datetime(encs.encounter_date)
    encs = encs[encs.ID.isin(data.ID)]
    to_keep = []
    for i in range(data.shape[0]):
        aux = encs[encs.ID==data.ID.iloc[i]]
        aux = aux[(encs.encounter_date<=(data.First_occurrence.iloc[i]-pd.Timedelta(days=30)))&
                 (encs.encounter_date>=(data.First_occurrence.iloc[i]-pd.Timedelta(days=365*2)))]
        aux = aux.drop_duplicates(subset=['encounter_date'])
        if aux.shape[0]>=2:
            to_keep.append(data.ID.iloc[i])

    diagnosis = data[data.ID.isin(to_keep)]
    return diagnosis

def get_controls(dems, encs, codes):
    
    #Filter out those patients that have relevant diagnosis codes
    print('Filtering controls based on IDs and codes..')
    encs = encs[encs.ID.isin(dems.ID)]
    encs_controls = encs[~encs.code.apply(lambda x: any(str(x).startswith(str(s)) for s in codes.Code))]
    controls = dems[dems.ID.isin(encs_controls.ID)]
    controls.last_recorded_encounter_date = pd.to_datetime(controls.last_recorded_encounter_date)
    encs.encounter_date = pd.to_datetime(encs.encounter_date)
    
    # Check that there are at least 2 encounters in the two years prior the last recorded encounter
    to_keep=[]
    for i in range(controls.shape[0]):
        aux = encs[encs.ID==controls.ID.iloc[i]]
        aux = aux[(encs.encounter_date<=(controls.last_recorded_encounter_date.iloc[i]-pd.Timedelta(days=30)))&
                 (encs.encounter_date>=(controls.last_recorded_encounter_date.iloc[i]-pd.Timedelta(days=365*2)))]
        aux = aux.drop_duplicates(subset=['encounter_date'])
        if aux.shape[0]>=2:
            to_keep.append(controls.ID.iloc[i])
    controls = controls[controls.ID.isin(to_keep)]
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
        codes = pd.read_csv("data/phenotyping_codes_BC.csv")
    elif phenotype=='CHD':
        full_dems = dems.copy()
        codes = pd.read_csv("data/phenotyping_codes_CHD.csv")
    else:
        print(f'Unknown phenotype {phenotype}. Exiting..')
        sys.exit(0)

    print('\n-----------------------------------------------')
    print(f'Started extraction of {phenotype} cases/controls')
    print('-----------------------------------------------\n')

    # Get first occurrence of diagnosis code in population of interest
    first_occurrence = get_first_occurrence(dems, encs, codes)

    print(f'Found {first_occurrence.shape[0]} individuals with {phenotype} diagnosis codes')
    print('Checking to see if they had encounters in two years prior..')
    # Check that they had two encounters in two years prior to first occurrence
    cases = check_presence_in_EHR_before_first_occurrence(encs, first_occurrence)
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
                
    #Force the ids to be int to handle merging
    encs.ID = encs.ID.astype(int)
    bmis.ID = bmis.ID.astype(int)
    dems.ID = dems.ID.astype(int)

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

    dems= dems[['ID', 'self_identified_sex', 'last_recorded_encounter_date', 'SIRE', 'ancestry']+
             [f'PC{i}' for i in range(1, 11)]+
             ['BC', 'BC_date', 'BC_BMI', 'BC_age']+
             ['CHD', 'CHD_date', 'CHD_BMI', 'CHD_age']]
    #Save dataframe to output 
    dems.to_csv(output_file, index=None, sep='\t')
    
    
if __name__=='__main__':
    warnings.simplefilter('ignore')
    fire.Fire()
