import pandas as pd 
import numpy as np
import fire
import warnings

map_race = {'Afghanistan':'White', #'Middle Eastern', 
            'African':'Black',#'Black or African American', 
            'African American':'Black',#'Black or African American',
            'Alaska Native':'Other',#'American or Alaska Native',
            'Alaska Native or American Indian or Indigenous: Not Listed':'Other',#'American or Alaska Native',
            'American Indian':'Other',#'American or Alaska Native', 
            'American Indian or Alaska Native':'Other',#'American or Alaska Native',
            'Arab':'White', #'Middle Eastern',
            'Armenian':'White', #'Middle Eastern', 
            'Asian':'Asian', 
            'Asian: Asian Indian':'Asian', #'South Asian',
            'Asian: Chinese':'Asian', #'East Asian', 
            'Asian: Filipino':'Asian', #'South Asian',
            'Asian: Indonesian':'Asian', #'South Asian',
            'Asian: Japanese':'Asian', #'East Asian', 
            'Asian: Korean':'Asian', #'East Asian', 
            'Asian: Not Listed':'Asian',
            'Asian: Other':'Asian',
            'Asian: Pakistani':'Asian', #'South Asian',
            'Asian: Taiwanese':'Asian', #'East Asian', 
            'Asian: Thai':'Asian', #'South Asian',
            'Asian: Vietnamese':'Asian', #'South Asian',
            'Assyrian':'White', #'Middle Eastern', 
            'Bangladeshi':'Asian', #'South Asian',
            'Black':'Black',#'Black or African American', 
            'Black or African American: Not Listed':'Black',#'Black or African American', 
            'Black: African American':'Black',#'Black or African American', 
            'Burmese':'Asian', #'South Asian', 
            'Cambodian':'Asian', #'South Asian', 
            'Caribbean/West Indian':'Other',#'Two or more',
            'Chinese':'Asian', #'East Asian', 
            'Choose Not to Answer':'Other',#'None of the above',
            'Declined to Specify':'Other',#'None of the above',
            'Do not identify with Race':'Other',#'None of the above',
            'Egyptian':'White', #''Middle Eastern', 
            'European':'White',
            'Fijian':'Asian', #'South Asian', 
            'Filipino':'Asian', #'South Asian', 
            'Guamanian/Chamorro':'Other',#'Native Hawaiian or Pacific Islander',
            'Haitian':'Black',#'Black or African American',
            'Hmong':'Asian',#'East Asian', 
            'Indian (India)':'Asian',#'South Asian',
            'Indigenous':'Other',#'American or Alaska Native',
            'Indonesian':'Asian',#'South Asian',
            'Iranian':'White', #'Middle Eastern', 
            'Iraqi':'White', #'Middle Eastern', 
            'Israeli':'White', #'Middle Eastern', 
            'Japanese':'Asian',#'East Asian', 
            'Korean':'Asian',#'East Asian',  
            'Lebanese':'White', #'Middle Eastern', 
            'Middle Eastern or North African: Not Listed':'White',# 'Middle Eastern', 
            'Native Hawaiian':'Other',#'Native Hawaiian or Pacific Islander',
            'Native Hawaiian or Other Pacific Islander':'Other',#'Native Hawaiian or Pacific Islander',
            'Not Listed':'Other',#'None of the above',
            'Other Race':'Other',#'None of the above',
            'Pacific Islander: Guamanian or Chamorro':'Other',#'Native Hawaiian or Pacific Islander',
            'Pacific Islander: Native Hawaiian':'Other',#'Native Hawaiian or Pacific Islander',
            'Pacific Islander: Not Listed':'Other',#'Native Hawaiian or Pacific Islander',
            'Pacific Islander: Other':'Other',#'Native Hawaiian or Pacific Islander',
            'Pacific Islander: Samoan':'Other',#'Native Hawaiian or Pacific Islander', 
            'Pakistani':'Asian',#''South Asian',
            'Palestinian':'White', #'Middle Eastern', 
            'Samoan/American Samoan':'Other',#'Native Hawaiian or Pacific Islander', 
            'Sri Lankan':'Asian',#'South Asian',
            'Syrian':'White', #'Middle Eastern', 
            'Taiwanese':'Asian',#''East Asian',  
            'Thai':'Asian',#'South Asian',
            'Unknown':'Other',#'None of the above',
            'Vietnamese':'Asian',#''South Asian',
            'White':"White", 
            'White or Caucasian':"White", 
            'White: Not Listed':"White",  
            'nan':'Other',
            'Black or African American':'Black',
            'Some Other Race':'Other',
            'American Indian or Alaskan Native':'Other',
            'HLW-Hispanic Latino/White':'Other',
            'Patient Declined':'Other',
            'East Indian':'Asian', 
            'HLB-Hispanic Latino/Black':'Other'}#"None of the above"}

map_ethnicity= {'Choose Not to Answer':'Unknown',
                'Cuban':'Hispanic',
                'Hispanic or Latino':'Hispanic',
                'Hispanic/Spanish origin Other':'Hispanic',
                'Mexican, Mexican American, Chicano/a':'Hispanic',
                'Not Hispanic or Latino':'Non-Hispanic',
                'Puerto Rican':'Hispanic',
                'Unknown':'Unknown',
                'nan':'Unknown',
               'Hispanic Latino':'Hispanic',
                'WHI':'Non-Hispanic',
                'Patient Declined':'Unknown',
                'WHITE':'Non-Hispanic',
                'Blk':'Non-Hispanic',
                'deactivate':'Unknown',
                'HLW':'Non-Hispanic',
                'Other':'Unknown',
                'Black':'Non-Hispanic'}

def f_SIRE(x):
    if x.ethnicity_grouped=='Hispanic':
        return 'Hispanic'
    elif str(x.race_grouped)=='nan':
        return 'Other'
    else:
        return x.race_grouped

def get_descriptive_characteristics(phenotype_file, output_file, biobank):
    full = pd.read_csv(phenotype_file, sep='\t')
    full['ancestry'] = full['ancestry'].apply(lambda x: x if x in ['EUR', 'AFR', 'AMR', 'EAS', 'SAS'] else 'Unclassified')

    if 'age' not in full.columns:
        full['birth_date'] = pd.to_datetime(full['birth_date'])
        today = pd.Timestamp('today')
        full['age'] = today.year - full['birth_date'].dt.year

    if 'SIRE' not in full.columns:
        full['race_grouped'] = full.race.astype(str).map(map_race)
        full['ethnicity_grouped'] = full.ethnicity.astype(str).map(map_ethnicity)
        full['SIRE'] = full.apply(f_SIRE, axis=1)

    age_brakets = [(0, 18), (18, 45), (45, 65), (66, 100)]
    bmi_brakets = [(0, 18.5), (18.5, 25), (25, 30), (30, 70)]
    res = [f"{full.shape[0]}"]
    # Age
    res.append(f"{np.round(np.mean(full.age), 1)} ({np.round(np.std(full.age), 1)})")
    for ages in age_brakets:
        aux = full[(full.age>ages[0])&(full.age<=ages[1])]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/full.shape[0]*100, 1)})")
    # BMI  
    res.append(f"{np.round(np.mean(full.BMI), 0)} ({np.round(np.std(full.BMI), 1)})")
    for bmis in bmi_brakets:
        aux = full[(full.BMI>bmis[0])&(full.BMI<=bmis[1])]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/full.shape[0]*100, 1)})")
    # Sex
    for sex in ['Male', 'Female']:
        aux = full[(full.self_identified_sex==sex)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/full.shape[0]*100, 1)})")
    # Ancestry
    for ancestry in ['EUR', 'AFR', 'AMR', 'EAS', 'SAS','Unclassified']:
        aux = full[(full.ancestry==ancestry)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/full.shape[0]*100, 1)})")
    # SIRE
    for sire in ['White', 'Black', 'Hispanic', 'Asian', 'Other']:
        aux = full[(full.SIRE==sire)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/full.shape[0]*100, 1)})")

    descr = pd.DataFrame([res], 
                         columns=['N', 
                                  'Age, mean(SD)', '<=18', '19-45', '46-65', '>=66', 
                                  'BMI, mean(SD)', '<18.5', '18.5-24.9', '25-29.9', '>=30', 
                                  'Male', 'Female', 
                                  'EUR', 'AFR', 'AMR', 'EAS', 'SAS','Unclassified',
                                  'White', 'Black', 'Hispanic', 'Asian', 'Other'],
                         index=[biobank]).T
    descr.to_csv(output_file, sep='\t')

if __name__=='__main__':
    fire.Fire()