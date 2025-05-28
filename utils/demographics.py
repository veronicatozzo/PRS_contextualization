import pandas as pd
import numpy as np
import fire
import warnings

def demographics(data, full):
    age_brakets = [(0, 18), (18, 45), (45, 65), (66, 100)]
    bmi_brakets = [(0, 18.5), (18.5, 25), (25, 30), (30, 70)]
    res = [f"{data.shape[0]} ({np.round(data.shape[0]/full.shape[0],2)}) ({full.shape[0]})"]
    res.append(f"{np.round(np.mean(data.age), 0)} ({np.round(np.std(data.age), 0)})")
    for ages in age_brakets:
        aux = data[(data.age>ages[0])&(data.age<=ages[1])]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/data.shape[0]*100, 0)})")
    res.append(f"{np.round(np.mean(data.BMI), 0)} ({np.round(np.std(data.BMI), 0)})")
    
    for bmis in bmi_brakets:
        aux = data[(data.BMI>bmis[0])&(data.BMI<=bmis[1])]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/data.shape[0]*100, 0)})")
    for sex in ['Male', 'Female']:
        aux = data[(data.self_identified_sex==sex)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/data.shape[0]*100, 0)})")
    for ancestry in ['EUR', 'AFR', 'AMR', 'EAS', 'SAS','Unclassified']:
        aux = data[(data.ancestry==ancestry)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/data.shape[0]*100, 0)})")
    for sire in ['White', 'Black', 'Hispanic', 'Asian', 'Other']:
        aux = data[(data.SIRE==sire)]
        res.append(f"{aux.shape[0]} ({np.round(aux.shape[0]/data.shape[0]*100, 0)})")
    return res 
    
def get_descriptive_characteristics(phenotype_file, output_file):
    data = pd.read_csv(phenotype_file, sep='\t')
    data['ancestry'] = data['ancestry'].apply(lambda x: x if x in ['EUR', 'AFR', 'AMR', 'EAS', 'SAS'] else 'Unclassified')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aux = data[data.BC==1]
        aux['age'] = aux.BC_age
        aux['BMI'] = aux.BC_BMI
        bc_cases = demographics(aux, data[data.self_identified_sex == 'Female'])

        aux = data[data.BC==0]
        aux['age'] = aux.BC_age
        aux['BMI'] = aux.BC_BMI
        bc_controls = demographics(aux, data[data.self_identified_sex == 'Female'])

        aux = data[data.CHD==1]
        aux['age'] = aux.CHD_age
        aux['BMI'] = aux.CHD_BMI
        chd_cases = demographics(aux, data)

        aux = data[data.CHD==0]
        aux['age'] = aux.CHD_age
        aux['BMI'] = aux.CHD_BMI
        chd_controls = demographics(aux, data)

    
    descr = pd.DataFrame([bc_cases, bc_controls, chd_cases, chd_controls], 
                         columns=['N %', 'Age', '<=18', '19-45', '46-65', '>=66', 
                                'BMI', '<18.5', '18.5-24.9', '25-29.9', '>=30', 
                                'Male', 'Female', 'EUR', 'AFR', 'AMR', 'EAS', 'SAS','Unclassified',
                                'White', 'Black', 'Hispanic', 'Asian', 'Other'],
                         index=['BC_cases', 'BC_controls', 'CHD_cases', 'CHD_controls']).T
    descr.to_csv(output_file, sep='\t')

if __name__=='__main__':
    fire.Fire()