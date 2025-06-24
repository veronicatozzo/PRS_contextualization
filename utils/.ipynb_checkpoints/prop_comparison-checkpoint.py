import pandas as pd
import numpy as np
import warnings
from os.path import join
from statsmodels.stats.proportion import proportions_ztest
import fire

def compare_proportions(pgs_file, data_file, path_to_results):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        pgs = pd.read_csv(pgs_file, sep='\t')
        data = pd.read_csv(data_file, sep='\t')
        data = data.merge(pgs, on='ID', how = 'left')
        
        data = pd.concat((data, pd.get_dummies(data['SIRE'], dtype=int), 
                          pd.get_dummies(data['self_identified_sex'], dtype=int)),axis=1)
    
        bmi_bins = ['(0.0, 18.5]','(18.5, 24.9]','(24.9, 30.0]', '(30.0, 50.0]']
        sire_bins = ['White', 'Black', 'Hispanic', 'Asian', 'Other']
        age_bins = ['(18, 45]', '(45, 65]', '(65, 100]']
        sex_categories = ['Female', 'Male']
        ancestry_categories = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
        
        PGS = ['BC_eMERGE_GIRA', 'CHD_eMERGE_GIRA', 'BC_PGS000507_GIRA', 'CHD_PGS003725_GIRA']
        diff_in_prop = {}
        for p in PGS:
            phe = p.split('_')[0]
            extracted = p.split('_')[1:]
            pgs_type = '_'.join(extracted)
            
            polygenic_threshold = {'BC':5, 'CHD':5}
            q =np.nanquantile(data[p], (100-polygenic_threshold[phe])/100) 
            data[f'High_Risk_{phe}_{pgs_type}'] = (data[f'{phe}_{pgs_type}']>=q).astype(int)
        
            data[f"{phe}_Age_binned"] = pd.cut(data[f"{phe}_age"],[18, 45, 65, 100,110]).astype(str)
            data[f"{phe}_BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50,80]).astype(str)
        
            if phe == 'BC':
                data_pgs = data[data.self_identified_sex == 'Female']
            else:
                data_pgs = data.copy()        
                
            for category in bmi_bins + age_bins + sex_categories + sire_bins + ancestry_categories:
                if category in ancestry_categories:
                    x1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1) & (data_pgs['ancestry'] == category)].shape[0]
                    n1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1)].shape[0]
                    p1 = x1/n1
                    x2 = data_pgs[(data_pgs['ancestry'] == category)].shape[0]
                    n2 = data_pgs.shape[0]
                    p2 = x2/n2
                elif category in bmi_bins:
                    x1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1) & (data_pgs[f"{phe}_BMI_binned"] == category)].shape[0]
                    n1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1)].shape[0]
                    p1 = x1/n1
                    x2 = data_pgs[(data_pgs[f"{phe}_BMI_binned"] == category)].shape[0]
                    n2 = data_pgs.shape[0]
                    p2 = x2/n2
                elif category in age_bins:
                    x1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1) & (data_pgs[f"{phe}_Age_binned"] == category)].shape[0]
                    n1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1)].shape[0]
                    p1 = x1/n1
                    x2 = data_pgs[(data_pgs[f"{phe}_Age_binned"] == category)].shape[0]
                    n2 = data_pgs.shape[0]
                    p2 = x2/n2
                elif category in sex_categories:
                    x1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1) & (data_pgs['self_identified_sex'] == category)].shape[0]
                    n1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1)].shape[0]
                    p1 = x1/n1
                    x2 = data_pgs[(data_pgs['self_identified_sex'] == category)].shape[0]
                    n2 = data_pgs.shape[0]
                    p2 = x2/n2
                elif category in sire_bins:
                    x1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1) & (data_pgs['SIRE'] == category)].shape[0]
                    n1 = data_pgs[(data_pgs[f'High_Risk_{phe}_{pgs_type}']==1)].shape[0]
                    p1 = x1/n1
                    x2 = data_pgs[(data_pgs['SIRE'] == category)].shape[0]
                    n2 = data_pgs.shape[0]
                    p2 = x2/n2
        
                diff_prop = (p1 - p2)
                SE_diff_prop = np.sqrt(((p1*(1-p1))/n1) + ((p2*(1-p2))/n2))
        
                se_highrisk = np.sqrt((p1 * (1 - p1))/n1)
                se_biobank = np.sqrt((p2 * (1 - p2))/n2)
                
                success_cnts = np.array([x1, x2])
                total_cnts = np.array([n1, n2])
                test_stat, pval = proportions_ztest(count=success_cnts, nobs=total_cnts, alternative='two-sided')

                if f"{phe}_{pgs_type}" not in diff_in_prop:
                    diff_in_prop[f"{phe}_{pgs_type}"] = []
                    
                diff_in_prop[f"{phe}_{pgs_type}"].append({
                    'Category': category,
                    'Prop. High Risk': p1,
                    'SE High Risk': se_highrisk,
                    'Prop. Biobank': p2,
                    'SE Biobank': se_biobank,
                    'Diff Prop': diff_prop,
                    'SE Diff Prop': SE_diff_prop,
                    'Zscore': test_stat,
                    'p-val': pval, 
                })
            pd.DataFrame(diff_in_prop[f"{p}"]).to_csv(join(path_to_results, p+'_prop_comparison.tsv'), sep='\t', index=None)

if __name__=='__main__':
    fire.Fire()