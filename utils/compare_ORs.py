import numpy as np
import pandas as pd
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt 
from matplotlib import rcParams 
from os.path import join
import warnings
from scipy.stats import norm
from itertools import combinations

def get_glm_coef_full_standardized(data, phe, p):
    results = smf.glm(
                formula=f'{phe} ~ {p}',
                data=data,
                family=sm.families.Binomial(),
                missing='drop'
            ).fit()
    betas = results.params[p]
    se = results.bse[p]
    ORs = np.exp(results.params[p])
    conf = results.conf_int().loc[p]
    lower_ci = np.exp(conf[0])
    upper_ci = np.exp(conf[1])
    pval = results.pvalues[p]
    return [betas, se, ORs, lower_ci, upper_ci, pval]

def get_significance_of_ORs(data, p, path_to_results):
    groups = {
        'Ancestry': ['All', 'EUR', 'AFR', 'AMR', 'EAS', 'SAS'],
        'SIRE': ['All', 'White', 'Black', 'Asian', 'Hispanic', 'Other'],
        'BMI': ['All', '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]'],
        'Age': ['All', '(18, 45]', '(45, 65]', '(65, 100]'],
        'Sex': ['All','Female', 'Male']
    }
    
    results = []
    
    total_comparisons = 0
    for group_list in groups.values():
        num_comparisons = len(list(combinations(group_list, 2)))
        total_comparisons += num_comparisons
        
    for group_title, group_list in groups.items():
        # Get unique pairwise combinations of groups within the group_list
        pairwise_combinations = combinations(group_list, 2)
        for base_group, compare_group in pairwise_combinations:
            try:
                beta1 = data.loc[data.data_partition == base_group, 'betas'].iloc[0]
                SE1 = data.loc[data.data_partition == base_group, 'se'].iloc[0]
    
                beta2 = data.loc[data.data_partition == compare_group, 'betas'].iloc[0]
                SE2 = data.loc[data.data_partition == compare_group, 'se'].iloc[0]
    
                # Calculate the Z-score and P-value
                delta = abs(beta1 - beta2)
                SE_delta = np.sqrt(SE1**2 + SE2**2)
                z = delta / SE_delta
                p_value = 2 * (1 - norm.cdf(z))
    
                nominal_significance = 'Yes' if p_value < 0.05 else 'No'
                # Bonferroni correction
                significance_threshold = 0.05 / total_comparisons
                bonferroni_significance = 'Yes' if p_value < significance_threshold else 'No'
    
                results.append({
                    'Group Category': group_title,
                    'Group 1': base_group,
                    'Group 2': compare_group,
                    'Z-score': z,
                    'P-value': p_value,
                    f'Nominal Significance: 0.05': nominal_significance,
                    f'Bonferroni Significance: {significance_threshold:.2e}': bonferroni_significance
                })
            except IndexError:
                print(f"Missing data for group: {base_group} or {compare_group} in {group_title}")
                continue
    results_df = pd.DataFrame(results)
    results_df.to_csv(join(path_to_results, f"{p}_signif_OR_differences.tsv"), index=None)

def get_gira_results(pgs_file, data_file, path_to_results):
    warnings.simplefilter('ignore')
    pgs = pd.read_csv(pgs_file, sep='\t')
    data = pd.read_csv(data_file, sep='\t')
    data = data.merge(pgs, on='ID')
    data = data[data.self_identified_sex.isin(['Male', 'Female'])]
    data = pd.concat((data,
                      pd.get_dummies(data['SIRE'], dtype=int),
                      pd.get_dummies(data['self_identified_sex'], dtype=int)), axis=1
                    )
    
    final_results = {}
    PGS = ['BC_eMERGE_GIRA', 'CHD_eMERGE_GIRA', 'BC_PGS000507_GIRA', 'CHD_PGS003725_GIRA']
    for p in PGS:
        phe = p.split('_')[0]
    
        data_to_analyze = data[pd.notnull(data[phe])]
        data_to_analyze[f"Age_binned"] = pd.cut(data[f"{phe}_age"], [18, 45, 65, 100]).astype(str)
        data_to_analyze[f"BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50]).astype(str)
        data_to_analyze[p] = data_to_analyze[p]
            
        res = get_glm_coef_full_standardized(data_to_analyze, phe, p) 
        
        res_temp = [['All']+res]
        for context in ['ancestry', 'SIRE', 'BMI_binned', 'Age_binned', 'self_identified_sex']:
            for cont in np.unique(data_to_analyze[context]):
                if str(cont)=='nan' or str(cont)=='Unknown':
                    continue
                try: #If there are no cases in the context the model fails
                    data_to_analyze_cont = data_to_analyze[data_to_analyze[context]==cont]
                    res = get_glm_coef_full_standardized(data_to_analyze_cont, phe, p) 
                    res_temp.append([cont]+res)
                except:
                    continue
        result_df = pd.DataFrame(res_temp, columns=['data_partition', 'betas', 'se', 'OR', 'lower_ci', 'upper_ci', 'pval'])
        result_df.to_csv(join(path_to_results, p+'_OR_results.tsv'), sep='\t', index=None)
        get_significance_of_ORs(result_df, p, path_to_results)