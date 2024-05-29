import pandas as pd 
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt 
from matplotlib import rcParams 
from os.path import join

import warnings
import fire

warnings.simplefilter('ignore')
def get_glm_coef_full_standardized(data, phe, p):
    mdl = sm.GLM(data[phe], 
                 data[p+'_standardized'], 
                 family=sm.families.Binomial(),
                 missing='drop')
    results = mdl.fit()
    res1 = results.summary2().tables[1]
    mdl = sm.GLM(data[phe], 
                 data[p+'_standardized_resid'], 
                 family=sm.families.Binomial(),
                 missing='drop')
    results = mdl.fit()
    res2 = results.summary2().tables[1]
    return [ res1['Coef.'].values[0], res1['Std.Err.'].values[0], res1['P>|z|'].values[0], 
            res2['Coef.'].values[0], res2['Std.Err.'].values[0], res2['P>|z|'].values[0]]

def get_glm_coef_context_standardized(data, phe, p, covs):
    data[p+'_standardized_cont'] = (data[p]-np.mean(data[p]))/np.std(data[p])
    data[p+'_standardized_resid_cont'] = sm.OLS(data[p], 
                                                              sm.add_constant(data[covs]),
                                                              missing='drop').fit().resid
    data[p+'_standardized_resid_cont'] = (data[p+'_standardized_resid_cont']-data[p+'_standardized_resid_cont'].mean())/data[p+'_standardized_resid_cont'].std()
    
    mdl = sm.GLM(data[phe], 
                 data[p+'_standardized_cont'], 
                 family=sm.families.Binomial(),
                 missing='drop')
    results = mdl.fit()
    res1 = results.summary2().tables[1]
    
    mdl = sm.GLM(data[phe], 
                 data[p+'_standardized_resid_cont'], 
                 family=sm.families.Binomial(),
                 missing='drop')
    results = mdl.fit()
    res2 = results.summary2().tables[1]
    return [ res1['Coef.'].values[0], res1['Std.Err.'].values[0], res1['P>|z|'].values[0], 
            res2['Coef.'].values[0], res2['Std.Err.'].values[0], res2['P>|z|'].values[0]]


def plot_effect_sizes(results, file_to_save):
    rcParams.update({'font.size':12})

    fig, ax = plt.subplots(1, 4, figsize=(15, 5))#, sharex=True)
    
    groups = {'Ancestry':['All','AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
             'SIRE':['All','Asian', 'Black', 'Hispanic', 'Other', 'White'],
             'BMI':['All','(0.0, 18.5]','(18.5, 24.9]','(24.9, 30.0]', '(30.0, 50.0]'],
             'Age/Sex':['All','(18, 45]','(45, 65]','(65, 100]','Female','Male']}
    
    for i, (title, l) in enumerate(groups.items()):
        res = results[results.data_partition.isin(l)]
        ax[i].errorbar(res['coef_sd'], 
                     np.arange(0, res.shape[0]*5, 5)[::-1],
                     xerr=res['sd_sd'],
                     ls='none',
                    marker='o', markersize=5, label='Standardized across all') 
        ax[i].axvline(0, linestyle='--', color='k')
        
        for y in np.arange(4, res.shape[0]*5, 5):
            ax[i].axhline(y, linestyle='--', color='lightgray')
        
        ax[i].set_yticks(np.arange(1.5, res.shape[0]*5, 5))
        ax[i].set_yticklabels(res['data_partition'][::-1])

        ax[i].errorbar(res['coef_resid_sd'], 
                     np.arange(1, res.shape[0]*5, 5)[::-1],
                     xerr=res['sd_resid_sd'],
                     ls='none',
                    marker='o', markersize=5, color='C1', label='Contexts removed and standardized across all')  

        ax[i].errorbar(res['coef_cont_sd'], 
                     np.arange(2, res.shape[0]*5, 5)[::-1],
                     xerr=res['sd_cont_sd'],
                     ls='none',
                    marker='o', markersize=5, color='C2', label='Standardized within context') 

        ax[i].errorbar(res['coef_cont_resid_sd'], 
                     np.arange(3, res.shape[0]*5, 5)[::-1],
                     xerr=res['sd_resid_cont_sd'],
                     ls='none',
                    marker='o', markersize=5, color='C3', label='Contexts removed and standardized within context')   
        ax[i].set_title(title)
        ax[i].set_xlabel('Effect size per 1 PGS SD')
    plt.tight_layout()

    ax[0].legend(bbox_to_anchor=(4, 1.3), ncol=2)
    plt.savefig(file_to_save, bbox_inches='tight', transparent=True, dpi=300)
    plt.show()
    
    
def calculate_effect_sizes(pgs_file, data_file, path_to_results):
    PGS = ['BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725']
    pgs = pd.read_csv(pgs_file, sep='\t')
    data = pd.read_csv(data_file, sep='\t')
    data = data.merge(pgs, on='ID')
#     data = data[data.self_identified_sex.isin(['Male', 'Female'])]
    data = pd.concat((data,
                     pd.get_dummies(data['SIRE']),
                    pd.get_dummies(data['self_identified_sex'])),
    axis=1)
    covs = ['Asian', 'Black', 'Hispanic', 'Other',
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8',
           'PC9', 'PC10', 'Female'] 

    final_results = {}
    for p in PGS:
        phe = p.split('_')[0]

        data_to_analyze = data[pd.notnull(data[phe])]
        data_to_analyze[f"Age_binned"] = pd.cut(data[f"{phe}_age"], [18, 45, 65, 100]).astype(str)
        data_to_analyze[f"BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50]).astype(str)
        covs_aux = covs+[f'{phe}_age', f'{phe}_BMI']

        data_to_analyze[p+'_standardized'] = (data_to_analyze[p]-np.mean(data_to_analyze[p]))/np.std(data_to_analyze[p])
        data_to_analyze[p+'_standardized_resid'] = sm.OLS(data_to_analyze[p], 
                                                          sm.add_constant(data_to_analyze[covs_aux]),
                                                          missing='drop').fit().resid
        data_to_analyze[p+'_standardized_resid'] = (data_to_analyze[p+'_standardized_resid']-data_to_analyze[p+'_standardized_resid'].mean())/data_to_analyze[p+'_standardized_resid'].std()
        res = get_glm_coef_full_standardized(data_to_analyze, phe, p) 
        res_temp = [['All']+res+[np.nan]*6]

        for context in ['ancestry', 'SIRE', 'BMI_binned', 'Age_binned', 'self_identified_sex']:
            for cont in np.unique(data_to_analyze[context]):
                if str(cont)=='nan' or str(cont)=='Unknown':
                    continue
                try: #If there are no cases in the context the model fails
                    data_to_analyze_cont = data_to_analyze[data_to_analyze[context]==cont]
                    res = get_glm_coef_full_standardized(data_to_analyze_cont, phe, p) 
                    res_context = get_glm_coef_context_standardized(data_to_analyze_cont, phe, p, covs_aux)
                    res_temp.append([cont]+res+res_context)
                except:
                    continue
        final_results[p] = pd.DataFrame(res_temp, columns=['data_partition', 
                                                           'coef_sd', 'sd_sd', 'p_sd',
                                                          'coef_resid_sd', 'sd_resid_sd', 'p_resid_sd',
                                                          'coef_cont_sd', 'sd_cont_sd', 'p_cont_sd',
                                                          'coef_cont_resid_sd', 'sd_resid_cont_sd', 'p_resid_cont_sd'])
        
        plot_effect_sizes(final_results[p], join(path_to_results, p+'.pdf'))
        final_results[p].to_csv(join(path_to_results, p+'.tsv'), sep='\t', index=None)  
        
        
        
if __name__=='__main__':
    fire.Fire()