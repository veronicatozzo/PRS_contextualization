import pandas as pd 
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt 
from matplotlib import rcParams 
from os.path import join

import warnings
import fire

warnings.simplefilter('ignore')
def get_glm_coef_emerge(data, phe, p):
    mean_pgs = np.mean(data[p+'_gira'])
    results = smf.glm(formula=f"{phe} ~ {p+'_gira'}",
                      data=data,
                      family=sm.families.Binomial(),
                      missing='drop').fit()
    conf = results.conf_int().loc[p+'_gira']
    return [mean_pgs, np.exp(results.params[p+'_gira']), np.exp(conf[0]), np.exp(conf[1]), results.pvalues[p+'_gira']]
    
def get_glm_coef_full_standardized(data, phe, p): 
    mean_pgs_stand = np.mean(data[p+'_standardized'])
    mdl = smf.glm(formula=f"{phe} ~ {p+'_standardized'}",
                  data=data,
                  family=sm.families.Binomial(),
                  missing='drop')
    res1 = mdl.fit()
    conf_res1 = res1.conf_int().loc[p+'_standardized']

    mean_pgs_stand_resid = np.mean(data[p+'_standardized_resid'])
    mdl = smf.glm(formula=f"{phe} ~ {p+'_standardized_resid'}",
                  data=data,
                  family=sm.families.Binomial(),
                  missing='drop')
    res2 = mdl.fit()
    conf_res2 = res2.conf_int().loc[p+'_standardized_resid']
    
    return [mean_pgs_stand, np.exp(res1.params[p+'_standardized']), np.exp(conf_res1[0]), np.exp(conf_res1[1]), res1.pvalues[p+'_standardized'],
            mean_pgs_stand_resid, np.exp(res2.params[p+'_standardized_resid']), np.exp(conf_res2[0]), np.exp(conf_res2[1]), res2.pvalues[p+'_standardized_resid']]

    
def get_glm_coef_context_standardized(data, phe, p, covs):
    data[p+'_standardized_cont'] = (data[p]-np.mean(data[p]))/np.std(data[p])
    mean_pgs_stand_cont = np.mean(data[p+'_standardized_cont'])
    data[p+'_standardized_resid_cont'] = sm.OLS(data[p], 
                                                sm.add_constant(data[covs]),
                                                missing='drop').fit().resid
    data[p+'_standardized_resid_cont'] = (data[p+'_standardized_resid_cont']-data[p+'_standardized_resid_cont'].mean())/data[p+'_standardized_resid_cont'].std()
    mean_pgs_stand_resid_cont = np.mean(data[p+'_standardized_resid_cont'])

    mdl = smf.glm(formula=f"{phe} ~ {p+'_standardized_cont'}",
                  data=data,
                  family=sm.families.Binomial(),
                  missing='drop')
    res1 = mdl.fit()
    conf_res1 = res1.conf_int().loc[p+'_standardized_cont']

    mdl = smf.glm(formula=f"{phe} ~ {p+'_standardized_resid_cont'}",
                  data=data,
                  family=sm.families.Binomial(),
                  missing='drop')
    res2 = mdl.fit()
    conf_res2 = res2.conf_int().loc[p+'_standardized_resid_cont']
    
    return [mean_pgs_stand_cont, np.exp(res1.params[p+'_standardized_cont']), np.exp(conf_res1[0]), np.exp(conf_res1[1]), res1.pvalues[p+'_standardized_cont'], 
            mean_pgs_stand_resid_cont, np.exp(res2.params[p+'_standardized_resid_cont']), np.exp(conf_res2[0]), np.exp(conf_res2[1]), res2.pvalues[p+'_standardized_resid_cont']]


def plot_effect_sizes(results, file_to_save):
    rcParams.update({'font.size':12})

    fig, ax = plt.subplots(1, 4, figsize=(15, 5), sharex=True)
    
    groups = {'Ancestry': ['All', 'EUR', 'AFR', 'AMR', 'EAS', 'SAS'],
              'SIRE': ['All', 'White', 'Black', 'Asian', 'Hispanic', 'Other'],
             'BMI':['All','(0.0, 18.5]','(18.5, 24.9]','(24.9, 30.0]', '(30.0, 50.0]'],
             'Age/Sex':['All','(18, 45]','(45, 65]','(65, 100]','Female','Male']}
    
    def custom_ordering(group, categories):
        order_map = {category: idx for idx, category in enumerate(categories)}
        return sorted(group, key=lambda x: order_map.get(x, float('inf')))
    
    for i, (title, l) in enumerate(groups.items()):
        res = results[results.data_partition.isin(l)]
        if title == 'Ancestry' or title == 'SIRE':
            ordered_labels = custom_ordering(res['data_partition'], l)
            res1 = res.set_index('data_partition').loc[ordered_labels].reset_index()
            ax[i].set_yticks(np.arange(1.5, len(ordered_labels) * 5, 5))
            ax[i].set_yticklabels(ordered_labels[::-1])
        else:
            ax[i].set_yticks(np.arange(1.5, res.shape[0] * 5, 5))
            ax[i].set_yticklabels(res['data_partition'][::-1]) 
            res1 = res  # For other groups, no change in order
        
        ax[i].errorbar(res1['OR'], 
                     np.arange(0, res1.shape[0]*7, 7)[::-1],
                     xerr=[res1['OR'] - res1['lower_ci'], 
                           res1['upper_ci'] - res1['OR']],
                     ls='none',
                    marker='o', markersize=5, label='Standardized across all') 
        ax[i].axvline(1, linestyle='--', color='k')
        
        for y in np.arange(6, res1.shape[0]*7, 7):
            ax[i].axhline(y, linestyle='--', color='lightgray')
        
        ax[i].set_yticks(np.arange(2.5, res1.shape[0]*7, 7))
        ax[i].set_yticklabels(res1['data_partition'][::-1])

        ax[i].errorbar(res1['OR_resid'], 
                     np.arange(1, res1.shape[0]*7, 7)[::-1],
                     xerr=[res1['OR_resid'] - res1['lower_ci_resid'], 
                           res1['upper_ci_resid'] - res1['OR_resid']],
                     ls='none',
                    marker='o', markersize=5, color='C1', label='Contexts removed and standardized across all')  

        ax[i].errorbar(res1['OR_cont'], 
                     np.arange(2, res1.shape[0]*7, 7)[::-1],
                     xerr=[res1['OR_cont'] - res1['lower_ci_cont'], 
                           res1['upper_ci_cont'] - res1['OR_cont']],
                     ls='none',
                    marker='o', markersize=5, color='C2', label='Standardized within context') 

        ax[i].errorbar(res1['OR_cont_resid'], 
                     np.arange(3, res1.shape[0]*7, 7)[::-1],
                     xerr=[res1['OR_cont_resid'] - res1['lower_ci_cont_resid'], 
                           res1['upper_ci_cont_resid'] - res1['OR_cont_resid']],
                     ls='none',
                    marker='o', markersize=5, color='C3', label='Contexts removed and standardized within context')  
        ax[i].errorbar(res1['OR_gira'], 
                     np.arange(4, res1.shape[0]*7, 7)[::-1],
                     xerr=[res1['OR_gira'] - res1['lower_ci_gira'], 
                           res1['upper_ci_gira'] - res1['OR_gira']],
                     ls='none',
                    marker='o', markersize=5, color='C4', label='eMERGE PGS Ancestry Calibration')
        ax[i].set_title(title)
        ax[i].set_xlabel('OR (95% CI)')
    plt.tight_layout()

    ax[0].legend(bbox_to_anchor=(5.3, 1.3), ncol=3)
    plt.savefig(file_to_save, bbox_inches='tight', transparent=True, dpi=300)
    plt.show()


def calculate_effect_sizes_gira(data):
    final_results = {}
    PGS = ['BC_eMERGE_GIRA', 'CHD_eMERGE_GIRA', 'BC_PGS000507_GIRA', 'CHD_PGS003725_GIRA']
    for p in PGS:
        phe = p.split('_')[0]
        pgs_type = p.split('_')[1]
        data_to_analyze = data[pd.notnull(data[phe])]
        data_to_analyze[f"Age_binned"] = pd.cut(data[f"{phe}_age"], [18, 45, 65, 100]).astype(str)
        data_to_analyze[f"BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50]).astype(str)
    
        data_to_analyze[f"{phe}_{pgs_type}_gira"] = data_to_analyze[p]
        res_gira = get_glm_coef_emerge(data_to_analyze, phe, f"{phe}_{pgs_type}")
    
        res_temp = [['All']+res_gira]
    
        for context in ['ancestry', 'SIRE', 'BMI_binned', 'Age_binned', 'self_identified_sex']:
            for cont in np.unique(data_to_analyze[context]):
                if str(cont)=='nan' or str(cont)=='Unknown':
                    continue
                try: #If there are no cases in the context the model fails
                    data_to_analyze_cont = data_to_analyze[data_to_analyze[context]==cont]
                    res = get_glm_coef_emerge(data_to_analyze_cont, phe, f"{phe}_{pgs_type}") 
                    res_temp.append([cont]+res)
                except:
                    continue
    
        final_results[f"{phe}_{pgs_type}"] = pd.DataFrame(res_temp, columns=['data_partition', 'mean_PGS_gira', 'OR_gira', 'lower_ci_gira', 'upper_ci_gira', 'pval_gira'])
    return final_results
    
def calculate_effect_sizes(pgs_file, data_file, path_to_results):
    PGS = ['BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725']
    pgs = pd.read_csv(pgs_file, sep='\t')
    data = pd.read_csv(data_file, sep='\t')
    data = data.merge(pgs, on='ID')
    data = data[data.self_identified_sex.isin(['Male', 'Female'])]
    data = pd.concat((data,
                     pd.get_dummies(data['SIRE'], dtype=int),
                    pd.get_dummies(data['self_identified_sex'], dtype=int)),
    axis=1)

    gira_final_results = calculate_effect_sizes_gira(data)
    
    covs = ['White', 'Asian', 'Black', 'Hispanic', 'Other',
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 
            'Female'] 

    final_results = {}
    merged_results = {}
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
        res_temp = [['All']+res+[np.nan]*8]
    
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
                                                           'mean_PGS','OR', 'lower_ci', 'upper_ci', 'pval',
                                                           'mean_PGS_resid', 'OR_resid', 'lower_ci_resid', 'upper_ci_resid', 'pval_resid',
                                                           'mean_PGS_cont', 'OR_cont', 'lower_ci_cont', 'upper_ci_cont', 'pval_cont',
                                                           'mean_PGS_cont_resid','OR_cont_resid', 'lower_ci_cont_resid', 'upper_ci_cont_resid', 'pval_cont_resid'])
        merged_results[p] = pd.merge(gira_final_results[p], final_results[p], on = 'data_partition', how = 'inner')
        
        plot_effect_sizes(merged_results[p], join(path_to_results, p+'.pdf'))
        merged_results[p].to_csv(join(path_to_results, p+'.tsv'), sep='\t', index=None)  
         
        

def residualize(data, prs, covs):
    res = sm.OLS(data[prs], 
                sm.add_constant(data[covs]),
                missing='drop').fit().resid
    res = (res-res.mean())/res.std()
    return res
    
    
def format_PGS_for_R2_analysis(phenotype_file, pgs_file, path_to_results):
    PGS = ['BC_eMERGE', 'BC_PGS000507', 'CHD_eMERGE', 'CHD_PGS003725']
    pgs = pd.read_csv(pgs_file, sep='\t')
    data = pd.read_csv(phenotype_file, sep='\t')
    data = data.merge(pgs, on='ID')
    data = pd.concat((data,
                     pd.get_dummies(data['SIRE'], dtype=int),
                     pd.get_dummies(data['self_identified_sex'], dtype=int)),
                    axis=1)
    covs = ['White', 'Asian', 'Black', 'Hispanic', 'Other',
        'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8',
       'PC9', 'PC10', 'Female'] 
    pcs = [ 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
    for p in PGS:
        phe = p.split('_')[0]

        data_to_analyze = data[pd.notnull(data[phe])]
        data_to_analyze[f"Age_binned"] = pd.cut(data[f"{phe}_age"], [18, 45, 65, 100]).astype(str)
        data_to_analyze[f"BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50]).astype(str)
        covs_aux = covs+[f'{phe}_age', f'{phe}_BMI']
        data_to_analyze = data_to_analyze[(data_to_analyze.self_identified_sex.isin(['Female', 'Male']))&
                                         (data_to_analyze.Age_binned!='nan')&
                                         (data_to_analyze.BMI_binned!='nan')&
                                         (data_to_analyze.ancestry!='Unclassified')
                                         ]

        data_to_analyze['PRS_PC'] = residualize(data_to_analyze, p, pcs)
        data_to_analyze['PRS_PC_sex'] = residualize(data_to_analyze, p, pcs+['Female'])
        data_to_analyze['PRS_PC_BMI'] = residualize(data_to_analyze, p, pcs+[f'{phe}_BMI'])
        data_to_analyze['PRS_PC_age'] = residualize(data_to_analyze, p, pcs+[ f'{phe}_age'])
        data_to_analyze['PRS_PC_all'] = residualize(data_to_analyze, p, covs_aux)

        data_to_analyze = data_to_analyze[['ID', 'self_identified_sex', 'SIRE', 'ancestry', 
                                           'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 
                                           'PC8', 'PC9', 'PC10', phe, f"{phe}_age",f"{phe}_BMI",
                                           'Age_binned', 'BMI_binned', 'PRS_PC', 'PRS_PC_sex', 
                                           'PRS_PC_BMI', 'PRS_PC_age', 'PRS_PC_all']]
        data_to_analyze.columns = ['ID', 'self_identified_sex', 'SIRE', 'ancestry', 'PC1',
                                   'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9',
                                   'PC10', 'phenotype', 'age', 'bmi','Age_binned', 'BMI_binned',
                                   'PRS_PC', 'PRS_PC_sex', 'PRS_PC_BMI', 'PRS_PC_age', 'PRS_PC_all']
        data_to_analyze.to_csv(join(path_to_results, f"{p}.csv"), index=None)

        
if __name__=='__main__':
    fire.Fire()