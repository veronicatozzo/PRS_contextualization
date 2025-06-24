import pandas as pd 
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt 
from matplotlib import rcParams 
from os.path import join
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D

import warnings
import fire

warnings.simplefilter('ignore')
def get_glm_coef_gira(data, phe, p):
    results = smf.glm(
                formula=f'{phe} ~ {p}',
                data=data,
                family=sm.families.Binomial(),
                missing='drop'
            ).fit()
    ORs = np.exp(results.params[p])
    conf = results.conf_int().loc[p]
    lower_ci = np.exp(conf[0])
    upper_ci = np.exp(conf[1])
    pval = results.pvalues[p]
    return [ORs, lower_ci, upper_ci, pval]

def get_plot_within_context(results, file_to_save):
    rcParams.update({'font.size':12})
    fig, ax = plt.subplots(2, 5, figsize=(18, 10), 
                           sharex=False, gridspec_kw={'hspace': 0.5, 'wspace': 0.7})
    
    groups = {'EUR':['All', 'Asian', 'Black', 'Hispanic', 'Other', 'White',
                       '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]',
                       '(18, 45]', '(45, 65]', '(65, 100]', 'Female', 'Male'],
             'AFR':['All', 'Asian', 'Black', 'Hispanic', 'Other', 'White',
                       '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]',
                       '(18, 45]', '(45, 65]', '(65, 100]', 'Female', 'Male'],
             'AMR':['All', 'Asian', 'Black', 'Hispanic', 'Other', 'White',
                       '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]',
                       '(18, 45]', '(45, 65]', '(65, 100]', 'Female', 'Male'],
             'EAS':['All', 'Asian', 'Black', 'Hispanic', 'Other', 'White',
                       '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]',
                       '(18, 45]', '(45, 65]', '(65, 100]', 'Female', 'Male'],
             'SAS':['All', 'Asian', 'Black', 'Hispanic', 'Other', 'White',
                       '(0.0, 18.5]', '(18.5, 24.9]', '(24.9, 30.0]', '(30.0, 50.0]',
                       '(18, 45]', '(45, 65]', '(65, 100]', 'Female', 'Male']}
    
    for i, (title, l) in enumerate(groups.items()):
        trans1 = Affine2D().translate(0.0, +2.5) + ax[0,i].transData
        trans2 = Affine2D().translate(0.0, +0.5) + ax[0,i].transData
        trans3 = Affine2D().translate(0.0, +2.5) + ax[1,i].transData
        trans4 = Affine2D().translate(0.0, +0.5) + ax[1,i].transData
        
        res1 = results['BC_eMERGE_GIRA'][title]
        res1 = res1[(res1.Strata.isin(l)) & (res1.upper_ci < 15)]
        res2 = results['BC_PGS000507_GIRA'][title]
        res2 = res2[(res2.Strata.isin(l)) & (res2.upper_ci < 15)]
    
        ax[0,i].errorbar(res1['OR'], 
                         np.arange(0, res1.shape[0]*5, 5)[::-1],
                         xerr=[res1['OR'] - res1['lower_ci'], 
                               res1['upper_ci'] - res1['OR']],
                         ls='none',
                         marker='o', markersize=5, label='BC_eMERGE', transform=trans1) 
        ax[0,i].errorbar(res2['OR'], 
                     np.arange(0, res2.shape[0]*5, 5)[::-1],
                     xerr=[res2['OR'] - res2['lower_ci'], 
                         res2['upper_ci'] - res2['OR']],
                     ls='none',
                    marker='o', markersize=5, label='BC_PGS000507', transform=trans2) 
    
        ax[0,i].axvline(1, linestyle='--', color='k')  
        ax[0,i].set_title(title)
        ax[0,i].set_xlabel('OR (95% CI)')
        for y in np.arange(4, res1.shape[0]*5, 5):
            ax[0,i].axhline(y, linestyle='--', color='lightgray')
        ax[0,i].set_yticks(np.arange(1.5, res1.shape[0]*5, 5))
        ax[0,i].set_yticklabels(res1['Strata'][::-1])
    
        res3 = results['CHD_eMERGE_GIRA'][title]
        res3 = res3[(res3.Strata.isin(l)) & (res3.upper_ci < 15)]
        res4 = results['CHD_PGS003725_GIRA'][title]
        res4 = res4[(res4.Strata.isin(l)) & (res4.upper_ci < 15)]
    
        ax[1,i].errorbar(res3['OR'], 
                         np.arange(0, res3.shape[0]*5, 5)[::-1],
                         xerr=[res3['OR'] - res3['lower_ci'], 
                               res3['upper_ci'] - res3['OR']],
                         ls='none',
                         marker='o', markersize=5, label='CHD_eMERGE', transform=trans3) 
        ax[1,i].errorbar(res4['OR'], 
                     np.arange(0, res4.shape[0]*5, 5)[::-1],
                     xerr=[res4['OR'] - res4['lower_ci'], 
                         res4['upper_ci'] - res4['OR']],
                     ls='none',
                    marker='o', markersize=5, label='CHD_PGS003725', transform=trans4) 
    
        ax[1,i].axvline(1, linestyle='--', color='k')  
        ax[1,i].set_title(title)
        ax[1,i].set_xlabel('OR (95% CI)')
        for y in np.arange(4, res3.shape[0]*5, 5):
            ax[1,i].axhline(y, linestyle='--', color='lightgray')
        ax[1,i].set_yticks(np.arange(1.5, res3.shape[0]*5, 5))
        ax[1,i].set_yticklabels(res3['Strata'][::-1])
        
    # plt.tight_layout()
    ax[0,0].legend(bbox_to_anchor=(5, 1.3), ncol=2)
    ax[1,0].legend(bbox_to_anchor=(5, 1.3), ncol=2)
    plt.subplots_adjust(hspace=0.5, wspace= 0.7)
    plt.savefig(file_to_save, bbox_inches='tight', transparent=True, dpi=300)

def variability_within_ancestry(pgs_file, data_file, path_to_results):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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

            anc_results = {}
            ancestry = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
            for anc in ancestry:
                phe = p.split('_')[0]
                data_to_analyze = data[data.ancestry == anc]
                data_to_analyze = data_to_analyze[pd.notnull(data_to_analyze[phe])]
                data_to_analyze[f"Age_binned"] = pd.cut(data[f"{phe}_age"], [18, 45, 65, 100]).astype(str)
                data_to_analyze[f"BMI_binned"] = pd.cut(data[f"{phe}_BMI"], [0, 18.5, 24.9, 30, 50]).astype(str)
            
                data_to_analyze = data_to_analyze[(data_to_analyze.self_identified_sex.isin(['Female', 'Male']))&
                                             (data_to_analyze.Age_binned!='nan')&
                                             (data_to_analyze.BMI_binned!='nan')&
                                             (data_to_analyze.ancestry!='Unclassified')
                                             ]
                res = get_glm_coef_gira(data_to_analyze, phe, p)         
                res_temp = [['All']+res]
            
                for context in ['SIRE', 'BMI_binned', 'Age_binned', 'self_identified_sex']:
                    for cont in np.unique(data_to_analyze[context]):
                        if str(cont)=='nan' or str(cont)=='Unknown':
                            continue
                        try: #If there are no cases in the context the model fails
                            data_to_analyze_cont = data_to_analyze[data_to_analyze[context]==cont]
                            res_context = get_glm_coef_gira(data_to_analyze_cont, phe, p) 
                            res_temp.append([cont]+res_context) 
                        except:
                            continue
                anc_results[anc] = pd.DataFrame(res_temp, columns=['Strata', 'OR', 'lower_ci', 'upper_ci', 'pval'])
            final_results[p] = anc_results
            
            get_plot_within_context(final_results[p], join(path_to_results, p+'_var_within_anc.pdf'))
            
            for anc, df in final_results[p].items():
                filename = f"{p}_var_within_{anc}.tsv"
                df.to_csv(join(path_to_results, filename), sep='\t', index=None)

            
            # final_results[p].to_csv(join(path_to_results, p+'_var_within_anc.tsv'), sep='\t', index=None)  

if __name__=='__main__':
    fire.Fire()