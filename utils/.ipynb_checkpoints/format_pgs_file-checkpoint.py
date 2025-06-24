import pandas as pd
from os.path import join
import fire

def get_formatted_pgs_file(pgs_emerge_file, pgs_catalog_file, biobankname, path_to_results):
    emerge_pgs = pd.read_csv(pgs_emerge_file, sep = '\t')
    emerge_pgs = emerge_pgs[emerge_pgs.sampleset==f"{biobankname}"]
    emerge_pgs = emerge_pgs[['IID', 'PGS', 'SUM', 'Z_norm2']]
    emerge_pgs = emerge_pgs.pivot(index='IID', columns='PGS', values=['SUM', 'Z_norm2']).reset_index()
    emerge_pgs.columns = [f'{col[1]}_{col[0]}' for col in emerge_pgs.columns]
    emerge_pgs.rename(columns={
        '_IID':'ID',
        'breast_cancer_SUM': 'BC_eMERGE',
        'chd_SUM': 'CHD_eMERGE',
        'breast_cancer_Z_norm2': 'BC_eMERGE_GIRA',
        'chd_Z_norm2': 'CHD_eMERGE_GIRA',
    }, inplace=True)
    
    
    pgs_catalog = pd.read_csv(pgs_catalog_file, sep = '\t')
    pgs_catalog = pgs_catalog[pgs_catalog.sampleset==f"{biobankname}"]
    pgs_catalog = pgs_catalog[['IID', 'PGS', 'SUM', 'Z_norm2']]
    pgs_catalog = pgs_catalog.pivot(index='IID', columns='PGS', values=['SUM', 'Z_norm2']).reset_index()
    pgs_catalog.columns = [f'{col[1]}_{col[0]}' for col in pgs_catalog.columns]
    pgs_catalog.rename(columns={
        '_IID':'ID',
        'PGS000507_hmPOS_GRCh38_SUM': 'BC_PGS000507',
        'PGS003725_hmPOS_GRCh38_SUM': 'CHD_PGS003725',
        'PGS000507_hmPOS_GRCh38_Z_norm2': 'BC_PGS000507_GIRA',
        'PGS003725_hmPOS_GRCh38_Z_norm2': 'CHD_PGS003725_GIRA',
    }, inplace=True)
    
    pgs = emerge_pgs.merge(pgs_catalog, on = 'ID', how = 'inner')
    pgs.to_csv(join(path_to_results, "calculated_PGS.tsv"), index=None, sep='\t')

if __name__=='__main__':
    fire.Fire()