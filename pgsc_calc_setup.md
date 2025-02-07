### 1. PGSC_CALC setup

- Install pgsc_calc following instructions at https://pgsc-calc.readthedocs.io/en/latest/getting-started.html 
- Test that it runs on your computational environment with 
      `nextflow run pgscatalog/pgsc_calc -profile test,<docker|singularity|conda>`
- Download 1000 genomes reference panel for ancestry normalization
        `wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_1000G_v1.tar.zst`
- Test that it runs on your computational environment with 
    
       `nextflow run pgscatalog/pgsc_calc -profile test,docker \
        --run_ancestry path/to/reference/pgsc_HGDP+1kGP_v1.tar.zst`
        
More informations are available at https://pgsc-calc.readthedocs.io/en/latest/how-to/ancestry.html

### 2. Samplesheet formatting
Documentation is provided at https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html. In python this can be done as follows:

    import pandas as pd
    
    def generate_samplesheet():
    out = pd.DataFrame({
                'sampleset': [‘<biobankname>']*22 ,
                       'path_prefix':[f'/path/to/bfiles/chr{CHROM} ' 
                                      for CHROM in range(1,23,1)],
                'chrom':[CHROM  for CHROM in range(1, 23, 1)],
                'format': ['bfile']*22})
    out.to_csv(f'PGS_samplesheet.csv', index = False)


### 3. Run pgsc_calc
You can run PGS calculation using the script `2_compute_PGS.sh` and specifying the correct paths to files. More information on this are provided in the README. To merge together the pgs files produced by pgsc_calc run the following code:

    import pandas as pd
    
    emerge_pgs = pd.read_csv("path/to/pgsc_calc/pgs_file.txt.gz", sep = '\t')
    emerge_pgs = emerge_pgs[emerge_pgs.sampleset==‘<biobankname>']
    emerge_pgs = emerge_pgs[['IID', 'PGS', 'SUM']]
    emerge_pgs = emerge_pgs.pivot(index='IID', columns='PGS', values=['SUM']).reset_index()
    emerge_pgs.rename(columns={
        'IID':'ID',
        'breast_cancer_SUM': 'BC_eMERGE',
        'chd_SUM': 'CHD_eMERGE'}, inplace=True)

    pgs_catalog = pd.read_csv("path/to/pgsc_calc/pgs_file.txt.gz", sep = '\t')
    pgs_catalog = pgs_catalog[pgs_catalog.sampleset==‘<biobankname>']
    pgs_catalog = pgs_catalog[['IID', 'PGS', 'SUM']]
    pgs_catalog = pgs_catalog.pivot(index='IID', columns='PGS', values=['SUM']).reset_index()
    pgs_catalog.rename(columns={
        'IID':'ID',
        'PGS000507_hmPOS_GRCh38_SUM': 'BC_PGS000507',
        'PGS003725_hmPOS_GRCh38_SUM': 'CHD_PGS003725',
    }, inplace=True)

    pgs = emerge_pgs.merge(pgs_catalog, on = 'ID', how = 'inner')
    pgs.to_csv("calculated_PGS.tsv", index=None, sep='\t')