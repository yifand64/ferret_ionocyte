## Grab the rare airway cell profiles from the Human Lung Cell Atlas. 
## Used to compare ferret rare cells to human for John Engelhardt collab. 
## Written by adam, september 2022

## Load and initialize scanpy
import numpy as np
import scanpy as sc
sc.settings.verbosity = 3


# Simekka et al (Human Lung Cell Atlas, 2022)
# Downloaded from 
# https://beta.fastgenomics.org/datasets/detail-dataset-427f1eee6dd44f50bae1ab13f0f3c6a9#Files
input_file = "/Users/ahaber/Dropbox/HSPH/projects/nasal_mvc/analysis/191218_scRNASeq/adam/revision/human_lca/data/HLCA_v1.h5ad"
adata = sc.read(input_file)


pd.crosstab(adata.obs['ann_level_3'], adata.obs['anatomical_region_level_3']).to_csv('cell_counts.txt', sep="\t")
pd.crosstab(adata.obs['ann_level_4'], adata.obs['anatomical_region_level_3']).to_csv('cell_counts_al4.txt', sep="\t")

# these are the same:
rcells = adata[adata.obs['ann_level_4'].isin(['Tuft', 'Ionocyte', 'Neuroendocrine'])]
rcells = adata[adata.obs['ann_level_3'].isin(['Rare'])]

rcells.write_csvs("rarecells", skip_data=False)


### Write out some data to use in R
sc.get.obs_df(rcells, keys=["sample", "study", "subject_ID", 'sex',
       'ethnicity', 'mixed_ethnicity', 'smoking_status', 'BMI', 'condition',
       'subject_type', 'sample_type', 'single_cell_platform', 
       'sequencing_platform', 'cell_ranger_version', 'fresh_or_frozen',
       'dataset', 'anatomical_region_level_1', 'anatomical_region_level_2',
       'anatomical_region_level_3', 'anatomical_region_highest_res', 'age',
       'ann_highest_res', 'n_genes', 'log10_total_counts', 'mito_frac',
       'ribo_frac', 'ann_finest_level', 'ann_level_1',
       'ann_level_2', 'ann_level_3', 'ann_level_4', 'ann_level_5',
       'ann_coarse_for_GWAS_and_modeling', "CYSLTR1", "CYSLTR2", "OXGR1", "LTC4S", "TRPM5", "EPCAM", "IL25", "LRMP", "IL17RB", "POU2F3", "PTGS1", "ALOX5", "ALOX5AP"]).to_csv('selected_metadata.txt', sep="\t")