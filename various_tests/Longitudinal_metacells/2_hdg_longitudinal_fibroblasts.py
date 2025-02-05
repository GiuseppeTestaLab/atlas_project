#%%
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)
from cell_labeller import assign_scores, actual_labeller, create_fibroblast_adata

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/Longitudinal/'
destDir = scriptsPath + '4_hdg/Tables/'

#%%
adata = sc.read(initDir + 'RnaAdataLongitudinalwithSubtypes_meta_cell_labelled.h5ad')

#%%
### I will take as var names the var names of this dataset
common_var_names = adata.var_names.tolist()
#%%
## Cell labelling strategy
# adata = assign_scores(adata)
# adata = actual_labeller(adata)
adata_fibroblasts = create_fibroblast_adata(adata)

adata_fibroblasts.write(initDir + 'longitudinal_adata_fibroblasts.h5ad')

#%%
dataDir = '/group/testa/Project/OvarianAtlas/'
dataName = ['longitudinal']

# common_var_names = pd.read_csv(scriptsPath + '4_hdg/Tables/common_varnames_datasets.csv', index_col=0)

dispersion_table = pd.DataFrame(index = common_var_names)
hvg_table = pd.DataFrame(index = common_var_names)

for j in dataName:
    path = (dataDir + j.capitalize() + "/{}_adata_fibroblasts.h5ad").format(j)
    print(path)
    adata_cancer = sc.read(path)
    dispersion_gene_xpatient = {}
    highly_variable_genes_per_patient = {}
    for patient_code in adata_cancer.obs["paper_ID"].unique():
        patient_anndata = adata_cancer[adata_cancer.obs["paper_ID"] == patient_code].copy()
        sc.pp.highly_variable_genes(patient_anndata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        highly_variable_genes_per_patient[patient_code] = patient_anndata.var.highly_variable
        dispersion_gene_xpatient[patient_code] = patient_anndata.var.dispersions_norm
    for i in dispersion_gene_xpatient:
        print(i)
        print(dispersion_gene_xpatient[i])
        dispersion_table[i] = dispersion_gene_xpatient[i]
        hvg_table[i] = highly_variable_genes_per_patient[i]

dispersion_table.to_csv(destDir + 'dispersion_table_fibroblasts_longitudinal.csv')
hvg_table.to_csv(destDir + 'hvg_table_fibroblasts_longitudinal.csv')

# 2_hvg_union_patients_dispersion_fibroblasts longitudinal
## Computing HDG by patient based on dispersion values

dir = scriptsPath + '4_hdg/Tables/'
#%%
dispersion_table = pd.read_csv(dir + 'dispersion_table_fibroblasts_longitudinal.csv', index_col=0)

#%%
dispersion_table

#%%
max_values_per_row = dispersion_table.max(axis=1)

print(max_values_per_row)

#%%
filtered_df = dispersion_table[dispersion_table > 1.40].dropna(how='all')  # 5035 genes

#%%
list_of_genes = filtered_df.index.tolist()
print(len(list_of_genes))

#%%
hdg_table = pd.DataFrame(index=dispersion_table.index)

#%%
hdg_table['highly_variable'] = hdg_table.index.isin(list_of_genes)

#%%
hdg_table.to_csv(dir + 'atlas_hdg_common_dispersion_patients_fibroblasts_longitudinal.csv')

#%%
hdg_true = pd.DataFrame(index=list_of_genes)
hdg_true['highly_variable'] = "True"
hdg_true.to_csv(dir + 'atlas_hdg_dispersion_patients_fibroblasts_longitudinal.csv')
