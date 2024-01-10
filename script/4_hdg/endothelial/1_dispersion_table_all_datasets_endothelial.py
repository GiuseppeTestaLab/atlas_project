#%%
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from cell_labeller import assign_scores, actual_labeller, create_endothelial_adata

#%%
#initialize directories
paths = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/filepaths.csv', sep = ';')

geistlinger = paths.Geistlinger[0]
loret = paths.Loret[0]
olbrecht = paths.Olbrecht[0]
qian = paths.Qian[0]
regner = paths.Regner[0]
ren = paths.Ren[0]
vasquez = paths.Vasquez[0]
xu = paths.Xu[0]
zhang = paths.Zhang[0]

#%%
## Geistlinger
adata = sc.read(geistlinger + 'geistlinger2020_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(geistlinger + 'geistlinger2020_adata_endothelial.h5ad')

#%%
## Loret
adata = sc.read(loret + 'loret2022_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(loret + 'loret2022_adata_endothelial.h5ad')

#%%
## Olbrecht
adata = sc.read(olbrecht + 'olbrecht2021_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(olbrecht + 'olbrecht2021_adata_endothelial.h5ad')
#%%
## Qian
adata = sc.read(qian + 'qian2020_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(qian + 'qian2020_adata_endothelial.h5ad')

#%%
## Regner
adata = sc.read(regner + 'regner2021_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(regner + 'regner2021_adata_endothelial.h5ad')

#%%
## Ren
adata = sc.read(ren + 'ren2022_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(ren + 'ren2022_adata_endothelial.h5ad')
#%%
## Vasquez
adata = sc.read(vasquez + 'vasquez2022_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(vasquez + 'vasquez2022_adata_endothelial.h5ad')

#%%
## Xu
adata = sc.read(xu + 'xu2022_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(xu + 'xu2022_adata_endothelial.h5ad')

#%%
## Zhang
adata = sc.read(zhang + 'zhang2022_embeddings.h5ad')

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_endothelial = create_endothelial_adata(adata)

adata_endothelial.write(zhang + 'zhang2022_adata_endothelial.h5ad')

#%%

tableDir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/"
dataDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/"
dataName = ['regner2021', 'geistlinger2020', 'qian2020', 'ren2022', 'loret2022', 'olbrecht2021', 'xu2022', 'zhang2022', 'vasquez2022']

common_var_names = pd.read_csv(tableDir + 'common_varnames_datasets.csv', index_col=0)

dispersion_table = pd.DataFrame(index = common_var_names.index)
hvg_table = pd.DataFrame(index = common_var_names.index)

for j in dataName:
    path = (dataDir + j.capitalize() + "/{}_adata_endothelial.h5ad").format(j)
    print(path)
    adata_endo = sc.read(path)
    dispersion_gene_xpatient = {}
    highly_variable_genes_per_patient = {}
    for patient_code in adata_endo.obs["paper_ID"].unique():
        patient_anndata = adata_endo[adata_endo.obs["paper_ID"] == patient_code].copy()
        sc.pp.highly_variable_genes(patient_anndata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        highly_variable_genes_per_patient[patient_code] = patient_anndata.var.highly_variable
        dispersion_gene_xpatient[patient_code] = patient_anndata.var.dispersions_norm
    for i in dispersion_gene_xpatient:
        print(i)
        print(dispersion_gene_xpatient[i])
        dispersion_table[i] = dispersion_gene_xpatient[i]
        hvg_table[i] = highly_variable_genes_per_patient[i]

dispersion_table.to_csv(tableDir + 'dispersion_table_endothelial.csv')
hvg_table.to_csv(tableDir + 'hvg_table_endothelial.csv')

# %%
