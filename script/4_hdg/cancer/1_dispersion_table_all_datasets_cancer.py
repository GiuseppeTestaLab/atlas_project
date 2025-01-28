# %%
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import sys
from cell_labeller import assign_scores, actual_labeller, create_cancer_adata
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)

# %%
# initialize directories
paths = pd.read_csv(
    scriptsPath+ "4_hdg/filepaths.csv",
    sep=";",
)

geistlinger = paths.Geistlinger[0]
loret = paths.Loret[0]
olbrecht = paths.Olbrecht[0]
qian = paths.Qian[0]
regner = paths.Regner[0]
ren = paths.Ren[0]
vasquez = paths.Vasquez[0]
xu = paths.Xu[0]
zhang = paths.Zhang[0]

# %%
## Geistlinger
adata = sc.read(geistlinger + "geistlinger2020_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(geistlinger + "geistlinger2020_adata_cancer.h5ad")

# %%
## Loret
adata = sc.read(loret + "loret2022_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(loret + "loret2022_adata_cancer.h5ad")

# %%
## Olbrecht
adata = sc.read(olbrecht + "olbrecht2021_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(olbrecht + "olbrecht2021_adata_cancer.h5ad")
# %%
## Qian
adata = sc.read(qian + "qian2020_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(qian + "qian2020_adata_cancer.h5ad")

# %%
## Regner
adata = sc.read(regner + "regner2021_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(regner + "regner2021_adata_cancer.h5ad")

# %%
## Ren
adata = sc.read(ren + "ren2022_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(ren + "ren2022_adata_cancer.h5ad")

# %%
## Vasquez
adata = sc.read(vasquez + "vasquez2022_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(vasquez + "vasquez2022_adata_cancer.h5ad")

# %%
## Xu
adata = sc.read(xu + "xu2022_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(xu + "xu2022_adata_cancer.h5ad")

# %%
## Zhang
adata = sc.read(zhang + "zhang2022_embeddings.h5ad")

## Cell labelling strategy
adata = assign_scores(adata)
adata = actual_labeller(adata)
adata_cancer = create_cancer_adata(adata)

adata_cancer.write(zhang + "zhang2022_adata_cancer.h5ad")

# %%
dataDir = rawPath + "original_anndata/"
dataName = [
    "geistlinger2020",
    "loret2022",
    "olbrecht2021",
    "qian2020",
    "regner2021",
    "ren2022",
    "vasquez2022",
    "xu2022",
    "zhang2022",
]

common_var_names = pd.read_csv(
    scriptsPath + "4_hdg/Tables/common_varnames_datasets.csv",
    index_col=0,
)

dispersion_table = pd.DataFrame(index=common_var_names.index)
hvg_table = pd.DataFrame(index=common_var_names.index)

for j in dataName:
    path = (dataDir + j.capitalize() + "/{}_adata_cancer.h5ad").format(j)
    print(path)
    adata_cancer = sc.read(path)
    dispersion_gene_xpatient = {}
    highly_variable_genes_per_patient = {}
    for patient_code in adata_cancer.obs["paper_ID"].unique():
        patient_anndata = adata_cancer[
            adata_cancer.obs["paper_ID"] == patient_code
        ].copy()
        sc.pp.highly_variable_genes(
            patient_anndata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )
        highly_variable_genes_per_patient[patient_code] = (
            patient_anndata.var.highly_variable
        )
        dispersion_gene_xpatient[patient_code] = patient_anndata.var.dispersions_norm
    for i in dispersion_gene_xpatient:
        print(i)
        print(dispersion_gene_xpatient[i])
        dispersion_table[i] = dispersion_gene_xpatient[i]
        hvg_table[i] = highly_variable_genes_per_patient[i]

dispersion_table.to_csv(
    scriptsPath + "4_hdg/Tables/dispersion_table_cancer.csv"
)
hvg_table.to_csv(
    scriptsPath + "4_hdg/Tables/hvg_table_cancer.csv"
)

# %%
