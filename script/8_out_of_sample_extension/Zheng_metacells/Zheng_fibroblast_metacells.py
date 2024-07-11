# Zheng fibroblasts metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/Zheng2023/Adata/'
destDir = '/group/testa/Project/OvarianAtlas/Zheng2023/Metacells/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/metacells/fibroblasts/"

## Load data
#%%
adata= sc.read("/group/testa/Project/OvarianAtlas/Zheng2023/Adata/zheng_fibroblast_filt_norm_nolog.h5ad")
genes = '/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv'

## Preprocessing
#%%
sc.pp.neighbors(adata)
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad
sc.pp.log1p(adata)
hdg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv',  index_col=0)
hdg[hdg.highly_variable]
adata.var['highly_variable']=hdg.highly_variable
adata.var.highly_variable = adata.var.highly_variable.fillna(False)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.tl.umap(adata)
#adata.write_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Adata/zheng_fibroblast_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'fibroblast_seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'fibroblast_seacells_assignment_hdg_patients.h5ad')
ad = create_mc_matrix(adata)
#%%
ad.obs = ad.obs.drop(columns=['n_genes', 'n_genes_by_counts', 
                    'total_counts', 'total_counts_mt', 'pct_counts_mt', 'EndothelialMSK', 'FibroblastsMSK', 'CancerMSK', 'HematopoieticMSK', 
                    'cell_labels_ratio', 'max', 'assignment', 'leiden-1.8'])

ad.write(destDir + 'fibroblasts_seacells_hdg_patients.h5ad')

#%%
## Compute embeddings and plot metacells
adata = sc.read(destDir + 'fibroblasts_seacells_hdg_patients.h5ad')
adata = preprocess_mc(adata, genes)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='zheng_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='zheng_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='zheng_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='zheng_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='zheng_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='zheng_seacells_HDG_anatomy.png')

#%%
adata.write(destDir + 'fibroblasts_seacells_hdg_patients_embeddings.h5ad')