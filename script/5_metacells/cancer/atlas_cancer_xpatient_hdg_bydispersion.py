# Atlas metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'
destDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'

## Load data
#%%
adata= sc.read(initDir + "atlas_cancer_filt_norm_nolog.h5ad")
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'
# meta = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/Metacell_gen/atlas_cancer_obsnames.csv', index_col=0)
# adata.obs = meta

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'atlas_cancer_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'seacells_assignment_hdg_patients.h5ad')
ad = create_mc_matrix(adata)

ad.obs.drop(columns=['ID', 'sample_name', 'cell_type', 'cell_subtype', 'sample_ID', 
                                    'patient_id', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 
                                    'pct_counts_mt', 'CancerMSK', 'EndothelialMSK', 'FibroblastsMSK', 'HematopoieticMSK', 
                                    'cell_labels_ratio', 'max', 'assignment', 'leiden-1.8'], inplace=True)

ad.write(destDir + 'seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/metacells/cancer/"

adata = sc.read(destDir + 'seacells_hdg_patients.h5ad')
adata = preprocess_mc(adata, genes)
#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='cancer_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='cancer_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='cancer_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='cancer_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='cancer_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='cancer_seacells_HDG_anatomy.png')

#%%
adata.write(destDir + 'seacells_hdg_patients_embeddings.h5ad')