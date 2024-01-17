# Atlas endothelial metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
#%%
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc


#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'
destDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/endothelial/'

## Load data
#%%
adata = sc.read(initDir + "atlas_endothelial_filt_norm_nolog.h5ad")
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_endothelial.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'atlas_endothelial_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'seacells_assignment_hdg_patients.h5ad')

#%%
ad = create_mc_matrix(adata)
ad.write(destDir + 'seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells

#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/metacells/endothelial/"

#%%
adata = sc.read(destDir + 'seacells_hdg_patients.h5ad')
adata = preprocess_mc(adata, genes)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='endothelial_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='endothelial_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='endothelial_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='endothelial_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='endothelial_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='endothelial_seacells_HDG_anatomy.png')

#%%
adata.write(destDir + 'seacells_hdg_patients_embeddings.h5ad')