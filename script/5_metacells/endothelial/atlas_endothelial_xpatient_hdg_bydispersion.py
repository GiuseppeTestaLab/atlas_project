# Atlas endothelial metacells generation

## Import libraries
#%%

import scanpy as sc
import sys
import configparser
import os
import numpy as np
np.random.seed(17)
# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

sys.path.insert(1, utilsPath)
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc # type: ignore


#%%
initDir = rawPath + 'atlas_annotated/'
destDir = rawPath + 'metacells/endothelial/'
if not os.path.exists(destDir):
    os.makedirs(destDir)
scriptsPath = config.get("DEFAULT", "scriptsPath")

## Load data
#%%
adata = sc.read(initDir + "atlas_endothelial_filt_norm_nolog.h5ad")
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_endothelial.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'atlas_endothelial_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'seacells_assignment_hdg_patients_seed_1.h5ad')

# Creating metacell matrix
# #%%
# adata = sc.read(destDir + 'seacells_assignment_hdg_patients.h5ad')

# #%%
# ad = create_mc_matrix(adata)
# ad.write(destDir + 'seacells_hdg_patients.h5ad')

# ## Compute embeddings and plot metacells

# #%%
# sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
# sc.settings.figdir = figPath + "metacells/endothelial/"

# #%%
# adata = sc.read(destDir + 'seacells_hdg_patients.h5ad')
# adata = preprocess_mc(adata, genes)

# #%%
# sc.pl.umap(adata, color=["treatment"], frameon=False, save='endothelial_seacells_HDG_treatm.png')
# sc.pl.umap(adata, color=["tissue"], frameon=False, save='endothelial_seacells_HDG_tissue.png')
# sc.pl.umap(adata, color=["dataset"], frameon=False, save='endothelial_seacells_HDG_dataset.png')
# sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='endothelial_seacells_HDG_patients.png')
# sc.pl.umap(adata, color=["phase"], frameon=False, save='endothelial_seacells_HDG_cellcycle.png')
# sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='endothelial_seacells_HDG_anatomy.png')

# #%%
# adata.write(destDir + 'seacells_hdg_patients_embeddings.h5ad')