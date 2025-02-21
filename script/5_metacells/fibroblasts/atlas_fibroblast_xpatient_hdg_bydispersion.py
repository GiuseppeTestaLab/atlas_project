# Atlas metacells generation

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
destDir = rawPath + 'metacells/fibroblasts/'
if not os.path.exists(destDir):
    os.makedirs(destDir)
scriptsPath = config.get("DEFAULT", "scriptsPath")

## Load data
#%%
adata= sc.read(initDir + 'atlas_fibroblast_filt_norm_nolog.h5ad')
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'atlas_fibroblasts_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'seacells_assignment_hdg_patients_seed_1.h5ad')

# Creating metacell matrix
# #%%
# adata = sc.read(destDir + 'seacells_assignment_hdg_patients.h5ad')
# ad = create_mc_matrix(adata)
# ad.write(destDir + 'seacells_hdg_patients.h5ad')

# ## Compute embeddings and plot metacells
# #%%
# sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
# sc.settings.figdir = figPath + "metacells/fibroblasts/"

# #%%
# adata = sc.read(destDir + 'seacells_hdg_patients.h5ad')

# #%%
# # inflammatory CAFs (iCAF)
# sc.tl.score_genes(adata, ['IL6', 'CXCL12', 'VTN', 'CADM3', 'PLIN2', 'SERPINB2', 'KRT19', 'DES', 'CALB2', 'WT1', 'KRT7'], 
# score_name = "iCAF", use_raw=False)

# # matrix CAFs (mCAF)
# sc.tl.score_genes(adata, ['ACTA2', 'COL11A1', 'MFAP5', 'SFRP2', 'ISLR', 'COL10A1', 'VIM', 'COL3A', 
#                        'MMP11', 'PTHLH', 'FGF1', 'WNT7B', 'WNT2', 'TGFB3', 'THRC1', 'POSTN', 'VCAN', 'ZEB1'], 
# score_name = "mCAF", use_raw=False)

# # vascular cancer-associated fibroblasts CAFs (vCAF)
# sc.tl.score_genes(adata, ['COX4I2', 'HIGD1B', 'PTP4A3', 'MCAM', 'PPP1R14A'], 
# score_name = "vCAF", use_raw=False)

# # STAR gene expressing CAFs (starCAF)
# sc.tl.score_genes(adata, ['STAR', 'IGFBP5', 'TSPAN8', 'C7', 'ALDH1A1', 'LGR5'], 
# score_name = "starCAF", use_raw=False)

# #%%
# adata.obs['cell_types'] = adata.obs[['iCAF', 'mCAF','vCAF', 'starCAF']].idxmax(axis=1)
# adata.obs

# adata = preprocess_mc(adata, genes)
# #%%
# sc.pl.umap(adata, color=["treatment"], frameon=False, save='fibroblasts_seacells_HDG_treatm.png')
# sc.pl.umap(adata, color=["tissue"], frameon=False, save='fibroblasts_seacells_HDG_tissue.png')
# sc.pl.umap(adata, color=["dataset"], frameon=False, save='fibroblasts_seacells_HDG_dataset.png')
# sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='fibroblasts_seacells_HDG_patients.png')
# sc.pl.umap(adata, color=["phase"], frameon=False, save='fibroblasts_seacells_HDG_cellcycle.png')
# sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='fibroblasts_seacells_HDG_anatomy.png')
# sc.pl.umap(adata, color=["cell_types"], frameon=False, save='fibroblasts_seacells_HDG_celltypes.png')

# #%%
# adata.write(destDir + 'seacells_hdg_patients_embeddings.h5ad')
