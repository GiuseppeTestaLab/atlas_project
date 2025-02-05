# Longitudinal fibroblasts metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/Longitudinal/'
destDir = '/group/testa/Project/OvarianAtlas/Longitudinal/Metacells/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/metacells/fibroblasts/"

## Load data
#%%
adata= sc.read(initDir + "longitudinal_fibroblast_filt_norm_nolog.h5ad")
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'longitudinal_fibroblasts_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'fibroblasts_seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'fibroblasts_seacells_assignment_hdg_patients.h5ad')
ad = create_mc_matrix(adata)
ad.write(destDir + 'fibroblasts_seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells

#%%
adata = sc.read(destDir + 'fibroblasts_seacells_hdg_patients.h5ad')
adata = preprocess_mc(adata, genes)

#%%
adata.write(destDir + 'fibroblasts_seacells_hdg_patients_embeddings.h5ad')
