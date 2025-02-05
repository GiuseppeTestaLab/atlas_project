# Longitudinal cancer metacells generation

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
figPath = config.get("DEFAULT", "figPath")

sys.path.insert(1, utilsPath)
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/Longitudinal/'
destDir = '/group/testa/Project/OvarianAtlas/Longitudinal/Metacells/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "metacells/cancer/"

## Load data
#%%
adata= sc.read(initDir + "longitudinal_cancer_filt_norm_nolog.h5ad")
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'longitudinal_cancer_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'cancer_seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'cancer_seacells_assignment_hdg_patients.h5ad')
ad = create_mc_matrix(adata)
ad.write(destDir + 'cancer_seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells

#%%
adata = sc.read(destDir + 'cancer_seacells_hdg_patients.h5ad')
adata = preprocess_mc(adata, genes)

#%%
adata.write(destDir + 'cancer_seacells_hdg_patients_embeddings.h5ad')