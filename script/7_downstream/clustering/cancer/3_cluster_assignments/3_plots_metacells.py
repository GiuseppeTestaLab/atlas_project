# Plotting cluster assignments in metacells

#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")

#%%
## Initialize directories
tissueDir = rawPath + 'downstream/clustering/cancer/'
sc.settings.figdir = figPath + 'cluster_assignments/cancer/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
## Load the data
primary = sc.read(tissueDir + 'adata_primary_embeddings.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings.h5ad')
# %%
primary.obs['cell_states'] = np.nan

col = {"0":"Unknown_primary", 
        "1":"Cellular_metabolism-extracellular_signaling", 
        "2":"Organelles_organization-metabolism", 
        "3":"Unknown_primary", 
        "4":"Unknown_primary", 
        "5":"Cycling_cells", 
        "6":"Cycling_cells",
        "7":"Extracellular_signaling", 
        "8":"Immunoreactive_cells", 
        "9":"Organelles_organization-cell_cycle", 
        "10":"Unknown_primary", 
        "11":"Organelles_organization-cell_movement", 
        "12":"Ciliated_cancer_cells",
        "13":"Immunoreactive_cells",
        "14":"Organelles_organization-metabolism",
        "15":"Unknown_primary",
        "16":"INF_mediated_signaling",
        "17":"Cellular_metabolism",
        "18":"Unknown_primary",
        "19":"Unknown_primary",
        "20":"Unknown_primary"}

primary.obs['cell_states'] = primary.obs['leiden-0.41'].replace(col)

sc.pl.umap(primary, color=['cell_states'], frameon=False, save='_primary_cell_states.png')
sc.pl.umap(primary, color=['leiden-0.41'], frameon=False, save='_priamry_leiden041.png')
sc.pl.umap(primary, color=['treatment'], frameon=False, save='_primary_treatment.png')
sc.pl.umap(primary, color=['phase'], frameon=False, save='_primary_cellcycle.png')

primary.write_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')

#%%
metastasis.obs['cell_states'] = np.nan

col = { "0":"Unknown_metastasis", 
        "1":"Cellular_metabolism", 
        "2":"Unknown_metastasis", 
        "3":"Unknown_metastasis", 
        "4":"Cycling_cells", 
        "5":"Unknown_metastasis",
        "6":"Immunoreactive_cells", 
        "7":"Organelles_organization-cell_cycle", 
        "8":"Cycling_cells",
        "9":"Response_to_extracellular_signals",
        "10":"ECM_shaping_cells",
        "11":"Ciliated_cancer_cells",
        "12":"Cellular_metabolism",
        "13":"Extracellular_signaling-immune_cells",
        "14":"Unknown_metastasis",
        "15":"RNA_metabolism"}

metastasis.obs['cell_states'] = metastasis.obs['leiden-0.41'].replace(col)

sc.pl.umap(metastasis, color=['cell_states'], frameon=False, save='_metastasis_cell_states.png')
sc.pl.umap(metastasis, color=['leiden-0.41'], frameon=False, save='_metastasis_leiden041.png')
sc.pl.umap(metastasis, color=['treatment'], frameon=False, save='_metastasis_treatment.png')
sc.pl.umap(metastasis, color=['phase'], frameon=False, save='_metastasis_cellcycle.png')

metastasis.write_h5ad(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')

#%%
ascites.obs['cell_states'] = np.nan

col = { "0":"Organelles_organization-metabolism", 
        "1":"Cycling_cells", 
        "2":"Unknown_ascites", 
        "3":"Unknown_ascites", 
        "4":"Unknown_ascites", 
        "5":"Unknown_ascites", 
        "6":"Cellular_metabolism",
        "7":"Response_to_extracellular_signals", 
        "8":"Unknown_ascites", 
        "9":"Response_to_stress", 
        "10":"Organelles_organization-metabolism"}

ascites.obs['cell_states'] = ascites.obs['leiden-0.31'].replace(col)

sc.pl.umap(ascites, color=['cell_states'], frameon=False, save='_ascites_cell_states.png')
sc.pl.umap(ascites, color=['leiden-0.31'], frameon=False, save='_ascites_leiden031.png')
sc.pl.umap(ascites, color=['treatment'], frameon=False, save='_ascites_treatment.png')
sc.pl.umap(ascites, color=['phase'], frameon=False, save='_ascites_cellcycle.png')

ascites.write_h5ad(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

# %%
