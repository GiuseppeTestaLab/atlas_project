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
tissueDir = rawPath + 'downstream/clustering/fibroblasts/'
sc.settings.figdir = figPath + 'cluster_assignments/fibroblasts/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
## Load the data
primary = sc.read(tissueDir + 'adata_primary_embeddings.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings.h5ad')
# %%
primary.obs['cell_states'] = np.nan

col = {"0":"Cellular_metabolism", 
        "1":"ECM_shaping_cells", 
        "2":"Unknown_primary", 
        "3":"Smooth_muscle_cells_development", 
        "4":"Immunoreactive_cells", 
        "5":"ECM_shaping_cells", 
        "6":"Cellular_metabolism",
        "7":"Protein_metabolism-cell_death", 
        "8":"Unknown_primary", 
        "9":"Epithelium_development", 
        "10":"Cycling_cells", 
        "11":"Angiogenesis", 
        "12":"Immunoreactive_cells",
        "13":"Cellular_metabolism",
        "14":"Response_to_stress-ROS"}

primary.obs['cell_states'] = primary.obs['leiden-0.41'].replace(col)

sc.pl.umap(primary, color=['cell_states'], frameon=False, save='_primary_cell_states.png')
sc.pl.umap(primary, color=['leiden-0.41'], frameon=False, save='_priamry_leiden041.png')
sc.pl.umap(primary, color=['treatment'], frameon=False, save='_primary_treatment.png')
sc.pl.umap(primary, color=['phase'], frameon=False, save='_primary_cellcycle.png')
sc.pl.umap(primary, color=['cell_types'], frameon=False, save='_primary_celltypes.png')

primary.write_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')

#%%
metastasis.obs['cell_states'] = np.nan

col = {"0":"Collagen_degradation", 
        "1":"Protein_catabolism", 
        "2":"Cellular_metabolism", 
        "3":"Unknown_metastasis", 
        "4":"Vascular_processes_regulation", 
        "5":"Immunoreactive_cells-T_cells",
        "6":"Cellular_metabolism", 
        "7":"ECM_shaping_cells", 
        "8":"Cellular_metabolism",
        "9":"ECM_shaping_cells",
        "10":"Angiogenesis",
        "11":"Unknown_metastasis",
        "12":"Epithelium_development-cell_division",
        "13":"Smooth_muscle_cells_development",
        "14":"Cycling_cells",
        "15":"Unknown_metastasis",
        "16":"Immunoreactive_cells-T_cells"} 

metastasis.obs['cell_states'] = metastasis.obs['leiden-0.51'].replace(col)

sc.pl.umap(metastasis, color=['cell_states'], frameon=False, save='_metastasis_cell_states.png')
sc.pl.umap(metastasis, color=['leiden-0.51'], frameon=False, save='_metastasis_leiden051.png')
sc.pl.umap(metastasis, color=['treatment'], frameon=False, save='_metastasis_treatment.png')
sc.pl.umap(metastasis, color=['phase'], frameon=False, save='_metastasis_cellcycle.png')
sc.pl.umap(metastasis, color=['cell_types'], frameon=False, save='_metastasis_celltypes.png')

metastasis.write_h5ad(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')

#%%
ascites.obs['cell_states'] = np.nan

col = {"0":"Extracellular_tissue_development", 
        "1":"Unknown_ascites", 
        "2":"Cellular_metabolism", 
        "3":"Cellular_metabolism-ECM", 
        "4":"Unknown_ascites", 
        "5":"Unknown_ascites", 
        "6":"Unknown_ascites",
        "7":"RNA_metabolism", 
        "8":"Unknown_ascites",
        "9":"Unknown_ascites",
        "10":"Cycling_cells",
        "11":"ECM_shaping_cells",
        "12":"Cellular_metabolism",
        "13":"Unknown_ascites",
        "14":"Angiogenesis",
        "15":"Unknown_ascites"}

ascites.obs['cell_states'] = ascites.obs['leiden-0.81'].replace(col)

sc.pl.umap(ascites, color=['cell_states'], frameon=False, save='_ascites_cell_states.png')
sc.pl.umap(ascites, color=['leiden-0.81'], frameon=False, save='_ascites_leiden081.png')
sc.pl.umap(ascites, color=['treatment'], frameon=False, save='_ascites_treatment.png')
sc.pl.umap(ascites, color=['phase'], frameon=False, save='_ascites_cellcycle.png')
sc.pl.umap(ascites, color=['cell_types'], frameon=False, save='_ascites_celltypes.png')

ascites.write_h5ad(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

