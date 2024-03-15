# Plotting cluster assignments in metacells

#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np

#%%
## Initialize directories
tissueDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/endothelial/'
sc.settings.figdir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/endothelial/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
## Load the data
primary = sc.read(tissueDir + 'adata_primary_embeddings.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings.h5ad')
# %%
primary.obs['cell_states'] = np.nan

col = {"0":"Immunoreactive_cells-T_cells", 
        "1":"Immunoreactive_cells-neutrophils", 
        "2":"Unknown_primary", 
        "3":"Cellular_metabolism", 
        "4":"Angiogenesis", 
        "5":"Unknown_primary", 
        "6":"Unknown_primary",
        "7":"RNA_metabolism", 
        "8":"Unknown_primary", 
        "9":"Immunoreactive_cells", 
        "10":"Cellular_metabolism", 
        "11":"Cycling_cells", 
        "12":"Unknown_primary"}

primary.obs['cell_states'] = primary.obs['leiden-0.51'].replace(col)

sc.pl.umap(primary, color=['cell_states'], frameon=False, save='_primary_cell_states.png')
sc.pl.umap(primary, color=['leiden-0.51'], frameon=False, save='_primary_leiden051.png')
sc.pl.umap(primary, color=['treatment'], frameon=False, save='_primary_treatment.png')
sc.pl.umap(primary, color=['phase'], frameon=False, save='_primary_cellcycle.png')

primary.write_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')

#%%
metastasis.obs['cell_states'] = np.nan

col = { "0":"Immunoreactive_cells-T_cells", 
        "1":"Immunoreactive_cells-neutrophils", 
        "2":"Angiogenesis", 
        "3":"Unknown_metastasis", 
        "4":"Immunoreactive_cells-B_cells", 
        "5":"Cellular_metabolism",
        "6":"Immunoreactive_cells-neutrophils", 
        "7":"Unknown_metastasis", 
        "8":"Immunoreactive_cells-T_cells",
        "9":"Immunoreactive_cells",
        "10":"Immunoreactive_cells-T_cells",
        "11":"Cycling_cells",
        "12":"RNA_metabolism",
        "13":"Immunoreactive_cells-neutrophils",
        "14":"Unknown_metastasis",
        "15":"Unknown_metastasis"} 

metastasis.obs['cell_states'] = metastasis.obs['leiden-0.51'].replace(col)

sc.pl.umap(metastasis, color=['cell_states'], frameon=False, save='_metastasis_cell_states.png')
sc.pl.umap(metastasis, color=['leiden-0.51'], frameon=False, save='_metastasis_leiden051.png')
sc.pl.umap(metastasis, color=['treatment'], frameon=False, save='_metastasis_treatment.png')
sc.pl.umap(metastasis, color=['phase'], frameon=False, save='_metastasis_cellcycle.png')

metastasis.write_h5ad(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')

#%%
ascites.obs['cell_states'] = np.nan

col = { "0":"Unknown_ascites", 
        "1":"Immunoreactive_cells-T_cells", 
        "2":"Immunoreactive_cells-neutrophils", 
        "3":"Immunoreactive_cells-neutrophils", 
        "4":"RNA_metabolism", 
        "5":"Immunoreactive_cells-B_cells", 
        "6":"Immunoreactive_cells-B_cells",
        "7":"Phagocytic_cells", 
        "8":"Angiogenesis"}

ascites.obs['cell_states'] = ascites.obs['leiden-0.31'].replace(col)

sc.pl.umap(ascites, color=['cell_states'], frameon=False, save='_ascites_cell_states.png')
sc.pl.umap(ascites, color=['leiden-0.31'], frameon=False, save='_ascites_leiden031.png')
sc.pl.umap(ascites, color=['leiden-0.21'], frameon=False, save='_ascites_leiden021.png')
sc.pl.umap(ascites, color=['treatment'], frameon=False, save='_ascites_treatment.png')
sc.pl.umap(ascites, color=['phase'], frameon=False, save='_ascites_cellcycle.png')

ascites.write_h5ad(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')
