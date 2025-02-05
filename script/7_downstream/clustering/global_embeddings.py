# Creating a combined UMAP embedding for all the cell types
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
cancerDir = rawPath + 'downstream/clustering/cancer/'
fibroDir = rawPath + 'downstream/clustering/fibroblasts/'
endoDir = rawPath + 'downstream/clustering/endothelial/'
immuneDir = rawPath + 'downstream/clustering/immune/'
finalDir = rawPath + 'downstream/clustering/'
sc.settings.figdir = figPath + 'cluster_assignments/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
## Load the data
cancer = sc.read(cancerDir + 'tissues_combined_cellstates_embeddings.h5ad')
cancer

#%%
fibro = sc.read(fibroDir + 'tissues_combined_cellstates_embeddings.h5ad')
fibro

#%%
endo = sc.read(endoDir + 'tissues_combined_cellstates_embeddings.h5ad')
endo

#%%
immune = sc.read(immuneDir + 'tissues_combined_cellstates_embeddings.h5ad')
immune

#%%
adata_combined = cancer.concatenate(fibro, endo, immune, batch_key='major_celltypes', batch_categories=['CancerMSK', 'FibroblastsMSK', 'EndothelialMSK', 'HematopoieticMSK'])

#%%
# Extract and shift UMAP embeddings
umap1 = cancer.obsm['X_umap_shifted']
umap2 = fibro.obsm['X_umap_shifted'] + np.array([70, 0])  
umap3 = endo.obsm['X_umap_shifted'] + np.array([0, 60]) 
umap4 = immune.obsm['X_umap_shifted'] + np.array([70, 60]) 

#%%
umap_combined = np.vstack([umap1, umap2, umap3, umap4])
adata_combined.obsm['X_umap_shifted'] = umap_combined

#%%
# Plot the combined UMAP embeddings
sc.pl.embedding(adata_combined, color=['major_celltypes'], basis="X_umap_shifted", frameon=False)

#%%
adata_combined.write(finalDir + "major_celltypes_combined_embeddings.h5ad")