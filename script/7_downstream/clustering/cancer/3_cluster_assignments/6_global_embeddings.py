# Creating a combined UMAP embedding for all the tissues
#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np

#%%
## Initialize directories
tissueDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/cancer/'
sc.settings.figdir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/cancer/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
## Load the data
primary = sc.read(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')
primary
#%%
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')
metastasis
#%%
ascites = sc.read(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')
ascites

#%%
adata_combined = primary.concatenate(metastasis, ascites, batch_key='tissue', batch_categories=['primary', 'metastasis', 'ascites'])

#%%
# Extract and shift UMAP embeddings
umap1 = primary.obsm['X_umap']
umap2 = metastasis.obsm['X_umap'] + np.array([25, 0])  
umap3 = ascites.obsm['X_umap'] + np.array([15, 25])  

#%%
umap_combined = np.vstack([umap1, umap2, umap3])
adata_combined.obsm['X_umap_shifted'] = umap_combined

#%%
# Plot the combined UMAP embeddings
sc.pl.embedding(adata_combined, color=['tissue'], basis="X_umap_shifted", frameon=False)

#%%
adata_combined.write(tissueDir + "tissues_combined_cellstates_embeddings.h5ad")