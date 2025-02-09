
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
from scipy import sparse
from anndata import AnnData
from typing import Optional, Union

import numpy as np
from sklearn.neighbors import KNeighborsClassifier

#%%
ooseDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/out_of_sample_extension_backup/out_of_sample_extension/immune/'
mcDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream_backup/downstream/clustering/immune/'

#%%
## Setting fig parameteres
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/oose/immune/"

#%%
adata = sc.read_h5ad(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat.h5ad')
adata

#%%
adata.obs

#%%
primary = sc.read_h5ad(mcDir + 'adata_primary_embeddings_cellstates.h5ad')
ascites = sc.read_h5ad(mcDir + 'adata_ascites_embeddings_cellstates.h5ad')
metastasis = sc.read_h5ad(mcDir + 'adata_metastasis_embeddings_cellstates.h5ad')

#%%
# Extract cell states
primary_states = primary.obs[['cell_states']].rename(columns={'cell_states': 'primary_state'})
ascites_states = ascites.obs[['cell_states']].rename(columns={'cell_states': 'ascites_state'})
metastasis_states = metastasis.obs[['cell_states']].rename(columns={'cell_states': 'metastasis_state'})

#%%
# Merge with adata.obs
adata.obs = adata.obs.join(primary_states, how='left').join(ascites_states, how='left').join(metastasis_states, how='left')

#%%
# Create a final 'cell_states' column combining all sources, prioritizing primary > ascites > metastasis
adata.obs['cell_states'] = adata.obs['primary_state'].combine_first(adata.obs['ascites_state']).combine_first(adata.obs['metastasis_state'])


#%%
sc.pl.umap(adata, color='cell_states')

#%%
sc.pp.neighbors(adata, use_rep="latent_corrected")
sc.tl.umap(adata)
sc.pl.umap(adata, color='cell_states')

#%%
atlas_adata = adata[adata.obs["reference_map"] == "reference"]
new_dataset_adata = adata[adata.obs["reference_map"] != "reference"]

#%%
# Extract latent space embeddings
X_train = atlas_adata.obsm["latent_corrected"]
y_train = atlas_adata.obs["cell_states"].to_numpy()

X_test = new_dataset_adata.obsm["latent_corrected"]

#%%
# Train kNN classifier
knn = KNeighborsClassifier(n_neighbors=5, metric="cosine")  # Use cosine similarity for high-dimensional data
knn.fit(X_train, y_train)

#%%
# Predict labels for new dataset
new_dataset_adata.obs["predicted_cell_states"] = knn.predict(X_test)

#%%
sc.pl.umap(new_dataset_adata, color=["predicted_cell_states"], frameon=False, save="_oose_predicted_cell_states.png")

#%%
# Save adata
adata.write_h5ad(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat_predicted_cellstates.h5ad')

#%%
