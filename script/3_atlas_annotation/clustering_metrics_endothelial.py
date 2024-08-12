
#%%
import numpy as np
from sklearn import metrics
import scanpy as sc
import pandas as pd

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'

#%%
adata = sc.read(initDir + "atlas_endothelial_embeddings.h5ad")
# %%
leidenTotal=[]
for i in np.arange(0.3, 1.5, 0.1):
    sc.tl.leiden(adata,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

#%%
for i in leidenTotal:
    sc.pl.umap(adata, color=i, frameon=False)

#%%
# Step 1: Access UMAP coordinates from the AnnData object
y_true = adata.obs['paper_ID'].to_numpy()

# Step 2: Access the Leiden cluster labels from the AnnData object
y_pred = adata.obs['leiden-1.0'].to_numpy()  # Predicted labels (Leiden clusters)

#%%
# Step 3: Compute the supervised clustering metrics
adjusted_rand_index = metrics.adjusted_rand_score(y_true, y_pred)
adjusted_mutual_info = metrics.adjusted_mutual_info_score(y_true, y_pred)
fowlkes_mallows = metrics.fowlkes_mallows_score(y_true, y_pred)

#%%
# Print the results
print(f'Adjusted Rand Index: {adjusted_rand_index}')
print(f'Adjusted Mutual Information: {adjusted_mutual_info}')
print(f'Fowlkes-Mallows Index: {fowlkes_mallows}')

# Adjusted Rand Index: 0.1663130205776094
# Adjusted Mutual Information: 0.368091311281982
# Fowlkes-Mallows Index: 0.2058092921044652

#%%
# Save the results to a file
finalDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotation/'

with open('clustering_metrics_endothelial.txt', 'w') as file:
    file.write(f"Adjusted Rand Index: {adjusted_rand_index}\n")
    file.write(f"Adjusted Mutual Information: {adjusted_mutual_info}\n")
    file.write(f"Fowlkes-Mallows Index: {fowlkes_mallows}\n")

# %%
leidenTotal=[]
for i in np.arange(1.4, 1.9, 0.1):
    sc.tl.leiden(adata,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

# %%
leidenTotal=[]
for i in np.arange(1.9, 2.3, 0.1):
    sc.tl.leiden(adata,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))
# %%
