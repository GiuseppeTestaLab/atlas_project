
#%%
import numpy as np
from sklearn import metrics
import scanpy as sc
import pandas as pd

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'

#%%
adata = sc.read(initDir + "atlas_endothelial_embeddings.h5ad")

#%%
leidenTotal=[]
scores = {}
for i in np.arange(0.3, 2.2, 0.2):
    sc.tl.leiden(adata,resolution = i,key_added="leiden-{}".format(round(i,2)))

    # Step 1: Access UMAP coordinates from the AnnData object
    y_true = adata.obs['paper_ID'].to_numpy()

    # Step 2: Access the Leiden cluster labels from the AnnData object
    y_pred = adata.obs[f'leiden-{round(i,2)}'].to_numpy()  # Predicted labels (Leiden clusters)

    # Step 3: Compute the supervised clustering metrics
    adjusted_rand_index = metrics.adjusted_rand_score(y_true, y_pred)
    adjusted_mutual_info = metrics.adjusted_mutual_info_score(y_true, y_pred)
    fowlkes_mallows = metrics.fowlkes_mallows_score(y_true, y_pred)
    scores[i] = [adjusted_rand_index, adjusted_mutual_info, fowlkes_mallows]

    print(f'Resolution: {round(i,2)}')

pd_scores = pd.DataFrame(scores).T
# %%
