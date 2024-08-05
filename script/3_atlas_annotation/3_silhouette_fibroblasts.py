# Compute Silhouette index for fibroblasts

#%%
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import scanpy as sc
import pandas as pd

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'

#%%
adata = sc.read(initDir + "atlas_fibroblasts_embeddings.h5ad")

X = adata.X
labels = adata.obs['paper_ID']
score = silhouette_score(X, labels)
print(f"Silhouette Score: {score}")

with open('silhouette_score_fibroblasts.txt', 'w') as file:
    file.write(f"Silhouette Score: {score}\n")

