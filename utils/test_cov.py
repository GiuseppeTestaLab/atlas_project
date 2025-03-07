#%%
import scanpy as sc
import os
import configparser
from sklearn.metrics import pairwise_distances
import numpy as np
#%%
config = configparser.ConfigParser()
config.read("config.ini")
rawPath = config.get("DEFAULT", "rawPath")
initDir = rawPath + "atlas_annotated/"

adata = sc.read_h5ad(initDir + "atlas_endothelial_filt_norm_nolog.h5ad")
#%%
ori_path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated_backup/atlas_annotated/atlas_endothelial_filt_norm_nolog.h5ad"
adata_ma = sc.read_h5ad(ori_path)

# %%
adata_array = adata.X.toarray()
adata_ma_array = adata_ma.X.toarray()
correlation_distance = pairwise_distances(adata_array, adata_ma_array, metric="correlation", n_jobs=-1)

# Convert back to correlation similarity if needed
correlation_similarity = 1 - correlation_distance
# %%
np.save("correlation_similarity.npy", correlation_similarity)