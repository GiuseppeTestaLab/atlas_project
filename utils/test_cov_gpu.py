#%%
import scanpy as sc
import os
import configparser
from sklearn.metrics import pairwise_distances
import cupy as cp
import numpy as np
# %%

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
adata_array = adata.X.toarray().astype(np.float32)
adata_ma_array = adata_ma.X.toarray().astype(np.float32)
#correlation_distance = pairwise_distances(adata_array, adata_ma_array, metric="correlation", n_jobs=-1)

# Convert back to correlation similarity if needed
#correlation_similarity = 1 - correlation_distance
# %%


def fast_correlation_gpu(expr1, expr2):
    """
    Compute Pearson correlation using GPU with CuPy.
    expr1, expr2: CuPy arrays with shape (N, G), where N is number of cells, G is number of genes.
    """
    # Normalize: zero mean, unit variance for each row (cell)
    expr1 = (expr1 - cp.mean(expr1, axis=1, keepdims=True)) / cp.std(expr1, axis=1, keepdims=True)
    expr2 = (expr2 - cp.mean(expr2, axis=1, keepdims=True)) / cp.std(expr2, axis=1, keepdims=True)

    # Compute correlation using dot product
    corr_matrix = cp.dot(expr1, expr2.T) / expr1.shape[1]
    
    return corr_matrix

expr1 = cp.asarray(adata_array)
expr2 = cp.asarray(adata_ma_array)
# Compute Pearson correlation matrix using CuPy
corr_matrix_gpu = fast_correlation_gpu(expr1, expr2)

# Convert to NumPy (if needed for further processing)
corr_matrix_cpu = cp.asnumpy(corr_matrix_gpu)
del corr_matrix_gpu
# %%
mean = cp.mean(expr1, axis=1, keepdims=True)
# %%
a = [1,2,3,4,5,6,7]
a = cp.asarray(a)
m = cp.mean(a)
# %%
def fast_correlation_gpu_chunked(expr1, expr2, chunk_size):
    """
    Compute Pearson correlation using GPU with CuPy in chunks to limit VRAM usage.
    expr1, expr2: CuPy arrays with shape (N, G), where N is number of cells, G is number of genes.
    chunk_size: Number of rows to process at a time.
    """
    n_chunks = int(np.ceil(expr1.shape[0] / chunk_size))
    corr_matrix = np.zeros((expr1.shape[0], expr2.shape[0]), dtype=cp.float32)
    expr2 = cp.asarray(expr2, dtype=cp.float32)
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, expr1.shape[0])
        expr1_chunk = cp.asarray(expr1[start:end], dtype=cp.float32)
        # Normalize: zero mean, unit variance for each row (cell)
        expr1_chunk = (expr1_chunk - cp.mean(expr1_chunk, axis=1, keepdims=True)) / cp.std(expr1_chunk, axis=1, keepdims=True)
        expr2 = (expr2 - cp.mean(expr2, axis=1, keepdims=True)) / cp.std(expr2, axis=1, keepdims=True)
        
        # Compute correlation using dot product
        corr_matrix[start:end] = cp.asnumpy(cp.dot(expr1_chunk, expr2.T) / expr1_chunk.shape[1])
    
    return corr_matrix

#%%
chunk_size = 10000  # Adjust chunk size based on available VRAM
#%%
# Compute Pearson correlation matrix using CuPy in chunks
corr_matrix_gpu = fast_correlation_gpu_chunked(expr1, expr2, chunk_size)

# Convert to NumPy (if needed for further processing)
corr_matrix_cpu = cp.asnumpy(corr_matrix_gpu)

# %%
ad_slice = adata_array[:1000,:]
expr1 = cp.asarray(ad_slice, dtype=cp.float32)
expr2 = cp.asarray(adata_ma_array, dtype=cp.float32)

# %%
