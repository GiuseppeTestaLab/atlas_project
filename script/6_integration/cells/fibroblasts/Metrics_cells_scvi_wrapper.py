# Testing scib metrcis on fibroblasts cells with scvi implementation of scib 

#%%
import scanpy as sc
import pandas as pd
import scvi
from scib_metrics.benchmark import Benchmarker

import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'integration/cells/fibroblasts/'

adata_list={'scgen_batch_corr_celltypes_embeddings.h5ad':'scGen_HVG',
            'scgen_batch_corr_celltypes_embeddings_HDG.h5ad':'scGen_HDG',
            'scvi_batch_corr_celltypes_2500.h5ad':'scVI_HVG',
            'scanvi_batch_corr_celltypes_2500.h5ad':'scANVI_HVG',
            'scvi_batch_corr_celltypes_HDG.h5ad':'scVI_HDG',
            'scanvi_batch_corr_celltypes_HDG.h5ad':'scANVI_HDG'}

# Load the first AnnData object to store all embeddings
first_file = list(adata_list.keys())[0]
adata_first = sc.read(initDir + first_file)

# Function to get the correct embedding key based on the filename
def get_embedding_key(filename):
    if 'scgen' in filename.lower():
        return 'corrected_latent'
    elif 'scvi' in filename.lower():
        return 'X_scVI'
    elif 'scanvi' in filename.lower():
        return 'X_scANVI'
    else:
        raise ValueError(f"Unknown embedding type for file: {filename}")
    
# Iterate over each file, load the AnnData object, and extract the embeddings
for file, label in adata_list.items():
    adata = sc.read_h5ad(initDir + file)
    embedding_key = get_embedding_key(file)
    embeddings = adata.obsm[embedding_key]
    adata_first.obsm[label] = embeddings

adata_first

#%%
bm = Benchmarker(
    adata_first,
    batch_key="paper_ID",
    label_key="cell_types",
    embedding_obsm_keys = ['scGen_HVG', 'scGen_HDG', 'scVI_HVG', 'scANVI_HVG', 'scVI_HDG', 'scANVI_HDG'],
    n_jobs=-1,
)
bm.benchmark()

bm.plot_results_table(min_max_scale=False)

bm.plot_results_table(min_max_scale=True)

df = bm.get_results(min_max_scale=False)
print(df)

adata_first.write_h5ad(initDir + 'cells_fibroblasts_batch_corr_metrics.h5ad')