# library imports
# %%
import scanpy as sc
import pandas as pd
import os
import numpy as np
import anndata
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")
rawPath = config.get("DEFAULT", "rawPath")

# initialize directory
dir = rawPath + "original_anndata/Xu2022/"

# Read adata
adata = sc.read(dir + "Xu2022_filt_norm_nolog.h5ad")
adata.obs["dataset"] = adata.obs.paper_ID.str.split("_").str[0]

# Filtering out genes not present in the other datasets, taking only the common genes
common_var_names = pd.read_csv(
    rawPath + "original_anndata/common_varnames_datasets.csv", index_col=0
)

# Computing embeddings
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

adata.write(dir + "xu2022_embeddings.h5ad")
