
#library imports
#%%
import scanpy as sc
import pandas as pd
import os
import numpy as np
import anndata

#initialize directory
dir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Olbrecht2021/'

#Read adata
adata = sc.read(dir + "olbrecht2021_filt_norm_nolog.h5ad")
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

#Filtering out genes not present in the other datasets, taking only the common genes
common_var_names = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/common_varnames_datasets.csv', index_col=0) #this file is derived from common_var_names.py

#Computing embeddings
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

adata.write(dir + 'olbrecht2021_embeddings.h5ad')


