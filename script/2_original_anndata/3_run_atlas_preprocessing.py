import numpy as np
import pandas as pd
import scanpy as sc
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read('../../utils/config.ini')
rowPath = config.get('DEFAULT', 'rowPath')



#%%
adata = sc.read(rowPath + 'original_counts/atlas_rawcounts.h5ad')
#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
#%%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save = 'pct_counts_mt.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = 'n_genes_by_counts.png')
#%%
adata = adata[adata.obs.n_genes_by_counts < 7000, :]
adata = adata[adata.obs.pct_counts_mt < 30, :]
#%%
sc.pp.normalize_total(adata, target_sum=1e4)
#%%
adata.write(rowPath + 'original_anndata/atlas_filt_norm_nolog_def.h5ad')
#%%
sc.pp.log1p(adata)
#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#%%
sc.pl.highly_variable_genes(adata)
#%%
adata.raw = adata
#%%
adata = adata[:, adata.var.highly_variable]
#%%
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#%%
sc.pp.scale(adata, max_value=10)
#%%
sc.tl.pca(adata, svd_solver='arpack')
#%%
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
#%%
sc.tl.umap(adata)
#%%
adata.write(rowPath + 'original_anndata/atlas_embeddings.h5ad')