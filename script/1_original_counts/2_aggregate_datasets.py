# Concat atlas raw counts

#%%
## Imports 
import scanpy as sc
import os
import subprocess
import logging

## Directories
init_dir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/'
#%%
## Aggregating datasets

## load the first
adata = sc.read(init_dir + 'Geistlinger2020/Adata/geistlinger2020_rawcounts.h5ad')

#%%
data = [init_dir + 'Loret2022/Adata/loret2022_rawcounts.h5ad',
        init_dir + 'Olbrecht2021/Adata/olbrecht2021_rawcounts.h5ad',
        init_dir + 'Xu2022/Adata/xu2022_rawcounts.h5ad',
        init_dir + 'Zhang2022/Adata/zhang2022_rawcounts.h5ad',
        init_dir + 'Regner2021/Adata/regner2021_rawcounts.h5ad',
        init_dir + 'Ren2022/Adata/ren2022_rawcounts.h5ad',
        init_dir + 'Qian2020/Adata/qian2020_rawcounts.h5ad',
        init_dir + 'Vasquez2022/Adata/vasquez2022_rawcounts.h5ad']

for i in data:
    adata = adata.concatenate(sc.read(i), index_unique=None)

adata.obs = adata.obs.drop(columns=['batch'])

adata.write_h5ad(init_dir + 'atlas_rawcounts.h5ad')