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
adata = sc.read(init_dir + 'Geistlinger2020/Adata/geistlinger2020_rawcounts.h5ad')
# #%%
# datasets=["loret2022",
#             "olbrecht2021",
#             "qian2020",
#             "regner2021",
#             "ren2022",
#             "vasquez2022",
#             "xu2022",
#             "zhang2022"]

# files_path = os.listdir(init_dir)
# listpath = [os.path.join(init_dir, file, 'Adata/') for file in files_path]

#%%
data = ['/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Adata/loret2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Olbrecht2021/Adata/olbrecht2021_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Xu2022/Adata/xu2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Zhang2022/Adata/zhang2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Regner2021/Adata/regner2021_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Ren2022/Adata/ren2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Adata/qian2020_rawcounts.h5ad'
        '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Vasquez2022/Adata/vasquez2022_rawcounts.h5ad']

for i in data:
    adata = adata.concatenate(sc.read(i), index_unique=None)

adata.obs = adata.obs.drop(columns=['batch'])

adata.write_h5ad('/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/atlas_rawcounts.h5ad')

# nb_fname = "atlas_rawcounts"

# get_ipython().run_cell_magic('bash', '-s "$nb_fname"', 'jupyter nbconvert "$1".ipynb --to="python"\njupyter nbconvert "$1".ipynb --to="html"')