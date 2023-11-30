# Concat atlas raw counts

import scanpy as sc


adata = sc.read('/group/testa/Project/OvarianAtlas/raw_data/original_counts/Geistlinger2020/Adata/geistlinger2020_rawcounts.h5ad')

data = ['/group/testa/Project/OvarianAtlas/raw_data/original_counts/Loret2022/Adata/loret2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Olbrecht2021/Adata/olbrecht2021_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Xu2022/Adata/xu2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Zhang2022/Adata/zhang2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Regner2021/Adata/regner2021_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Ren2022/Adata/ren2022_rawcounts.h5ad',
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Qian2020/Adata/qian2020_rawcounts.h5ad'
        '/group/testa/Project/OvarianAtlas/raw_data/original_counts/Vasquez2022/Adata/vasquez2022_rawcounts.h5ad']

for i in data:
    adata = adata.concatenate(sc.read(i), index_unique=None)

adata.obs = adata.obs.drop(columns=['batch'])

adata.write_h5ad('/group/testa/Project/OvarianAtlas/raw_data/original_counts/atlas_rawcounts.h5ad')

nb_fname = "atlas_rawcounts"

get_ipython().run_cell_magic('bash', '-s "$nb_fname"', 'jupyter nbconvert "$1".ipynb --to="python"\njupyter nbconvert "$1".ipynb --to="html"')