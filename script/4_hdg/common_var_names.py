# Find all the common var_names across datasets

#%%
import scanpy as sc
import pandas as pd
import numpy as np

#%%
init_dir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/'
out_dir = '/home/marta.sallese/ov_cancer_atlas/atlas_project/4_hdg/Tables/'

adata1 = sc.read("Qian2020/Adata/qian2020_filt_norm_nolog.h5ad")
adata2 = sc.read("Regner2021/Adata/regner2021_filt_norm_nolog.h5ad")
adata3 = sc.read("Ren2022/Adata/ren2022_filt_norm_nolog.h5ad")
adata4 = sc.read("Geistlinger2020/Adata/geistlinger2020_filt_norm_nolog.h5ad")
adata5 = sc.read("Loret2022/Adata/loret2022_filt_norm_nolog.h5ad")
adata6 = sc.read("Olbrecht2021/Adata/olbrecht2021_filt_norm_nolog.h5ad")
adata7 = sc.read("Xu2022/Adata/xu2022_filt_norm_nolog.h5ad")
adata8 = sc.read("Zhang2022/Adata/zhang2022_filt_norm_nolog.h5ad")
adata9 = sc.read("Vasquez2022/Adata/vasquez2022_filt_norm_nolog.h5ad")

#%%
common_var_names = (adata1.var_names) & (adata2.var_names) & (adata3.var_names) & (adata4.var_names) & (adata5.var_names) & (adata6.var_names) & (adata7.var_names) & (adata8.var_names) & (adata9.var_names) 

pd.DataFrame(index = common_var_names).to_csv('common_varnames_datasets.csv')

