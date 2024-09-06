# Find all the common var_names across datasets

#%%
import scanpy as sc
import pandas as pd
import numpy as np

#%%
init_dir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/'
out_dir = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/'

adata1 = sc.read(init_dir + "Qian2020/Qian2020_filt_norm_nolog.h5ad")
adata2 = sc.read(init_dir + "Regner2021/Regner2021_filt_norm_nolog.h5ad")
adata3 = sc.read(init_dir + "Ren2022/Ren2022_filt_norm_nolog.h5ad")
adata4 = sc.read(init_dir + "Geistlinger2020/Geistlinger2020_filt_norm_nolog.h5ad")
adata5 = sc.read(init_dir + "Loret2022/Loret2022_filt_norm_nolog.h5ad")
adata6 = sc.read(init_dir + "Olbrecht2021/Olbrecht2021_filt_norm_nolog.h5ad")
adata7 = sc.read(init_dir + "Xu2022/Xu2022_filt_norm_nolog.h5ad")
adata8 = sc.read(init_dir + "Zhang2022/Zhang2022_filt_norm_nolog.h5ad")
adata9 = sc.read(init_dir + "Vasquez2022/Vasquez2022_filt_norm_nolog.h5ad")

#%%
common_var_names = (adata1.var_names) & (adata2.var_names) & (adata3.var_names) & (adata4.var_names) & (adata5.var_names) & (adata6.var_names) & (adata7.var_names) & (adata8.var_names) & (adata9.var_names) 

pd.DataFrame(index = common_var_names).to_csv(out_dir + 'common_varnames_datasets.csv')

