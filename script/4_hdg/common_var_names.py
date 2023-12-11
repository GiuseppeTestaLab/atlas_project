# Find all the common var_names across datasets

#%%
import scanpy as sc
import pandas as pd
import numpy as np

#%%
adata1 = sc.read("/group/testa/Project/OvarianAtlas/Qian2020/Adata/qian2020_filt_norm_nolog.h5ad")
adata2 = sc.read("/group/testa/Project/OvarianAtlas/Regner2021/Adata/regner2021_filt_norm_nolog.h5ad")
adata3 = sc.read("/group/testa/Project/OvarianAtlas/Ren2022/Adata/ren2022_filt_norm_nolog.h5ad")
adata4 = sc.read("/group/testa/Project/OvarianAtlas/Geistlinger2020/Adata/geistlinger2020_filt_norm_nolog.h5ad")
adata5 = sc.read("/group/testa/Project/OvarianAtlas/Loret2022/Adata/loret2022_filt_norm_nolog.h5ad")
adata6 = sc.read("/group/testa/Project/OvarianAtlas/Olbrecht2021/Adata/olbrecht2021_filt_norm_nolog.h5ad")
adata7 = sc.read("/group/testa/Project/OvarianAtlas/Xu2022_CCR/Adata/xu2022_filt_norm_nolog.h5ad")
adata8 = sc.read("/group/testa/Project/OvarianAtlas/Zhang2022/Adata/zhang2022_filt_norm_nolog.h5ad")
adata9 = sc.read("/group/testa/Project/OvarianAtlas/Vasquez2022/Adata/vasquez2022_filt_norm_nolog.h5ad")

#%%
common_var_names = (adata1.var_names) & (adata2.var_names) & (adata3.var_names) & (adata4.var_names) & (adata5.var_names) & (adata6.var_names) & (adata7.var_names) & (adata8.var_names) & (adata9.var_names) 

pd.DataFrame(index = common_var_names).to_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/common_varnames_datasets.csv')

