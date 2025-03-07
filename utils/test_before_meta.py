#%%
import scanpy as sc
import os
import configparser
# %%
config = configparser.ConfigParser()
config.read("config.ini")
rawPath = config.get("DEFAULT", "rawPath")
initDir = rawPath + "atlas_annotated/"

adata = sc.read_h5ad(initDir + "atlas_endothelial_filt_norm_nolog.h5ad")
#%%
ori_path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated_backup/atlas_annotated/atlas_endothelial_filt_norm_nolog.h5ad"
adata_ma = sc.read_h5ad(ori_path)
# %%
diff = adata_ma.X  - adata.X
diff = diff.toarray()
diff = diff.reshape(-1)
# %%
import seaborn as sns
sns.displot(diff)
# %%
import numpy as np
diff_by_gene = np.true_divide(diff.sum(1),(diff!=0).sum(1))
sns.displot(diff_by_gene)
# %%
