#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import configparser
import anndata

# Read configuration file
config = configparser.ConfigParser()
config.read("config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")
rawPath = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/"
outdir="/group/testa/Project/OvarianAtlasTest/atlas_project/cellxgene/"
ovcaPath="/group/testa/Project/OvarianAtlas/atlas_project/cellxgene/OvCA_umap_tiled.h5ad"

#%%
ovca = sc.read(ovcaPath)

#%%
tissues=['cancer','fibroblasts','endothelial','immune']
"seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad"
adatas = []
for tissue in tissues:
    tissueRawDir = '{}metacells_backup/metacells/{}/'.format(rawPath,tissue)
    adataRaw = sc.read(tissueRawDir + 'seacells_hdg_patients_embeddings.h5ad')
    adataRaw.obs_names = adataRaw.obs_names.str.cat(adataRaw.obs.tissue.str.lower(), sep="-")
    obsm_to_removed=['X_pca', 'corrected_latent', 'latent']
    adatas.append(adataRaw)

#%%
adata_concat = adatas[0].concatenate(adatas[1], adatas[2], adatas[3], batch_key='major_celltypes', batch_categories=['CancerMSK', 'FibroblastsMSK', 'EndothelialMSK', 'HematopoieticMSK'])

common_obs = adata_concat.obs_names.intersection(ovca.obs_names)
assert len(common_obs) == adata_concat.shape[0] 
assert len(common_obs) == ovca.shape[0]
adata_concat = adata_concat[ovca.obs_names]
new_adata = anndata.AnnData(X=adata_concat.X, obs=ovca.obs, var=adata_concat.var, obsm=ovca.obsm, uns=ovca.uns, varm=adata_concat.varm)

#%%
# Plot the combined UMAP embeddings
#%%
sc.pl.embedding(new_adata, color=['02_tissue', "01_major_celltypes"], basis="X_umap_shifted", frameon=False)
#%%
new_adata.write_h5ad(outdir + 'ovca_plus_raw.h5ad')


# %%
