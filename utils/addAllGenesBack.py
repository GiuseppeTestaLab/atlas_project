#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")
rawPath = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/"
outdir="/group/testa/Project/OvarianAtlasTest/atlas_project/cellxgene/"
#%%

tissues=['cancer','fibroblasts','endothelial', "immune"]
"seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad"
adatas = []
for tissue in tissues:
    tissueDir = '{}/integration_backup/integration/metacells/{}/'.format(rawPath,tissue)
    if(tissue == 'fibroblasts' or tissue == 'immune'):
       adata = sc.read(tissueDir + "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad")
    else:
        adata = sc.read(tissueDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad')
    tissueRawDir = '{}metacells_backup/metacells/{}/'.format(rawPath,tissue)
    adataRaw = sc.read(tissueRawDir + 'seacells_hdg_patients_embeddings.h5ad')
    assert (adataRaw.obs_names == adata.obs_names).all()

    obsm_to_removed=['X_pca', 'corrected_latent', 'latent']
    for obsm in obsm_to_removed:
        adata.obsm.pop(obsm)


    adataRaw.obsm = adata.obsm
    adataRaw.obs = adata.obs
    adataRaw.write_h5ad(outdir + '{}_raw_with_scGen_batch_corrected_hdg_projection.h5ad'.format(tissue))
    adatas.append(adataRaw)

#%%
adatas[0].obsm["X_umap_shifted"] = adatas[0].obsm["X_umap"]
adatas[1].obsm["X_umap_shifted"] = adatas[1].obsm["X_umap"] + np.array([40, 0])
adatas[2].obsm["X_umap_shifted"] = adatas[2].obsm["X_umap"] + np.array([-10, 40])
adatas[3].obsm["X_umap_shifted"] = adatas[3].obsm["X_umap"] + np.array([40, 40])
adata_combined = adatas[0].concatenate(adatas[1], adatas[2], adatas[3],
                                        batch_key="major_celltypes", batch_categories=["cancer", "fibroblasts", "endothelial", "immune"],
                                        index_unique="_")
adata_combined.obsm.pop("X_umap")




#%%
# Plot the combined UMAP embeddings
#%%
sc.pl.embedding(adata_combined, color=['tissue'], basis="X_umap_shifted", frameon=False)
sc.pl.embedding(adata_combined, color=['major_celltypes'], basis="X_umap_shifted", frameon=False)
#%%
adata_combined.write_h5ad(outdir + 'combined_raw_with_scGen_batch_corrected_hdg_projection.h5ad'.format(tissue))

# %%
