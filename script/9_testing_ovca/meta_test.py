#%%
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
# %%
adata_path = "/group/testa/Project/NeuroCOV/PBMC/seurat_PBMCs_neurocov/PBMC_DZNE_gex.h5ad"
adata = sc.read_h5ad(adata_path)
metadata_obs = pd.read_csv('/group/testa/Project/NeuroCOV/PBMC/seurat_PBMCs_neurocov/metadata_obs_26.03.2025.csv', sep = ' ')
metadata_obs
adata = adata[adata.obs.index.isin(metadata_obs.index), :]
adata.obs = adata.obs.merge(metadata_obs[['harm_snn_res.0.2', 'Annotation', 'gran_anno']], left_index = True, right_index = True)
adata
#%%
metadata = pd.read_csv('/group/testa/Users/alessia.valenti/NeuroCov/Repo/NeuroCOV/PBMC/NeuroCOVMetadata_sc_outer.csv')
metadata
# %%
metadata.set_index("u_sample_label", inplace = True)
# %%
meta_join = adata.obs.join(metadata, how="left", on="Sample_ID", validate="many_to_one")
adata.obs = meta_join
# %%
adata.write_h5ad("/group/testa/Users/vittorio.aiello/singularity/ovarian/cellxgene/PBMC_DZNE_gex_meta.h5ad")
# %%
harmony_umap = adata.obsm['H.UMAP']
harmony_umap = np.array(harmony_umap)
# %%
adata.obsm['X_umap_harmony'] = harmony_umap
# %%
umap = adata.obsm['UMAP']
umap = np.array(umap)
adata.obsm['X_umap'] = umap
#%%
umap_adt = pd.read_csv("/group/testa/Project/NeuroCOV/PBMC/seurat_PBMCs_neurocov/joint_umap_coords.csv",sep = ' ')
#%%
umap_adt = np.array(umap_adt)
adata.obsm['X_umap_adt'] = umap_adt

# %%
adata.write_h5ad("/group/testa/Users/vittorio.aiello/singularity/ovarian/cellxgene/PBMC_DZNE_gex_meta.h5ad")


# %%
