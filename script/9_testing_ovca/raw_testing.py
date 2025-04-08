#%%
import scanpy as sc
# %%
adata_target_path = "/group/testa/Project/OvarianAtlasTestStep0/raw_data/metacells_step0/fibroblasts/" + 'seacells_hdg_patients.h5ad'
adata = sc.read_h5ad(adata_target_path)
# %%
x = adata.X.todense()
# %%
xraw = adata.layers['raw'].todense()
# %%
ad_atlas = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")

# %%
ad_atlas
# %%
atlasRaw = ad_atlas.raw.X.todense()
atlasRaw = atlasRaw[:,ad_atlas.raw.var_names.isin(ad_atlas.var_names)]
# %%
atlas = ad_atlas.X
# %%
