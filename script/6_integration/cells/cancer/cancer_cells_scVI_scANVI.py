# Adata atlas cancer integration scVI

#%%
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
from rich import print
from scvi.model.utils import mde

#%%
adata = sc.read('/group/testa/Project/OvarianAtlas/atlas_cancer_seacells_assignment.h5ad')
adata

#%%
# adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

#%%
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2500,
    batch_key="paper_ID",
    subset=True,
)

#%%
scvi.model.SCVI.setup_anndata(adata, batch_key="paper_ID")

#%%
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

#%%
vae.train()

#%%
adata.obsm["X_scVI"] = vae.get_latent_representation()

#%% 
## visualize the latent space
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)
## visualize data
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)

#%%
## embeddings visualization
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "figures/"

#%%
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata.obs
sc.pl.umap(adata, color=['treatment', 'phase'], frameon=False)


#%%
#%%
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["paper_ID", "leiden"],
    frameon=False,
    ncols=1,
)

sc.pl.embedding(adata, basis="X_mde", color=["tissue"], frameon=False, ncols=1)

#%%#%%
adata.write_h5ad('/group/testa/Project/OvarianAtlas/Integrated_data/atlas_cancer_batch_corr_scVI.h5ad')

## Integration with scANVI for label key preservation (Treatment)

#%% 
adata.obs['tissue-treatment'] = adata.obs['tissue'].astype('str') + '_' + adata.obs['treatment'].astype('str')
#%%
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="tissue-treatment",
    unlabeled_category="Unknown",
)

#%%
lvae.train(max_epochs=20, n_samples_per_label=100)

#%%
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
## visualize the latent space
adata.obsm["X_mde_scanvi"] = mde(adata.obsm["X_scANVI"])

sc.pl.embedding(
    adata, basis="X_mde_scanvi", color=["tissue-treatment"], ncols=1, frameon=False
)

## visualize data
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['tissue-treatment', 'paper_ID'], frameon=False)

# %%
adata.write_h5ad('/group/testa/Project/OvarianAtlas/Integrated_data/atlas_cancer_batch_corr_scANVI.h5ad')