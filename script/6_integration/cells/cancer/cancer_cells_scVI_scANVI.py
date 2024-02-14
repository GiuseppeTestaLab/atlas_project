# Adata atlas cancer integration scVI

#%%
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
from rich import print
from scvi.model.utils import mde
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
#from integration import 

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/cells/cancer/'

#%%
adata = sc.read(initDir + 'seacells_assignment_hdg_patients.h5ad')
adata
batch = "paper_ID"
genes = 2500

#%%
scvi.model.SCVI.setup_anndata(adata, batch_key=batch)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()

#%%
adata.obsm["X_scVI"] = vae.get_latent_representation()

#%% 
## visualize the latent space (optional)
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