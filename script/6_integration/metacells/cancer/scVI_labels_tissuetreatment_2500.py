# Cancer seacells integration with scVI (without specifying the label key to preserve)

#%%
import scanpy as sc
import scvi
from rich import print
from scvi.model.utils import mde
import pandas as pd
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")
CCGenes = config.get("DEFAULT", "CCGenes")

sys.path.insert(1, utilsPath)
from integration import preprocess_scVI

#%%
initDir = rawPath + 'metacells/cancer/'
outDir = rawPath + 'integration/metacells/cancer/'

sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "integration/metacells/cancer/"

#%%
ad = sc.read(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')
ad
batch = "paper_ID"
genes = 2500
ad = preprocess_scVI(ad, batch, genes)

#%%
scvi.model.SCVI.setup_anndata(ad, batch_key=batch)
vae = scvi.model.SCVI(ad, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()

#%%
ad.obsm["X_scVI"] = vae.get_latent_representation()

#%% 
## visualize the latent space (optional)
sc.pp.neighbors(ad, use_rep="X_scVI")
sc.tl.leiden(ad)
## visualize data
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)

#%%
## embeddings visualization
ad.obsm["X_mde"] = mde(ad.obsm["X_scVI"])

#%%
cell_cycle_genes = [x.strip() for x in open(CCGenes)]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in ad.var_names]
sc.tl.score_genes_cell_cycle(ad, s_genes=s_genes, g2m_genes=g2m_genes, use_raw = True)
ad.obs
sc.pl.umap(ad, color=['treatment'], frameon=False, save='_treatment_scVI_2500.png')
sc.pl.umap(ad, color=['phase'], frameon=False, save='_cellcycle_scVI_2500.png')
sc.pl.umap(ad, color=['tissue-treatment'], frameon=False, save='_tissue-treatment_scVI_2500.png')
sc.pl.umap(ad, color=['paper_ID'], frameon=False, save='_patient_scVI_2500.png')
sc.pl.umap(ad, color=['tissue'], frameon=False, save='_tissue_scVI_2500.png')

#%%
sc.pl.embedding(
    ad,
    basis="X_mde",
    color=["paper_ID", "leiden"],
    frameon=False,
    ncols=1,
)

sc.pl.embedding(ad, basis="X_mde", color=["tissue"], frameon=False, ncols=1)

#%%
ad.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scVI_tissuetreat_embeddings_2500.h5ad')


# Cancer seacells integration with scANVI (specifying the label key to preserve)

## Since we’ve already trained an scVI model on our data, we will use it to initialize scANVI. When initializing scANVI, we provide it the labels_key. As scANVI can also be used for datasets with partially-observed annotations, we need to give it the name of the category that corresponds to unlabeled cells. As we have no unlabeled cells, we can give it any random name that is not the name of an exisiting cell type.
#%%
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=ad,
    labels_key="tissue-treatment",
    unlabeled_category="Unknown",
)

#%%
lvae.train(max_epochs=20, n_samples_per_label=100)

#%%
ad.obsm["X_scANVI"] = lvae.get_latent_representation(ad)
## visualize the latent space
ad.obsm["X_mde_scanvi"] = mde(ad.obsm["X_scANVI"])

sc.pl.embedding(
    ad, basis="X_mde_scanvi", color=["tissue-treatment"], ncols=1, frameon=False
)

## visualize data
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)
sc.pl.umap(ad, color=['treatment'], frameon=False, save='_treatment_scANVI_2500.png')
sc.pl.umap(ad, color=['phase'], frameon=False, save='_cellcycle_scANVI_2500.png')
sc.pl.umap(ad, color=['tissue-treatment'], frameon=False, save='_tissue-treatment_scANVI_2500.png')
sc.pl.umap(ad, color=['paper_ID'], frameon=False, save='_patient_scANVI_2500.png')
sc.pl.umap(ad, color=['tissue'], frameon=False, save='_tissue_scANVI_2500.png')

# %%
ad.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scANVI_tissuetreat_embeddings_2500.h5ad')
# %%
