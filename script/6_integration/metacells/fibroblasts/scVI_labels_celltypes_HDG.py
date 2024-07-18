# Cancer seacells integration with scVI (without specifying the label key to preserve)
## Here I will train the model already within the space of the HDG

#%%
import scanpy as sc
import scvi
from rich import print
from scvi.model.utils import mde
import pandas as pd
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from integration import preprocess_scVI_genes

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/fibroblasts/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/fibroblasts/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv'

sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/group/testa/Project/OvarianAtlas/atlas_project/plots_def/integration/metacells/fibroblasts/"

#%%
adata = sc.read(initDir + "seacells_hdg_patients.h5ad")
adata.obs['tissue-treatment'] = adata.obs['tissue'].astype('str') + '_' + adata.obs['treatment'].astype('str')
adata

#%%
# inflammatory CAFs (iCAF)
sc.tl.score_genes(adata, ['IL6', 'CXCL12', 'VTN', 'CADM3', 'PLIN2', 'SERPINB2', 'KRT19', 'DES', 'CALB2', 'WT1', 'KRT7'], 
score_name = "iCAF", use_raw=False)

# matrix CAFs (mCAF)
sc.tl.score_genes(adata, ['ACTA2', 'COL11A1', 'MFAP5', 'SFRP2', 'ISLR', 'COL10A1', 'VIM', 'COL3A', 
                       'MMP11', 'PTHLH', 'FGF1', 'WNT7B', 'WNT2', 'TGFB3', 'THRC1', 'POSTN', 'VCAN', 'ZEB1'], 
score_name = "mCAF", use_raw=False)

# vascular cancer-associated fibroblasts CAFs (vCAF)
sc.tl.score_genes(adata, ['COX4I2', 'HIGD1B', 'PTP4A3', 'MCAM', 'PPP1R14A'], 
score_name = "vCAF", use_raw=False)

# STAR gene expressing CAFs (starCAF)
sc.tl.score_genes(adata, ['STAR', 'IGFBP5', 'TSPAN8', 'C7', 'ALDH1A1', 'LGR5'], 
score_name = "starCAF", use_raw=False)

#%%
adata.obs['cell_types'] = adata.obs[['iCAF', 'mCAF','vCAF', 'starCAF']].idxmax(axis=1)
adata.obs

batch = "paper_ID"
adata = preprocess_scVI_genes(adata, batch, genes)

#%%
scvi.model.SCVI.setup_anndata(adata, batch_key="paper_ID")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
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
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw = True)
adata.obs

sc.pl.umap(adata, color=['treatment'], frameon=False, save='_treatment_scVI_HDG.png')
sc.pl.umap(adata, color=['phase'], frameon=False, save='_cellcycle_scVI_HDG.png')
sc.pl.umap(adata, color=['tissue-treatment'], frameon=False, save='_tissue-treatment_scVI_HDG.png')
sc.pl.umap(adata, color=['paper_ID'], frameon=False, save='_patient_scVI_HDG.png')
sc.pl.umap(adata, color=['tissue'], frameon=False, save='_tissue_scVI_HDG.png')

#%%
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["paper_ID", "leiden"],
    frameon=False,
    ncols=1,
)

sc.pl.embedding(adata, basis="X_mde", color=["tissue"], frameon=False, ncols=1)

#%%
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scVI_celltypes_embeddings_HDG.h5ad')


# fibroblasts seacells integration with scANVI (specifying the label key to preserve)

## Since we’ve already trained an scVI model on our data, we will use it to initialize scANVI. When initializing scANVI, we provide it the labels_key. As scANVI can also be used for datasets with partially-observed annotations, we need to give it the name of the category that corresponds to unlabeled cells. As we have no unlabeled cells, we can give it any random name that is not the name of an exisiting cell type.
#%%
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="cell_types",
    unlabeled_category="Unknown",
)

#%%
lvae.train(max_epochs=20, n_samples_per_label=100)

#%%
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
## visualize the latent space
adata.obsm["X_mde_scanvi"] = mde(adata.obsm["X_scANVI"])

sc.pl.embedding(
    adata, basis="X_mde_scanvi", color=["cell_types"], ncols=1, frameon=False
)

## visualize data
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['cell_types'], frameon=False, save='_celltypes_scANVI_HDG.png')
sc.pl.umap(adata, color=['treatment'], frameon=False, save='_treatment_scANVI_HDG.png')
sc.pl.umap(adata, color=['phase'], frameon=False, save='_cellcycle_scANVI_HDG.png')
sc.pl.umap(adata, color=['tissue-treatment'], frameon=False, save='_tissue-treatment_scANVI_HDG.png')
sc.pl.umap(adata, color=['paper_ID'], frameon=False, save='_patient_scANVI_HDG.png')
sc.pl.umap(adata, color=['tissue'], frameon=False, save='_tissue_scANVI_HDG.png')

# %%
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scANVI_celltypes_embeddings_HDG.h5ad')
