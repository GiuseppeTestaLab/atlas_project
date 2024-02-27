# Cancer seacells integration with scVI (without specifying the label key to preserve)
## Here I will train the model already within the space of the HDG

#%%
import scanpy as sc
import scvi
from rich import print
from scvi.model.utils import mde
import pandas as pd

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

#%%
adata = sc.read(initDir + "seacells_hdg_patients.h5ad")
adata.obs['tissue-treatment'] = adata.obs['tissue'].astype('str') + '_' + adata.obs['treatment'].astype('str')
adata

#%%
hvg = pd.read_csv(genes,  index_col=0)
hvg[hvg.highly_variable]
adata.var['highly_variable'] = hvg.highly_variable
adata.var.highly_variable = adata.var.highly_variable.fillna(False)

### Here, adata.X should have log transformed scran normalized expression.
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.raw = adata # keep full dimension safe
adata = adata[:, adata.var.highly_variable]
adata.var_names # HDG genes

#%%
# adata.raw = adata  # keep full dimension safe
# sc.pp.highly_variable_genes(
#     adata,
#     flavor="seurat_v3",
#     n_top_genes=adata.var_names,
#     batch_key="paper_ID",
#     subset=True,
# )

#%%
adata = adata.copy()

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

# #%%
# sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
# sc.settings.figdir = "figures/"

#%%
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw = True)
adata.obs
sc.pl.umap(adata, color=['treatment', 'phase'], frameon=False)

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
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scVI_tissuetreat_embeddings_HDG.h5ad')


# Cancer seacells integration with scANVI (specifying the label key to preserve)

## Since we’ve already trained an scVI model on our data, we will use it to initialize scANVI. When initializing scANVI, we provide it the labels_key. As scANVI can also be used for datasets with partially-observed annotations, we need to give it the name of the category that corresponds to unlabeled cells. As we have no unlabeled cells, we can give it any random name that is not the name of an exisiting cell type.
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
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scANVI_tissuetreat_embeddings_HDG.h5ad')
# %%
