# Immune seacells integration with scVI (without specifying the label key to preserve)
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
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/immune/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/immune/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_immune.csv'

sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/group/testa/Project/OvarianAtlas/atlas_project/plots_def/integration/metacells/immune/"

#%%
adata = sc.read(initDir + "seacells_hdg_patients.h5ad")
adata.obs['tissue-treatment'] = adata.obs['tissue'].astype('str') + '_' + adata.obs['treatment'].astype('str')
adata

#%%
# Plasma cells markers
sc.tl.score_genes(adata, ['IGKC','IGHG1','CD79A','IGHG2','IGLC2','IGLC3','IGHG3','IGHG4','JCHAIN','MZB1','XBP1'], 
score_name = "Plasma_cells", use_raw=False)

# T cells markers
sc.tl.score_genes(adata, ['CD2','CD3D','TRAC','GZMA','NKG7','CD3E','CD3G','CD4','TCF7','CD8A','PRF1','GZMB','CCL5','CCL4','IL32','CD52'], 
score_name = "T_cells", use_raw=False)

# Mast cells markers
sc.tl.score_genes(adata, ['KIT','CPA3','CTSG','MS4A2','TPSAB1','TPSB2','HPGD','HPGDS','GATA2'], 
score_name = "Mast_cells", use_raw=False)

# B cells markers
sc.tl.score_genes(adata, ['MS4A1', 'CD79A', 'CD19', 'BANK1', 'IGKC', 'IGHM'], score_name = "B_cells", use_raw=False)

# Myeloid cells markers
sc.tl.score_genes(adata, 
['CD14','FCER1G','FCGR3A','LYZ','CTSS','CD33','CD68','CD163','ITGAX','ITGAM','CD4','MRC1',
'VSIG4','SPP1','APOE','C1QA','C1QB','C1QC','APOC1','FTL','S100A9','TYROBP','AIF1','CD74','PSAP','CTSB'], 
score_name = "Myeloid_cells", use_raw=False)

# Dendritic cells markers
sc.tl.score_genes(adata, 
['IL3RA','IRF7','IRF8','GZMB','CD4','CLEC4C','JCHAIN',
'PTGDS','PLAC8','PLD4','TCF4','BCL11A','GPR183','CCDC50','LILRA4','TSPAN13','CLIC3','MPEG1'], 
score_name = "Dendritic_cells", use_raw=False)

adata.obs['cell_types'] = adata.obs[['Plasma_cells', 'T_cells','Mast_cells', 'B_cells', 'Myeloid_cells', 'Dendritic_cells' ]].idxmax(axis=1)

#%%
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

# sc.pl.umap(adata, color=['treatment'], frameon=False, save='_treatment_scVI_HDG.png')
# sc.pl.umap(adata, color=['phase'], frameon=False, save='_cellcycle_scVI_HDG.png')
# sc.pl.umap(adata, color=['cell_types'], frameon=False, save='_celltypes_scVI_HDG.png')
# sc.pl.umap(adata, color=['paper_ID'], frameon=False, save='_patient_scVI_HDG.png')
# sc.pl.umap(adata, color=['tissue'], frameon=False, save='_tissue_scVI_HDG.png')

#%%
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["paper_ID", "leiden"],
    frameon=False,
    ncols=1,
)

sc.pl.embedding(adata, basis="X_mde", color=["cell_types"], frameon=False, ncols=1)

#%%
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scVI_celltypes_embeddings_HDG.h5ad')


# immune seacells integration with scANVI (specifying the label key to preserve)

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
# sc.pl.umap(adata, color=['treatment'], frameon=False, save='_treatment_scANVI_HDG.png')
# sc.pl.umap(adata, color=['phase'], frameon=False, save='_cellcycle_scANVI_HDG.png')
# sc.pl.umap(adataad, color=['tissue-treatment'], frameon=False, save='_tissue-treatment_scANVI_HDG.png')
# sc.pl.umap(adata, color=['paper_ID'], frameon=False, save='_patient_scANVI_HDG.png')
# sc.pl.umap(adata, color=['tissue'], frameon=False, save='_tissue_scANVI_HDG.png')
# sc.pl.umap(adata, color=['cell_types'], frameon=False, save='_celltypes_scANVI_HDG.png')

# %%
adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scANVI_celltypes_embeddings_HDG.h5ad')