# Adata atlas fibroblasts integration

#%%
import scanpy as sc
import scgen
import pandas as pd
import sys
import configparser

import os
# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")
CCGenes = config.get("DEFAULT", "CCGenes")

sys.path.insert(1, utilsPath)
from integration import preprocess_scgen

#%%
initDir = rawPath + 'metacells/fibroblasts/'
outDir = rawPath + 'integration/cells/fibroblasts/'
if not os.path.exists(outDir):
    os.makedirs(outDir)

genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_fibroblasts.csv'
#%%
ad = sc.read(initDir + 'seacells_assignment_hdg_patients.h5ad')
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
# inflammatory CAFs (iCAF)
sc.tl.score_genes(ad, ['IL6', 'CXCL12', 'VTN', 'CADM3', 'PLIN2', 'SERPINB2', 'KRT19', 'DES', 'CALB2', 'WT1', 'KRT7'], 
score_name = "iCAF", use_raw=False)

# matrix CAFs (mCAF)
sc.tl.score_genes(ad, ['ACTA2', 'COL11A1', 'MFAP5', 'SFRP2', 'ISLR', 'COL10A1', 'VIM', 'COL3A', 
                       'MMP11', 'PTHLH', 'FGF1', 'WNT7B', 'WNT2', 'TGFB3', 'THRC1', 'POSTN', 'VCAN', 'ZEB1'], 
score_name = "mCAF", use_raw=False)

# vascular cancer-associated fibroblasts CAFs (vCAF)
sc.tl.score_genes(ad, ['COX4I2', 'HIGD1B', 'PTP4A3', 'MCAM', 'PPP1R14A'], 
score_name = "vCAF", use_raw=False)

# STAR gene expressing CAFs (starCAF)
sc.tl.score_genes(ad, ['STAR', 'IGFBP5', 'TSPAN8', 'C7', 'ALDH1A1', 'LGR5'], 
score_name = "starCAF", use_raw=False)

#%%
ad.obs['cell_types'] = ad.obs[['iCAF', 'mCAF','vCAF', 'starCAF']].idxmax(axis=1)
ad.obs

#%%
sc.pp.log1p(ad)
hdg = pd.read_csv(genes, index_col=0)
ad.var['highly_variable'] = hdg.highly_variable
ad.var.highly_variable = ad.var.highly_variable.fillna(False)
ad.raw = ad
ad = ad[:, ad.var.highly_variable]
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="cell_types")

#%%
model = scgen.SCGEN(ad)
model.save(outDir + "saved_models/model_scgen_batch_removal_celltypes_HDG.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'scgen_batch_corr_celltypes_HDG.h5ad')

## Visualization of the corrected adata

#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "integration/cells/fibroblasts/"

#%%
adata = sc.read(outDir + 'scgen_batch_corr_celltypes_HDG.h5ad')

#%%
cell_cycle_genes = [x.strip() for x in open(CCGenes)]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=True)

#%%
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='fibroblasts_cells_HDG_treatm_HDG.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='fibroblasts_cells_HDG_tissue_HDG.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='fibroblasts_cells_HDG_dataset_HDG.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='fibroblasts_cells_HDG_patients_HDG.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='fibroblasts_cells_HDG_cellcycle_HDG.png')
sc.pl.umap(adata, color=["cell_types"], frameon=False, save='fibroblasts_cells_HDG_celltypes_HDG.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='fibroblasts_cells_HDG_tissue-treat_HDG.png')

#%%
adata.write(outDir + 'scgen_batch_corr_celltypes_embeddings_HDG.h5ad')

