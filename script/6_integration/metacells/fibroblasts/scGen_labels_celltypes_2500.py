# Integration (preserving tissues and treatment) of fibroblasts seacells generated in the HDG space

#%%
import scanpy as sc
import scgen
import pandas as pd
import numpy as np
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from integration import preprocess_scgen_genes

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/fibroblasts/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/fibroblasts/'

#%%
ad = sc.read(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
batch = 'dataset'
genes = 2500
ad = preprocess_scgen_genes(ad, batch, genes)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="cell_types")

#%%
model = scgen.SCGEN(ad)
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/fibroblasts_batch_removal_celltypes_2500.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_celltypes_2500.h5ad')

## Processing of integrated metacells in the same HDG space used to generate metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/integration/metacells/fibroblasts/"

#%%
adata = sc.read(outDir + 'seacells_hdg_patients_batch_corr_scgen_celltypes_2500.h5ad')

#%%
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#%%
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='fibroblasts_seacells_treatm_2500.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='fibroblasts_seacells_tissue_2500.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='fibroblasts_seacells_dataset_2500.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='fibroblasts_seacells_patients_2500.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='fibroblasts_seacells_cellcycle_2500.png')
sc.pl.umap(adata, color=["cell_types"], frameon=False, save='fibroblasts_seacells_celltypes_2500.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='fibroblasts_seacells_tissuetreat_2500.png')

#%%
adata.write(outDir + 'seacells_hdg_patients_batch_corr_scgen_celltypes_embeddings_2500.h5ad')