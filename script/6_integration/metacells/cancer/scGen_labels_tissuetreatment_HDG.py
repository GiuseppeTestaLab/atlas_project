# Integration (preserving tissues and treatment) of cancer seacells generated in the HDG space with prior HDG space selection

#%%
import scanpy as sc
import scgen
import pandas as pd
import numpy as np

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

#%%
ad = sc.read(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
sc.pp.log1p(ad)
hdg = pd.read_csv(genes, index_col=0)
ad.var['highly_variable'] = hdg.highly_variable
ad.var.highly_variable = ad.var.highly_variable.fillna(False)
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(ad)
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_HDG.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_HDG.h5ad')

## Processing of integrated metacells in the same HDG space used to generate metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/integration/metacells/cancer/"

#%%
adata = sc.read(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_HDG.h5ad')

#%%
hvg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
adata.var['highly_variable']=hvg.highly_variable

adata.var.highly_variable = adata.var.highly_variable.fillna(False)

#%%
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#%%
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='cancer_seacells_HDG_treatm_HDG.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='cancer_seacells_HDG_tissue_HDG.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='cancer_seacells_HDG_dataset_HDG.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='cancer_seacells_HDG_patients_HDG.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='cancer_seacells_HDG_cellcycle_HDG.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='cancer_seacells_HDG_tissue-treat_HDG.png')

#%%
adata.write(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad')