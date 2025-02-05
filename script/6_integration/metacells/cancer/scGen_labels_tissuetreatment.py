# Integration (preserving tissues and treatment) of cancer seacells generated in the HDG space

#%%
import scanpy as sc
import scgen
import pandas as pd
import numpy as np
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

sys.path.insert(1, utilsPath)
from integration import preprocess_scgen

#%%
initDir = rawPath + 'metacells/cancer/'
outDir = rawPath + 'integration/metacells/cancer/'

#%%
ad = sc.read(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
batch = 'dataset'
genes = 2500
ad = preprocess_scgen(ad, batch, genes)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(ad)
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat.h5ad')

## Processing of integrated metacells in the same HDG space used to generate metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "integration/metacells/cancer/"

#%%
adata = sc.read(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat.h5ad')

#%%
hvg = pd.read_csv(scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
adata.var['highly_variable']=hvg.highly_variable

adata.var.highly_variable = adata.var.highly_variable.fillna(False)

#%%
cell_cycle_genes = [x.strip() for x in open(CCGenes)]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#%%
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='cancer_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='cancer_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='cancer_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='cancer_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='cancer_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='cancer_seacells_HDG_tissue-treat.png')

#%%
adata.write(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad')
