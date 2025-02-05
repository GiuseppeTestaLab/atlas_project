# Adata atlas cancer integration

#%%
import scanpy as sc
import scgen
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
from integration import preprocess_scgen

#%%
initDir = rawPath + 'metacells/cancer/'
outDir = rawPath + 'integration/cells/cancer/'
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'
#%%
ad = sc.read(initDir + 'seacells_assignment_hdg_patients.h5ad')
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

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
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(ad)
model.save(outDir + "saved_models/model_scgen_batch_removal_tissue-treatment_HDG.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'scgen_batch_corr_tissue-treatment_HDG.h5ad')

## Visualization of the corrected adata

#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "integration/cells/cancer/"

#%%
adata = sc.read(outDir + 'scgen_batch_corr_tissue-treatment_HDG.h5ad')

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
sc.pl.umap(adata, color=["treatment"], frameon=False, save='cancer_cells_HDG_treatm_HDG.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='cancer_cells_HDG_tissue_HDG.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='cancer_cells_HDG_dataset_HDG.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='cancer_cells_HDG_patients_HDG.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='cancer_cells_HDG_cellcycle_HDG.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='cancer_cells_HDG_tissue-treat_HDG.png')

#%%
adata.write(outDir + 'cells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad')

