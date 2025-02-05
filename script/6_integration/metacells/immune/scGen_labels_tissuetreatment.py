# Integration (preserving tissues and treatment) of immune seacells generated in the HDG space

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

initDir = rawPath + 'metacells/immune/'
outDir = rawPath + 'integration/metacells/immune/'

#%%
ad = sc.read(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
# Plasma cells markers
sc.tl.score_genes(ad, ['IGKC','IGHG1','CD79A','IGHG2','IGLC2','IGLC3','IGHG3','IGHG4','JCHAIN','MZB1','XBP1'], 
score_name = "Plasma_cells", use_raw=False)

# T cells markers
sc.tl.score_genes(ad, ['CD2','CD3D','TRAC','GZMA','NKG7','CD3E','CD3G','CD4','TCF7','CD8A','PRF1','GZMB','CCL5','CCL4','IL32','CD52'], 
score_name = "T_cells", use_raw=False)

# Mast cells markers
sc.tl.score_genes(ad, ['KIT','CPA3','CTSG','MS4A2','TPSAB1','TPSB2','HPGD','HPGDS','GATA2'], 
score_name = "Mast_cells", use_raw=False)

# B cells markers
sc.tl.score_genes(ad, ['MS4A1', 'CD79A', 'CD19', 'BANK1', 'IGKC', 'IGHM'], score_name = "B_cells", use_raw=False)

# Myeloid cells markers
sc.tl.score_genes(ad, 
['CD14','FCER1G','FCGR3A','LYZ','CTSS','CD33','CD68','CD163','ITGAX','ITGAM','CD4','MRC1',
'VSIG4','SPP1','APOE','C1QA','C1QB','C1QC','APOC1','FTL','S100A9','TYROBP','AIF1','CD74','PSAP','CTSB'], 
score_name = "Myeloid_cells", use_raw=False)

# Dendritic cells markers
sc.tl.score_genes(ad, 
['IL3RA','IRF7','IRF8','GZMB','CD4','CLEC4C','JCHAIN',
'PTGDS','PLAC8','PLD4','TCF4','BCL11A','GPR183','CCDC50','LILRA4','TSPAN13','CLIC3','MPEG1'], 
score_name = "Dendritic_cells", use_raw=False)

#%%
ad.obs
ad
ad.obs['cell_types'] = ad.obs[['Plasma_cells', 'T_cells','Mast_cells', 'B_cells', 'Myeloid_cells', 'Dendritic_cells' ]].idxmax(axis=1)

#%%
batch = 'dataset'
genes = 2500
ad = preprocess_scgen(ad, batch, genes)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(ad)
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/immune_batch_removal_tissue-treatment.pt", overwrite=True)

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
sc.settings.figdir = figPath + "integration/metacells/immune/"

#%%
adata = sc.read(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat.h5ad')
adata.obs

#%% 
hvg = pd.read_csv(scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_immune.csv', index_col=0)
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
sc.pl.umap(adata, color=["treatment"], frameon=False, save='immune_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='immune_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='immune_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='immune_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='immune_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='immune_seacells_HDG_tissue-treat.png')
sc.pl.umap(adata, color=["cell_types"], frameon=False, save='immune_seacells_HDG_celltypes.png')

#%%
adata.write(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad')