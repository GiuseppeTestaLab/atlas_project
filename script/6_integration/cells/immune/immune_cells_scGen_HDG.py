# Adata atlas immune integration

#%%
import scanpy as sc
import scgen
import pandas as pd
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from integration import preprocess_scgen

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/immune/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/cells/immune/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_immune.csv'
#%%
ad = sc.read(initDir + 'seacells_assignment_hdg_patients.h5ad')
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
sc.pp.log1p(ad)
hdg = pd.read_csv(genes, index_col=0)
ad.var['highly_variable'] = hdg.highly_variable
ad.var.highly_variable = ad.var.highly_variable.fillna(False)
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

