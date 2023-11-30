# Adata atlas immune seacells generated in the HDGxpatients_by dispersion space integration

#%%
import scanpy as sc
import scgen

#%%
train = sc.read("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/atlas_immune_seacells_hdg_patients.h5ad")
train.obs['tissue-treatment'] = train.obs['tissue'].astype('str') + '_' + train.obs['treatment'].astype('str')

#%%
# Plasma cells markers
sc.tl.score_genes(train, ['IGKC','IGHG1','CD79A','IGHG2','IGLC2','IGLC3','IGHG3','IGHG4','JCHAIN','MZB1','XBP1'], 
score_name = "Plasma_cells", use_raw=False)

# T cells markers
sc.tl.score_genes(train, ['CD2','CD3D','TRAC','GZMA','NKG7','CD3E','CD3G','CD4','TCF7','CD8A','PRF1','GZMB','CCL5','CCL4','IL32','CD52'], 
score_name = "T_cells", use_raw=False)

# Mast cells markers
sc.tl.score_genes(train, ['KIT','CPA3','CTSG','MS4A2','TPSAB1','TPSB2','HPGD','HPGDS','GATA2'], 
score_name = "Mast_cells", use_raw=False)

# B cells markers
sc.tl.score_genes(train, ['MS4A1', 'CD79A', 'CD19', 'BANK1', 'IGKC', 'IGHM'], score_name = "B_cells", use_raw=False)

# Myeloid cells markers
sc.tl.score_genes(train, 
['CD14','FCER1G','FCGR3A','LYZ','CTSS','CD33','CD68','CD163','ITGAX','ITGAM','CD4','MRC1',
'VSIG4','SPP1','APOE','C1QA','C1QB','C1QC','APOC1','FTL','S100A9','TYROBP','AIF1','CD74','PSAP','CTSB'], 
score_name = "Myeloid_cells", use_raw=False)

# Dendritic cells markers
sc.tl.score_genes(train, 
['IL3RA','IRF7','IRF8','GZMB','CD4','CLEC4C','JCHAIN',
'PTGDS','PLAC8','PLD4','TCF4','BCL11A','GPR183','CCDC50','LILRA4','TSPAN13','CLIC3','MPEG1'], 
score_name = "Dendritic_cells", use_raw=False)

#%%
train.obs
train
train.obs['cell_types'] = train.obs[['Plasma_cells', 'T_cells','Mast_cells', 'B_cells', 'Myeloid_cells', 'Dendritic_cells' ]].idxmax(axis=1)

#%%
#sc.pp.normalize_total(ad)
sc.pp.log1p(train)
sc.pp.highly_variable_genes(train, batch_key='dataset', n_top_genes=2500, inplace=True)

#%%
sc.tl.pca(train, use_highly_variable=True)
sc.pp.neighbors(train, use_rep='X_pca')
sc.tl.umap(train)

#%%
scgen.SCGEN.setup_anndata(train, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(train)
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/model_seacells_immune_batch_removal_tissuetreat.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad('/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/seacells_immune_hdg_patients_batch_corr_scgen_tissuetreat.h5ad')