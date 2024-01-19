# Adata atlas cancer integration

#%%
import scanpy as sc
import scgen

#%%
ad = sc.read('/group/testa/Project/OvarianAtlas/atlas_cancer_seacells_assignment.h5ad')
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
#sc.pp.normalize_total(ad)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, batch_key='dataset', n_top_genes=2500, inplace=True)

#%%
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment")

#%%
model = scgen.SCGEN(ad)
model.save("/home/marta.sallese/ov_cancer_atlas/Integration/saved_models/model_atlas_cancer_batch_removal_tissue-treatment.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad('/group/testa/Project/OvarianAtlas/Integrated_data/atlas_cancer_batch_corr.h5ad')

#%%
# sc.pp.neighbors(corrected_adata)
# sc.tl.umap(corrected_adata)

