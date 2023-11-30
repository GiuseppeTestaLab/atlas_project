# Adata atlas cancer seacells generated in the HDGxpatients_by dispersion space integration with prior HDG space selection

#%%
import scanpy as sc
import scgen
import pandas as pd

#%%
ad = sc.read("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/atlas_cancer_seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
#sc.pp.normalize_total(ad)
sc.pp.log1p(ad)

hvg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/hdg/cancer/atlas_cancer_hdg_dispersion_patients.csv',  index_col=0)
hvg[hvg.highly_variable]
ad.var['highly_variable']=hvg.highly_variable
ad.var.highly_variable = ad.var.highly_variable.fillna(False)

ad.raw = ad
ad = ad[:, ad.var.highly_variable]
ad.var_names #9288 genes
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
model.save("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/model_seacells_cancer_batch_removal.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad('/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/seacells_cancer_hdg_patients_batch_corr_scgen_tissuetreat_HDG.h5ad')