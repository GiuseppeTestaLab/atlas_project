# Reference mapping with scarches

#%%
import scanpy as sc
import pandas as pd
import numpy as np
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
from scipy import sparse
# from scvi import REGISTRY_KEYS
# from scvi.data import AnnDataManager
# from scvi.data.fields import CategoricalObsField, LayerField


## Initialize folders
#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

## Setting fig parameteres
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/integration/metacells/cancer/"

## Loading reference data and creating integration model
#%%
ad = sc.read_h5ad(initDir + "seacells_hdg_patients.h5ad")
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

from scipy import sparse
import scanpy as sc

def remove_sparsity(adata):
    """
    If `adata.X` is a sparse matrix, this will convert it to a normal (dense) matrix.

    Parameters:
    adata : AnnData
        Annotated dataset.

    Returns:
    AnnData
        Annotated dataset with dense matrix if it was originally sparse.
    """
    if sparse.issparse(adata.X):
        # Use toarray() method to convert sparse matrix to dense
        new_adata = sc.AnnData(X=adata.X.toarray(), obs=adata.obs.copy(deep=True), var=adata.var.copy(deep=True))
        return new_adata
    return adata

ad = remove_sparsity(ad) # remove sparsity
ad.X = ad.X.astype('float32')

#%%
epoch = 100

early_stopping_kwargs = {
    "early_stopping_metric": "val_loss",
    "patience": 25,
    "threshold": 0,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

network = sca.models.scgen(adata = ad, hidden_layer_sizes=[800,100])

# %%
network.train(n_epochs=epoch, early_stopping_kwargs = early_stopping_kwargs)
#%%
corrected_reference_adata = network.batch_removal(ad, batch_key="paper_ID", cell_label_key="tissue-treatment",return_latent=True)
#%%
corrected_reference_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_scarches_tissuetreat.h5ad')
#%%
sc.pp.neighbors(corrected_reference_adata, use_rep="latent_corrected")
sc.tl.umap(corrected_reference_adata)
sc.pl.umap(corrected_reference_adata,
           color=['paper_ID', 'tissue-treatment', 'tissue', 'treatment'],
           frameon=False,
           wspace=0.6,
           )
#%%
ref_path = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_scarches'
network.save(ref_path, overwrite=True)

## Project query on top of reference

### Query data needs to be preprocessed same way as reference data with same genes
#%%
adata_target = sc.read_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Metacells/cancer_seacells_hdg_patients.h5ad')
adata_target.obs['tissue-treatment'] = adata_target.obs['tissue'].astype('str') + '_' + adata_target.obs['treatment'].astype('str')
#%%
sc.pp.log1p(adata_target)
hdg = pd.read_csv(genes, index_col=0)
adata_target.var['highly_variable'] = hdg.highly_variable
adata_target.var.highly_variable = adata_target.var.highly_variable.fillna(False)
adata_target.raw = adata_target
adata_target = adata_target[:, adata_target.var.highly_variable]
sc.tl.pca(adata_target, use_highly_variable=True)
sc.pp.neighbors(adata_target, use_rep='X_pca')
sc.tl.umap(adata_target)
#%%
genes = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
missing_gene = genes[~genes.index.isin(adata_target.var_names)].index
missing_gene
missing_gene = 'ZBTB20-AS2'
#%%
new_adata = adata_target.copy()
adataX=pd.DataFrame(adata_target.X.todense().T,index=adata_target.var_names,columns=adata_target.obs_names)
# adata_target
# list(adata_target.var_names) + ['ZBTB20-AS2']
adataX=adataX.T
adataX['ZBTB20-AS2']=adataX['MT-CYB']*0
adata_new=sc.AnnData(adataX)
adata_new.obs=adata_target.obs
adata_new

### This function need pretrained reference model, corrected gene expression from reference data and incorrected query data

#%%
integrated_query = sca.models.scgen.map_query_data(reference_model = network,
                                                   corrected_reference = corrected_reference_adata,
                                                   query = adata_new,
                                                   batch_key = 'paper_ID',
                                                   return_latent=True)

## Plot the latent space of integrated query and reference
sc.pp.neighbors(integrated_query, use_rep="latent_corrected")
sc.tl.umap(integrated_query)
sc.pl.umap(
    integrated_query,
    color=['paper_ID', 'tissue-treatment', 'tissue', 'treatment', 'dataset'],
    frameon=False,
    wspace=0.6)

## Plot corrected gene expression space of integrated query and reference
sc.pp.neighbors(integrated_query)
sc.tl.umap(integrated_query)
sc.pl.umap(
    integrated_query,
    color=['paper_ID', 'tissue-treatment', 'tissue', 'treatment', 'dataset'],
    frameon=False,
    wspace=0.6)