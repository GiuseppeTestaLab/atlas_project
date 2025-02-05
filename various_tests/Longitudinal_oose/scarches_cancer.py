# Reference mapping with scarches

#%%
import scanpy as sc
import pandas as pd
import numpy as np
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
from scipy import sparse
from anndata import AnnData
from typing import Optional, Union

## Initialize folders
#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
ooseDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/out_of_sample_extension/cancer/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

## Setting fig parameteres
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/oose/cancer/"

#%%
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

## Loading corrected reference data and integration model
#%%
corrected_reference_adata = sc.read_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_scarches_tissuetreat.h5ad')

#%%
ref_path = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_scarches'
network = sca.models.scgen.load(ref_path, corrected_reference_adata)

## Project query on top of reference

### Query data needs to be preprocessed same way as reference data with same genes
#%%
adata_target = sc.read_h5ad('/group/testa/Project/OvarianAtlas/Longitudinal/Metacells/cancer_seacells_hdg_patients.h5ad')
adata_target.obs['tissue-treatment'] = adata_target.obs['tissue'].astype('str') + '_' + adata_target.obs['treatment'].astype('str')
adata_target.obs['cell_type'] = adata_target.obs.loc[:, 'tissue-treatment']
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
missing_gene = ['FAM19A3', 'FAM129A', 'FAM84A', 'C2orf70', 'C2orf40', 'ALS2CR12',
       'XCR1', 'IL2', 'GAPT', 'TIFAB', 'ANKHD1-EIF4EBP3', 'C6orf48', 'TREML1',
       'GUCA1A', 'RAET1L', 'SEPT6', 'CD40LG', 'C8orf44', 'HRASLS5', 'RARRES3',
       'HRASLS2', 'FAM19A2', 'RNF219', 'KDELC1', 'HIF1A-AS2', 'C15orf53',
       'TPSAB1', 'SEPT1', 'MTSS1L', 'P2RX5-TAX1BP3', 'RNASEK-C17orf49',
       'LGALS9B', 'CCL13', 'RETN', 'FAM129C', 'CEACAM4']
#%%
new_adata = adata_target.copy()
adataX = pd.DataFrame(adata_target.X.todense().T, index = adata_target.var_names, columns = adata_target.obs_names)
adataX = adataX.T
for i in missing_gene:
   adataX[i] = adataX['MT-CYB']*0
adata_new = sc.AnnData(adataX)
adata_new.obs = adata_target.obs
adata_new

adata_new = remove_sparsity(adata_new) # remove sparsity
adata_new.X = adata_new.X.astype('float32')

### This function need pretrained reference model, corrected gene expression from reference data and incorrected query data

#%%
integrated_query = sca.models.scgen.map_query_data(reference_model = network,
                                                   corrected_reference = corrected_reference_adata,
                                                   query = adata_new,
                                                   batch_key = 'paper_ID',
                                                   return_latent=True)

#%%
## Plot the latent space of integrated query and reference

integrated_query.obs['dataset_tt'] = integrated_query.obs['tissue-treatment'].astype('str') + '_' + integrated_query.obs['dataset'].astype('str')

#%%
sc.pp.neighbors(integrated_query, use_rep="latent_corrected")
sc.tl.umap(integrated_query)

### palettes definition
category_of_interest = 'Testa'
unique_categories = integrated_query.obs['dataset'].cat.categories
colors = ['#cccccc' if category != category_of_interest else '#C20078' for category in unique_categories]

#%%
sc.pl.umap(integrated_query, color='dataset', palette=colors, frameon=False, save='oose_long_dataset_latent.png')
# sc.pl.umap(integrated_query, color='paper_ID', frameon=False, save='oose_long_patient_latent.png')
# sc.pl.umap(integrated_query, color='tissue', frameon=False, save='oose_long_tissue_latent.png')
# sc.pl.umap(integrated_query, color='treatment', frameon=False, save='oose_long_treatment_latent.png')
# sc.pl.umap(integrated_query, color='tissue-treatment', frameon=False, save='oose_long_tissue-treatment_latent.png')

#%%
## Plot corrected gene expression space of integrated query and reference
sc.pp.neighbors(integrated_query)
sc.tl.umap(integrated_query)

#%%
integrated_query.write_h5ad(ooseDir + 'integrated_query_longitudinal_seacells_scarches_tissuetreat.h5ad')
