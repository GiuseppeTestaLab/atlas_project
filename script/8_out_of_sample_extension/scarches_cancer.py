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
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
ooseDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/out_of_sample_extension/cancer/'
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'

## Setting fig parameteres
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/oose/cancer/"

## Loading reference data and creating integration model
#%%
ad = sc.read_h5ad(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')
ad.obs['cell_type'] = ad.obs.loc[:, 'tissue-treatment']
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
condition_key = 'paper_ID'
cell_type_key = 'cell_type'

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

#%%
network.train(n_epochs=epoch, early_stopping_kwargs = early_stopping_kwargs)
#%%
corrected_reference_adata = network.batch_removal(ad, batch_key="paper_ID", cell_label_key="cell_type",return_latent=True)
#%%
corrected_reference_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_scarches_tissuetreat.h5ad')
### Plot latent space
#%%
sc.pp.neighbors(corrected_reference_adata, use_rep="latent_corrected")
sc.tl.umap(corrected_reference_adata)
sc.pl.umap(corrected_reference_adata,
           color=['paper_ID', 'tissue-treatment', 'tissue', 'treatment'],
           frameon=False,
           wspace=0.6,
           )

### Plot corrected data
#%%
sc.pp.neighbors(corrected_reference_adata)
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
missing_gene = 'ZBTB20-AS2'
#%%
new_adata = adata_target.copy()
adataX = pd.DataFrame(adata_target.X.todense().T, index = adata_target.var_names, columns = adata_target.obs_names)
# adata_target
# list(adata_target.var_names) + ['ZBTB20-AS2']
adataX = adataX.T
adataX['ZBTB20-AS2'] = adataX['MT-CYB']*0
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

# integrated_query = sc.read(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat.h5ad')
integrated_query.obs['dataset_tt'] = integrated_query.obs['tissue-treatment'].astype('str') + '_' + integrated_query.obs['dataset'].astype('str')

#%%
sc.pp.neighbors(integrated_query, use_rep="latent_corrected")
sc.tl.umap(integrated_query)

### palettes definition
category_of_interest = 'Zheng'
unique_categories = integrated_query.obs['dataset'].cat.categories
colors = ['#cccccc' if category != category_of_interest else '#C20078' for category in unique_categories]

palette_dataset_tt = {'Ascites_CHT_Ren':'#1f77b433', 
                      'Ascites_NACT_Loret':'#ff7f0e33', 
                      'Ascites_Naive_Loret':'#279e6833',
                      'Ascites_Naive_Vasquez':'#279e6833',
                      'Ascites_Naive_Zheng':'#279e68',
                      'Metastasis_CHT_Geistlinger':'#d6272833',
                      'Metastasis_NACT_Loret':'#aa40fc33',
                      'Metastasis_NACT_Zhang':'#aa40fc33',
                      'Metastasis_Naive_Loret':'#8c564b33',
                      'Metastasis_Naive_Olbrecht':'#8c564b33',
                      'Metastasis_Naive_Qian':'#8c564b33',
                      'Metastasis_Naive_Vasquez':'#8c564b33',
                      'Metastasis_Naive_Zhang':'#8c564b33',
                      'Metastasis_Naive_Zheng':'#8c564b',
                      'Primary_CHT_Ren':'#e377c233',
                      'Primary_NACT_Loret':'#b5bd6133',
                      'Primary_Naive_Loret':'#17becf33',
                      'Primary_Naive_Olbrecht':'#17becf33',
                      'Primary_Naive_Regner':'#17becf33',
                      'Primary_Naive_Vasquez':'#17becf33',
                      'Primary_Naive_Xu':'#17becf33',
                      'Primary_Naive_Zheng':'#17becf'}
#%%
sc.pl.umap(integrated_query, color='dataset', palette=colors, frameon=False, save='oose_dataset_latent.png')
sc.pl.umap(integrated_query, color='paper_ID', frameon=False, save='oose_patient_latent.png')
sc.pl.umap(integrated_query, color='tissue', frameon=False, save='oose_tissue_latent.png')
sc.pl.umap(integrated_query, color='treatment', frameon=False, save='oose_treatment_latent.png')
sc.pl.umap(integrated_query, color='tissue-treatment', frameon=False, save='oose_tissue-treatment_latent.png')
sc.pl.umap(integrated_query, color='dataset_tt', frameon=False, save='oose_tissue-treatment-dataset_latent.png')
#%%
## Plot corrected gene expression space of integrated query and reference
sc.pp.neighbors(integrated_query)
sc.tl.umap(integrated_query)
#%%
sc.pl.umap(integrated_query, color='dataset', palette=colors, frameon=False, save='oose_dataset.png')
sc.pl.umap(integrated_query, color='paper_ID', frameon=False, save='oose_patient.png')
sc.pl.umap(integrated_query, color='tissue', frameon=False, save='oose_tissue.png')
sc.pl.umap(integrated_query, color='treatment', frameon=False, save='oose_treatment.png')
sc.pl.umap(integrated_query, color='tissue-treatment', frameon=False, save='oose_tissue-treatment.png')

#%%
integrated_query.write_h5ad(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat.h5ad')

#%%
# integrated_query = sc.read(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat.h5ad')
sc.pp.neighbors(integrated_query, use_rep="latent_corrected")
sc.tl.umap(integrated_query)
sc.tl.embedding_density(integrated_query, basis='umap', groupby='dataset')
# %%
sc.pl.embedding_density(integrated_query, basis='umap', groupby='dataset', save='density_dataset.pdf')

#%%
integrated_query_sub = integrated_query[integrated_query.obs['tissue-treatment'].isin(['Ascites_Naive', 'Metastasis_Naive', 'Primary_Naive'])]
sc.pl.umap(integrated_query_sub, color='tissue-treatment', frameon=False)

#%%
sc.tl.leiden(integrated_query_sub, resolution=0.03, key_added="leiden-0.03")
sc.pl.umap(integrated_query_sub, color='leiden-0.03', frameon=False)
# %% ## check the accuracy of the mapping: how many cells from zheng 'Ascites_Naive_Zheng' fall into leiden 0.03 cluster 2
integrated_query_sub.obs[integrated_query_sub.obs['dataset_tt'] == 'Ascites_Naive_Zheng'].isin(['leiden-0.03' == 2]).sum()
#%% ## Compute percentage of cells from zheng naive primary assigned to cluster 0

df = pd.crosstab(integrated_query_sub.obs['dataset_tt'], integrated_query_sub.obs['leiden-0.03'])
#%%
primary_zheng = df.loc['Primary_Naive_Zheng']
primary_percentage = (primary_zheng[0] / primary_zheng.sum()) * 100

#%%
ascites_zheng = df.loc['Ascites_Naive_Zheng']
ascites_percentage = (ascites_zheng[2] / ascites_zheng.sum()) * 100

#%%
metastasis_zheng = df.loc['Metastasis_Naive_Zheng']
metastasis_percentage = (metastasis_zheng[1] / metastasis_zheng.sum()) * 100

#%%
# Print the results
print(f"Primary Naive Zheng assigned to cluster 0: {primary_percentage:.2f}%")
print(f"Ascites Naive Zheng assigned to cluster 2: {ascites_percentage:.2f}%")
print(f"Metastasis Naive Zheng assigned to cluster 1: {metastasis_percentage:.2f}%")
# %%
# Plotting
import matplotlib.pyplot as plt
#%%
# Example data based on previous calculations
percentages = {
    'Primary Naive Zheng': 98.12,
    'Ascites Naive Zheng': 79.55,
    'Metastasis Naive Zheng': 98.71
}

# Extract keys and values
labels = list(percentages.keys())
values = list(percentages.values())

# Plotting the data
barWidth = 0.40
plt.figure(figsize=(5, 6))
plt.bar(labels, values, color=['#17becf', '#279e68', '#8c564b'], width=0.2)
plt.grid(False)

# Adding titles and labels
plt.title('Accuracy % of reference mapping')
plt.ylabel('Percentage (%)')
plt.xticks(rotation=45, ha='right', fontsize=10)

# Adding the percentage values on top of the bars
for i, value in enumerate(values):
    plt.text(i, value + 1, f'{value:.2f}%', ha='center')

# Show the plot
plt.show()
# %%
