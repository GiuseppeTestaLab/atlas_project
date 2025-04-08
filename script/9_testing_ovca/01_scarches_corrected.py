# Reference mapping with scarches
#%%
import os

import scanpy as sc
import pandas as pd
import numpy as np
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
from scipy import sparse
from anndata import AnnData
from typing import Optional, Union

## Initialize folders

import configparser
#%%
# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

major_cell_types = ["fibroblasts", "endothelial", "cancer"]
cell_type_names = {"fibroblasts": "FibroblastsMSK", "immune": "ImmuneMSK", "endothelial": "EndothelialMSK", "cancer": "CancerMSK"}
integrated_data_path = {"fibroblasts": "seacells_hdg_patients_batch_corr_scgen_celltypes_embeddings_HDG.h5ad",
                        "endothelial": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",
                        "cancer": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",}

for major_cell_type in major_cell_types:

    initDir = rawPath + f'metacells_step0/{major_cell_type}/'
    outDir = rawPath + f'integration/metacells/{major_cell_type}_testing_ovca_corrected/'
    ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca_corrected/'
    genes_path = scriptsPath + f'4_hdg/Tables/atlas_hdg_dispersion_patients_{major_cell_type}.csv'

    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    ## Setting fig parameteres
    sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
    sc.settings.figdir = figPath + f"oose/{major_cell_type}_testing_ovca_corrected/"

    ## Loading reference data and creating integration model
    

    ad = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
    major_cell_type_obs = cell_type_names[major_cell_type]
    ad = ad[ad.obs["01_major_celltypes"] == major_cell_type_obs]
    ad.obs['tissue-treatment'] = ad.obs['02_tissue'].astype('str') + '_' + ad.obs['03_treatment'].astype('str')
    ad.obs['cell_type'] = ad.obs.loc[:, 'tissue-treatment']

    
    # sc.pp.log1p(ad)
    
    # ad = ad[~np.any(np.isnan(ad.X), axis=1)]
    # ad = ad.raw.to_adata()
    # hdg = pd.read_csv(genes, index_col=0)
    # ad.var['highly_variable'] = hdg.highly_variable
    # ad.var.highly_variable = ad.var.highly_variable.fillna(False)
    # ad.raw = ad.copy()
    # ad = ad[:, ad.var.highly_variable]
    sc.tl.pca(ad)
    sc.pp.neighbors(ad, use_rep='X_pca')
    sc.tl.umap(ad)

    
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

    
    network.train(n_epochs=epoch, early_stopping_kwargs = early_stopping_kwargs)
    
    corrected_reference_adata = network.batch_removal(ad, batch_key="14_paper_ID", cell_label_key="cell_type",return_latent=True)
    
    corrected_reference_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_scarches_tissuetreat.h5ad')
    ### Plot latent space
    
    sc.pp.neighbors(corrected_reference_adata, use_rep="latent_corrected")
    sc.tl.umap(corrected_reference_adata)
    sc.pl.umap(corrected_reference_adata,
            color=['14_paper_ID', 'tissue-treatment', '02_tissue', '03_treatment', '07_cell_states'],
            frameon=False,
            wspace=0.6,
            )

    ### Plot corrected data
    
    sc.pp.neighbors(corrected_reference_adata)
    sc.tl.umap(corrected_reference_adata)
    sc.pl.umap(corrected_reference_adata,
            color=['14_paper_ID', 'tissue-treatment', '02_tissue', '03_treatment'],
            frameon=False,
            wspace=0.6,
            )
    
    #ref_path = rawPath + f'integration/metacells/{major_cell_type}/seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad'
    ref_path = rawPath + f'integration/metacells/saved_models/{major_cell_type}_batch_removal_tissuetreatment_scarches'
    network.save(ref_path, overwrite=True)

    ## Project query on top of reference

    ### Query data needs to be preprocessed same way as reference data with same genes
    
    #adata_target_path = f"/group/testa/Project/OvarianAtlasTestStep0/raw_data/metacells_step0/{major_cell_type}/" + 'seacells_hdg_patients.h5ad'
    adata_target_path = f'{rawPath}integration/metacells/{major_cell_type}/{integrated_data_path[major_cell_type]}'
    adata_target = sc.read_h5ad(adata_target_path)
    adata_target.obs['tissue-treatment'] = adata_target.obs['tissue'].astype('str') + '_' + adata_target.obs['treatment'].astype('str')
    adata_target.obs['cell_type'] = adata_target.obs.loc[:, 'tissue-treatment']
    
    columns = adata_target.obs.columns

    columns_to_rename = {
        'tissue': '02_tissue',
        'developmental_stage': '12_developmental_stage',
        'treatment': '03_treatment',
        'recurrence': '11_recurrence',
        'tumor_stage': '10_tumor_stage',
        'paper_ID': '14_paper_ID',
        'anatomical_location': '09_anatomical_location',
        'dataset': '13_dataset',
    #    'total_counts': '14_paper_ID'
    }
    valid_columns_to_rename = {k: v for k, v in columns_to_rename.items() if k in columns}

    
    for column in valid_columns_to_rename.keys():
        adata_target.obs[valid_columns_to_rename[column]] = adata_target.obs[column]
        adata_target.obs.pop(column)
    
    #sc.pp.log1p(adata_target)
    
    hdg = pd.read_csv(genes_path, index_col=0)
    adata_target.var['highly_variable'] = hdg.highly_variable
    adata_target.var.highly_variable = adata_target.var.highly_variable.fillna(False)
    adata_target.raw = adata_target
    
    adata_target = adata_target[:, adata_target.var.highly_variable]
    
    sc.tl.pca(adata_target, use_highly_variable=True)
    sc.pp.neighbors(adata_target, use_rep='X_pca')
    sc.tl.umap(adata_target)
    
    #genes = pd.read_csv(scriptsPath + f'4_hdg/Tables/atlas_hdg_dispersion_patients_{major_cell_type}.csv', index_col=0)
    # missing_gene = genes[~genes.index.isin(adata_target.var_names)].index
    # missing_gene
    # missing_gene = ['ZBTB20-AS2', 'OTUD6A']
    
    new_adata = adata_target.copy()
    adataX = pd.DataFrame(adata_target.X.T, index = adata_target.var_names, columns = adata_target.obs_names)
    adataX = adataX.T
    # adataX['ZBTB20-AS2'] = adataX['MT-CO1']*0
    # adataX['OTUD6A'] = adataX['MT-CO1']*0
    adata_new = sc.AnnData(adataX)
    adata_new.obs = adata_target.obs

    adata_new = remove_sparsity(adata_new) # remove sparsity
    adata_new.X = adata_new.X.astype('float32')
    ### This function need pretrained reference model, corrected gene expression from reference data and incorrected query data

    
    adata_new.obs_names = adata_new.obs_names.str.replace("SEACell", "new_SEACell")
    integrated_query = sca.models.scgen.map_query_data(reference_model = network,
                                                    corrected_reference = corrected_reference_adata,
                                                    query = adata_new,
                                                    batch_key = '14_paper_ID',
                                                    return_latent=True)

    raw_concat = adata_target.raw.to_adata().concatenate(corrected_reference_adata.raw.to_adata())
    integrated_query.raw = raw_concat

    ## Plot the latent space of integrated query and reference
    sc.pp.neighbors(integrated_query, use_rep="latent_corrected")
    sc.tl.umap(integrated_query)

    
    sc.pl.umap(integrated_query, color='14_paper_ID', frameon=False, save=f'oose_patient_latent_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='02_tissue', frameon=False, save=f'oose_tissue_latent_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='03_treatment', frameon=False, save=f'oose_treatment_latent_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='tissue-treatment', frameon=False, save=f'oose_tissue-treatment_latent_{major_cell_type}.png')
    
    ## Plot corrected gene expression space of integrated query and reference
    sc.pp.neighbors(integrated_query)
    sc.tl.umap(integrated_query)
    
    sc.pl.umap(integrated_query, color='14_paper_ID', frameon=False, save=f'oose_patient_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='02_tissue', frameon=False, save=f'oose_tissue_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='03_treatment', frameon=False, save=f'oose_treatment_{major_cell_type}.png')
    sc.pl.umap(integrated_query, color='tissue-treatment', frameon=False, save=f'oose_tissue-treatment_{major_cell_type}.png')

    
    if not os.path.exists(ooseDir):
        os.makedirs(ooseDir)

    integrated_query.write_h5ad(ooseDir + f'integrated_query_seacells_scarches_tissuetreat_{major_cell_type}_raw_corrected.h5ad')

    # %%
