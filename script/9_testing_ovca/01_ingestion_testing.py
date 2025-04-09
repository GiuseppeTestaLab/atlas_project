# Reference mapping with scarches
#%%
import os
import scanpy as sc
import pandas as pd
from scipy import sparse
import anndata as an
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
inverted_cell_type_names = {value: key for key, value in cell_type_names.items()}
integrated_data_path = {"fibroblasts": "seacells_hdg_patients_batch_corr_scgen_celltypes_embeddings_HDG.h5ad",
                        "endothelial": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",
                        "cancer": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",}


def adata_by_celltype_tissue(adata):
    adata_dict = {}
    for celltype in adata.obs["01_major_celltypes"].unique():
        adata_celltype = adata[adata.obs["01_major_celltypes"] == celltype]
        celltype_dict = {}
        for tissue in adata_celltype.obs["02_tissue"].unique():
            celltype_dict[tissue] = adata_celltype[adata_celltype.obs["02_tissue"] == tissue]
        adata_dict[celltype] = celltype_dict
    return adata_dict

adata_ref = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
mask = adata_ref.obs['01_major_celltypes'] == "HematopoieticMSK"
ad = adata_ref[~mask]
ref_by = adata_by_celltype_tissue(adata_ref)

#%%
for tissue_dict in ref_by.values():
    for adata in tissue_dict.values():
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.umap(adata)
#%%
ref_by = {key.replace("MSK","").lower():value for key, value in ref_by.items()}

#%%
def rename_query_columns(adata_target):
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
    return adata_target

tissues = ["metastasis", "ascites", "primary"]
query_by = {}
for major_cell_type in major_cell_types:
    celltype_by = {}
    adata_path = f'{rawPath}integration/metacells/{major_cell_type}/{integrated_data_path[major_cell_type]}'
    adata = sc.read_h5ad(adata_path)
    adata = rename_query_columns(adata)
    adata.obs["02_tissue"] = adata.obs["02_tissue"].str.lower()
    for tissue in tissues:
        adata_tissue = adata[adata.obs["02_tissue"] == tissue]
        sc.tl.pca(adata_tissue)
        sc.pp.neighbors(adata_tissue, use_rep='X_pca')
        sc.tl.umap(adata_tissue)
        celltype_by[tissue] = adata_tissue
    query_by[major_cell_type] = celltype_by

#%%
for major_cell_type in major_cell_types:
    for tissue in tissues:
        adata_ref = ref_by[major_cell_type][tissue]
        adata_query = query_by[major_cell_type][tissue]
        adata_query = adata_query[:,adata_ref.var_names]
        sc.tl.ingest(adata_query, adata_ref, obs="07_cell_states")
        ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca_ingestion/'
        if not os.path.exists(ooseDir):
            os.makedirs(ooseDir)
        integrated = an.concat([adata_ref, adata_query], label="origin", keys=["ref", "new"])
        integrated.write_h5ad(f'{ooseDir}{tissue}_ingest.h5ad')

#%%
for major_cell_type in major_cell_types:
    for tissue in tissues:
        adata_ref = ref_by[major_cell_type][tissue]
        adata_query = query_by[major_cell_type][tissue]
        adata_query = adata_query[:,adata_ref.var_names]
        sc.tl.ingest(adata_query, adata_ref, obs="07_cell_states")
        ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca_ingestion/'
        if not os.path.exists(ooseDir):
            os.makedirs(ooseDir)
        integrated = an.concat([adata_ref, adata_query], label="origin", keys=["ref", "new"])
        integrated.write_h5ad(f'{ooseDir}{tissue}_ingest.h5ad')
