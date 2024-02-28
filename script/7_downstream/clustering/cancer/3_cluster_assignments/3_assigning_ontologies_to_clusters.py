#

## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/cancer/cluster_assignments/'

## Load the data
adata = sc.read(initDir + 'atlas_cell_labelled_cancer_clusters_from_mc.h5ad')
adata
adata.obs

## Remove NA values
# adata.obs.dropna(subset=['cell_types', 'major_celltypes'], inplace=True)
# adata.obs

## Create a new column with the ontology of the cancer cells
adata.obs['cell_ontologies'] = np.nan

col = {"Primary_0":"Cellular_metabolism", 
                              "Primary_1":"Unknown_primary_1", 
                              "Primary_2":"Unknown_primary_2", 
                              "Primary_3":"Immunoreactive_cells", 
                              "Primary_4":"Unknown_primary_4", 
                              "Primary_5":"Cycling_cells", 
                              "Primary_6":"Immunoreactive_cells",
                              "Primary_7":"EMT_cells", 
                              "Primary_8":"Ciliated_cancer_cells", 
                              "Primary_9":"Cellular_metabolism", 
                              "Primary_10":"ECM_shaping_cells", 
                              "Primary_11":"RNA_metabolism", 
                              "Primary_12":"ECM_shaping_cells", 
                              "Ascites_0":"Cellular_metabolism", 
                              "Ascites_1":"Unknown_ascites_1", 
                              "Ascites_2":"Unknown_ascites_2", 
                              "Ascites_3":"Unknown_ascites_3", 
                              "Ascites_4":"Unknown_ascites_4", 
                              "Ascites_5":"Unknown_ascites_5", 
                              "Ascites_6":"Cellular_metabolism",
                              "Ascites_7":"Immunoreactive_cells", 
                              "Ascites_8":"Cycling_cells", 
                              "Ascites_9":"RNA_metabolism", 
                              "Ascites_10":"ECM_shaping_cells", 
                              "Ascites_11":"Cellular_metabolism", 
                                "Ascites_12":"Unknown_ascites_12",
                              "Metastasis_0":"Cellular_metabolism", 
                              "Metastasis_1":"Protein_metabolism", 
                               "Metastasis_2":"Unknown_metastasis_2", 
                              "Metastasis_3":"Cycling_cells", 
                              "Metastasis_4":"Unknown_metastasis_4", 
                              "Metastasis_5":"Unknown_metastasis_5",
                              "Metastasis_6":"Immunoreactive_cells", 
                              "Metastasis_7":"Immunoreactive_cells", 
                              "Metastasis_8":"ECM_shaping_cells",
                              "Metastasis_9":"Organelles_organization",
                              "Metastasis_10":"Unknown_metastasis_10",
                              "Metastasis_11":"Cellular_metabolism",
                              "Metastasis_12":"Cellular_metabolism"}

adata.obs['cell_ontologies'] = adata.obs.cell_type.replace(col)

adata.obs

adata.obs['major_ontologies'] = np.nan

values = []

for index, row in adata.obs.iterrows():
    if row['cell_ontologies'] == 'Unknown_primary_1':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_2':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_4':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_ascites_1':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_2':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_3':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_4':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_5':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_12':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_metastasis_2':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_4':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_5':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_10':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'T_CD4_naive':
        values.append('T_CD4_cells')
    elif row['cell_ontologies'] == 'T_CD4_CXCL13':
        values.append('T_CD4_cells')
    elif row['cell_ontologies'] == 'T_CD4_reg':
        values.append('T_CD4_cells')
    elif row['cell_ontologies'] == 'T_CD8_cytotoxic':
        values.append('T_CD8_cells')
    elif row['cell_ontologies'] == 'T_CD8_ISG':
        values.append('T_CD8_cells')
    elif row['cell_ontologies'] == 'T_CD8_CXCL13':
        values.append('T_CD8_cells')
    elif row['cell_ontologies'] == 'NK_CD56':
        values.append('NK_cells')
    elif row['cell_ontologies'] == 'NK_cytotoxic':
        values.append('NK_cells')
    else:
        values.append(row['cell_ontologies'])

adata.obs['major_ontologies'] = values

set(adata.obs.major_ontologies)

adata.obs.drop(adata.obs[adata.obs['major_ontologies'] == 'nan'].index, inplace = True)

set(adata.obs.major_ontologies)

adata.obs['treatment_ontology'] = adata.obs['treatment'].astype('str') + '_' + adata.obs['major_ontologies'].astype('str')

## Save the data

adata.write_h5ad(initDir + 'atlas_cell_labelled_cancer_clusters_from_mc_ontologies.h5ad')

## Extra code to check for the correct assignment of the ontologies
adata_primary = adata[(adata.obs['tissue'] == 'Primary')]
adata_primary.obs

adata_ascites = adata[(adata.obs['tissue'] == 'Ascites')]
adata_ascites.obs

adata_metastasis = adata[(adata.obs['tissue'] == 'Metastasis')]
adata_metastasis.obs

