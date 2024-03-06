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
       "Primary_1":"Cellular_metabolism-extracellular_signaling", 
       "Primary_2":"Organelles_organization-metabolism", 
       "Primary_3":"Unknown_primary_3", 
       "Primary_4":"Unknown_primary_4", 
        "Primary_5":"Cycling_cells", 
        "Primary_6":"Cycling_cells",
        "Primary_7":"Extracellular_signaling", 
        "Primary_8":"Immunoreactive_cells", 
        "Primary_9":"Organelles_organization-cell_cycle", 
        "Primary_10":"Unknown_primary_10", 
        "Primary_11":"Organelles_organization-cell_movement", 
        "Primary_12":"Ciliated_cancer_cells",
        "Primary_13":"Immunoreactive_cells",
        "Primary_14":"Organelles_organization-metabolism",
        "Primary_15":"Unknown_primary_15",
        "Primary_16":"INF_mediated_signaling",
        "Primary_17":"Cellular_metabolism",
        "Primary_18":"Unknown_primary_18",
        "Primary_19":"Unknown_primary_19",
        "Primary_20":"Unknown_primary_20", 
        "Ascites_0":"Organelles_organization-metabolism", 
        "Ascites_1":"Cycling_cells", 
        "Ascites_2":"Unknown_ascites_2", 
        "Ascites_3":"Unknown_ascites_3", 
        "Ascites_4":"Unknown_ascites_4", 
        "Ascites_5":"Unknown_ascites_5", 
        "Ascites_6":"Cellular_metabolism",
        "Ascites_7":"Response_to_extracellular_signals", 
        "Ascites_8":"Unknown_ascites_8", 
        "Ascites_9":"Response_to_stress", 
        "Ascites_10":"Organelles_organization-metabolism",
        "Metastasis_0":"Unknown_metastasis_0", 
        "Metastasis_1":"Cellular_metabolism", 
        "Metastasis_2":"Unknown_metastasis_2", 
        "Metastasis_3":"Unknown_metastasis_3", 
        "Metastasis_4":"Cycling_cells", 
        "Metastasis_5":"Unknown_metastasis_5",
        "Metastasis_6":"Immunoreactive_cells", 
        "Metastasis_7":"Organelles_organization-cell_cycle", 
        "Metastasis_8":"Cycling_cells",
        "Metastasis_9":"Response_to_extracellular_signals",
        "Metastasis_10":"ECM_shaping_cells",
        "Metastasis_11":"Ciliated_cancer_cells",
        "Metastasis_12":"Cellular_metabolism",
        "Metastasis_13":"Extracellular_signaling-immune_cells",
        "Metastasis_14":"Unknown_metastasis_14",
        "Metastasis_15":"RNA_metabolism"}

adata.obs['cell_ontologies'] = adata.obs.cell_type.replace(col)

adata.obs

adata.obs['major_ontologies'] = np.nan

values = []

for index, row in adata.obs.iterrows():
    if row['cell_ontologies'] == 'Unknown_primary_0':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_3':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_4':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_10':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_15':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_18':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_19':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_20':
        values.append('Unknown_primary')
    elif row['cell_ontologies'] == 'Unknown_ascites_2':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_3':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_4':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_5':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_8':
        values.append('Unknown_ascites')
    elif row['cell_ontologies'] == 'Unknown_metastasis_0':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_2':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_3':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_5':
        values.append('Unknown_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_14':
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

