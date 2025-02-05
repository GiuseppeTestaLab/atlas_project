# Assigning ontologies to the clusters identified in the metacells

## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

## Initialize directories
clustersDir = rawPath + 'downstream/clustering/cancer/cluster_assignments/'

## Load the data
clusters = sc.read(clustersDir + 'atlas_cancer_clusters_from_seacells.h5ad')

clusters
clusters.obs

## Assigning ontologies

clusters.obs['cell_ontologies'] = np.nan

col = {"Primary_0":"Unknown_primary_0", 
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

clusters.obs['cell_ontologies'] = clusters.obs.cluster_from_seacells.replace(col)

clusters.obs
set(clusters.obs['cell_ontologies'])

clusters.obs['major_ontologies'] = np.nan

values = []

values = []

for index, row in clusters.obs.iterrows():
    if row['cell_ontologies'] == 'Unknown_primary_0':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_3':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_4':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_10':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_15':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_18':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_19':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_20':
        values.append('Unknown_cancer_primary')
    elif row['cell_ontologies'] == 'Unknown_ascites_2':
        values.append('Unknown_cancer_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_3':
        values.append('Unknown_cancer_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_4':
        values.append('Unknown_cancer_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_5':
        values.append('Unknown_cancer_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_8':
        values.append('Unknown_cancer_ascites')
    elif row['cell_ontologies'] == 'Unknown_metastasis_0':
        values.append('Unknown_cancer_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_2':
        values.append('Unknown_cancer_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_3':
        values.append('Unknown_cancer_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_5':
        values.append('Unknown_cancer_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_14':
        values.append('Unknown_cancer_metastasis')
    else:
        values.append(row['cell_ontologies'])

clusters.obs['major_ontologies'] = values

clusters.obs['major_ontologies'] = values

clusters.obs
set(clusters.obs['major_ontologies'])

clusters.obs.isna().sum()

## Save the data
clusters.write_h5ad(clustersDir + 'atlas_cancer_ontologies_from_seacells.h5ad')