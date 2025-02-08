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
clustersDir = rawPath + 'downstream/clustering/fibroblasts/cluster_assignments/'

## Load the data
clusters = sc.read(clustersDir + 'atlas_fibroblasts_clusters_from_seacells.h5ad')

clusters
clusters.obs

## Assigning ontologies

clusters.obs['cell_ontologies'] = np.nan

col = {"Primary_0":"Cellular_metabolism", 
       "Primary_1":"ECM_shaping_cells", 
       "Primary_2":"Unknown_primary_2", 
       "Primary_3":"Smooth_muscle_cells_development", 
        "Primary_4":"Immunoreactive_cells", 
        "Primary_5":"ECM_shaping_cells", 
        "Primary_6":"Cellular_metabolism",
        "Primary_7":"Protein_metabolism-cell_death", 
        "Primary_8":"Unknown_primary_8", 
        "Primary_9":"Epithelium_development", 
        "Primary_10":"Cycling_cells", 
        "Primary_11":"Angiogenesis", 
        "Primary_12":"Immunoreactive_cells",
        "Primary_13":"Cellular_metabolism",
        "Primary_14":"Response_to_stress-ROS",
        "Ascites_0":"Extracellular_tissue_development", 
        "Ascites_1":"Unknown_ascites_1", 
        "Ascites_2":"Cellular_metabolism", 
        "Ascites_3":"Cellular_metabolism-ECM", 
        "Ascites_4":"Unknown_ascites_4", 
        "Ascites_5":"Unknown_ascites_5", 
        "Ascites_6":"Unknown_ascites_6",
        "Ascites_7":"RNA_metabolism", 
        "Ascites_8":"Unknown_ascites_8",
        "Ascites_9":"Unknown_ascites_9",
        "Ascites_10":"Cycling_cells",
        "Ascites_11":"ECM_shaping_cells",
        "Ascites_12":"Cellular_metabolism",
        "Ascites_13":"Unknown_ascites_13",
        "Ascites_14":"Angiogenesis",
        "Ascites_15":"Unknown_ascites_15",
        "Metastasis_0":"Collagen_degradation", 
        "Metastasis_1":"Protein_catabolism", 
        "Metastasis_2":"Cellular_metabolism", 
        "Metastasis_3":"Unknown_metastasis_3", 
        "Metastasis_4":"Vascular_processes_regulation", 
        "Metastasis_5":"Immunoreactive_cells-T_cells",
        "Metastasis_6":"Cellular_metabolism", 
        "Metastasis_7":"ECM_shaping_cells", 
        "Metastasis_8":"Cellular_metabolism",
        "Metastasis_9":"ECM_shaping_cells",
        "Metastasis_10":"Angiogenesis",
        "Metastasis_11":"Unknown_metastasis_11",
        "Metastasis_12":"Epithelium_development-cell_division",
        "Metastasis_13":"Smooth_muscle_cells_development",
        "Metastasis_14":"Cycling_cells",
        "Metastasis_15":"Unknown_metastasis_15",
        "Metastasis_16":"Immunoreactive_cells-T_cells"}

clusters.obs['cell_ontologies'] = clusters.obs.cluster_from_seacells.replace(col)

clusters.obs
set(clusters.obs['cell_ontologies'])

clusters.obs['major_ontologies'] = np.nan

values = []

for index, row in clusters.obs.iterrows():
    if row['cell_ontologies'] == 'Unknown_primary_2':
        values.append('Unknown_fibroblasts_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_8':
        values.append('Unknown_fibroblasts_primary')
    elif row['cell_ontologies'] == 'Unknown_ascites_1':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_4':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_5':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_6':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_8':
        values.append('Unknown_fibroblasts_ascites') 
    elif row['cell_ontologies'] == 'Unknown_ascites_9':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_13':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_ascites_15':
        values.append('Unknown_fibroblasts_ascites')
    elif row['cell_ontologies'] == 'Unknown_metastasis_3':
        values.append('Unknown_fibroblasts_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_11':
        values.append('Unknown_fibroblasts_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_15':
        values.append('Unknown_fibroblasts_metastasis')
    else:
        values.append(row['cell_ontologies'])

clusters.obs['major_ontologies'] = values

clusters.obs
set(clusters.obs['major_ontologies'])

clusters.obs.isna().sum()

## Save the data
clusters.write_h5ad(clustersDir + 'atlas_fibroblasts_ontologies_from_seacells.h5ad')
