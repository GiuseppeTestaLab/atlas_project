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
clustersDir = rawPath + 'downstream/clustering/endothelial/cluster_assignments/'

## Load the data
clusters = sc.read(clustersDir + 'atlas_endothelial_clusters_from_seacells.h5ad')

clusters
clusters.obs

## Assigning ontologies

clusters.obs['cell_ontologies'] = np.nan

col = {"Primary_0":"Immunoreactive_cells-T_cells", 
        "Primary_1":"Immunoreactive_cells-neutrophils", 
        "Primary_2":"Unknown_primary_2", 
        "Primary_3":"Cellular_metabolism", 
        "Primary_4":"Angiogenesis", 
        "Primary_5":"Unknown_primary_5", 
        "Primary_6":"Unknown_primary_6",
        "Primary_7":"RNA_metabolism", 
        "Primary_8":"Unknown_primary_8", 
        "Primary_9":"Immunoreactive_cells", 
        "Primary_10":"Cellular_metabolism", 
        "Primary_11":"Cycling_cells", 
        "Primary_12":"Unknown_primary_12",
        "Ascites_0":"Unknown_ascites_0", 
        "Ascites_1":"Immunoreactive_cells-T_cells", 
        "Ascites_2":"Immunoreactive_cells-neutrophils", 
        "Ascites_3":"Immunoreactive_cells-neutrophils", 
        "Ascites_4":"RNA_metabolism", 
        "Ascites_5":"Immunoreactive_cells-B_cells", 
        "Ascites_6":"Immunoreactive_cells-B_cells",
        "Ascites_7":"Phagocytic_cells", 
        "Ascites_8":"Angiogenesis",
        "Metastasis_0":"Immunoreactive_cells-T_cells", 
        "Metastasis_1":"Immunoreactive_cells-neutrophils", 
        "Metastasis_2":"Angiogenesis", 
        "Metastasis_3":"Unknown_metastasis_3", 
        "Metastasis_4":"Immunoreactive_cells-B_cells", 
        "Metastasis_5":"Cellular_metabolism",
        "Metastasis_6":"Immunoreactive_cells-neutrophils", 
        "Metastasis_7":"Unknown_metastasis_7", 
        "Metastasis_8":"Immunoreactive_cells-T_cells",
        "Metastasis_9":"Immunoreactive_cells",
        "Metastasis_10":"Immunoreactive_cells-T_cells",
        "Metastasis_11":"Cycling_cells",
        "Metastasis_12":"RNA_metabolism",
        "Metastasis_13":"Immunoreactive_cells-neutrophils",
        "Metastasis_14":"Unknown_metastasis_14",
        "Metastasis_15":"Unknown_metastasis_15"}

clusters.obs['cell_ontologies'] = clusters.obs.cluster_from_seacells.replace(col)

clusters.obs
set(clusters.obs['cell_ontologies'])

clusters.obs['major_ontologies'] = np.nan

values = []

for index, row in clusters.obs.iterrows():
    if row['cell_ontologies'] == 'Unknown_primary_2':
        values.append('Unknown_endothelial_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_5':
        values.append('Unknown_endothelial_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_6':
        values.append('Unknown_endothelial_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_8':
        values.append('Unknown_endothelial_primary')
    elif row['cell_ontologies'] == 'Unknown_primary_12':
        values.append('Unknown_endothelial_primary')
    elif row['cell_ontologies'] == 'Unknown_ascites_0':
        values.append('Unknown_endothelial_ascites')
    elif row['cell_ontologies'] == 'Unknown_metastasis_3':
        values.append('Unknown_endothelial_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_7':
        values.append('Unknown_endothelial_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_14':
        values.append('Unknown_endothelial_metastasis')
    elif row['cell_ontologies'] == 'Unknown_metastasis_15':
        values.append('Unknown_endothelial_metastasis')
    else:
        values.append(row['cell_ontologies'])

clusters.obs['major_ontologies'] = values

clusters.obs
set(clusters.obs['major_ontologies'])

clusters.obs.isna().sum()

## Save the data
clusters.write_h5ad(clustersDir + 'atlas_endothelial_ontologies_from_seacells.h5ad')
