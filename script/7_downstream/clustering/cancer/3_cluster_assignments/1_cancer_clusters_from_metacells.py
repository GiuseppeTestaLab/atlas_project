# Assigning clusters from cancer metacells to the original cancer cells

## Objective: Here I will project the clusters from the cancer metacells on the original cancer cells. I will use the cluster level of resolution identified through the analysis of the Sankey plots.

## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import gc

## Initialize directories
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'atlas_annotated/'
metacellsDir = rawPath + 'metacells/cancer/'
tissueDir = rawPath + 'downstream/clustering/cancer/'
figDir = rawPath + 'downstream/clustering/cancer/figures/'
outDir = rawPath + 'downstream/clustering/cancer/cluster_assignments/'

## Load the data
cells = sc.read(initDir + 'atlas_cancer_filt_norm_nolog.h5ad')
adata = sc.read(metacellsDir + 'seacells_assignment_hdg_patients.h5ad')

## Appending the metacells to the cells they belong to
cells.obs.index.equals(adata.obs.index)

cells.obs.drop(columns=['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 
                          'sample_ID', 'cell_labels_ratio', 
                          'assignment', 'leiden-1.8'], inplace = True)

cells.obs = pd.concat([cells.obs, adata.obs.SEACell], axis='columns')

cells.obs['SEACell_patient_tissue'] = cells.obs['SEACell'].astype('str') + '_' + cells.obs['paper_ID'].astype('str') + '_' + cells.obs['tissue'].astype('str')

## Cancer primary: projecting clusters from metacells on the original cells

primary = sc.read(tissueDir + 'adata_primary_embeddings.h5ad')
primary.obs

primary.obs['leiden-0.41'] # 0.xx is the resolution of the clustering chosen for the primary cells

cells.obs['cluster_from_seacells'] = np.nan
cluster_name = 'leiden-0.41'
primary.obs['total_counts_seacell'] = 0

for index, row in primary.obs.iterrows():
    cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Primary'), 'cluster_from_seacells'] = 'Primary_' + row[cluster_name]
    primary.obs.loc[index, 'total_counts_seacell'] = cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Primary'), 'total_counts'].sum()

cells.obs

pd.crosstab(cells.obs.cluster_from_seacells, cells.obs.tissue)

primary.obs

sc.pl.umap(primary, color=['leiden-0.41', 'total_counts_seacell'], frameon=False)

cells.obs.cluster_from_seacells = cells.obs.cluster_from_seacells.astype('str')

del primary
gc.collect()

## Cancer metastasis: projecting clusters from metacells on the original cells  

metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings.h5ad')
metastasis.obs

metastasis.obs['leiden-0.41'] # 0.xx is the resolution of the clustering chosen for the primary cells
cluster_name = 'leiden-0.41'
metastasis.obs['total_counts_seacell'] = 0

for index, row in metastasis.obs.iterrows():
    cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Metastasis'), 'cluster_from_seacells'] = 'Metastasis_' + row[cluster_name]
    metastasis.obs.loc[index, 'total_counts_seacell'] = cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Metastasis'), 'total_counts'].sum()

cells.obs

pd.crosstab(cells.obs.cluster_from_seacells, cells.obs.tissue)

metastasis.obs

sc.pl.umap(metastasis, color=['leiden-0.51', 'total_counts_seacell'], frameon=False)

cells.obs.cluster_from_seacells = cells.obs.cluster_from_seacells.astype('str')

del metastasis
gc.collect()

## Cancer ascites: projecting clusters from metacells on the original cells

ascites = sc.read(tissueDir + 'adata_ascites_embeddings.h5ad')
ascites.obs

ascites.obs['leiden-0.31'] # 0.xx is the resolution of the clustering chosen for the primary cells
cluster_name = 'leiden-0.31'
ascites.obs['total_counts_seacell'] = 0

for index, row in ascites.obs.iterrows():
    cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Ascites'), 'cluster_from_seacells'] = 'Ascites_' + row[cluster_name]
    ascites.obs.loc[index, 'total_counts_seacell'] = cells.obs.loc[(cells.obs['SEACell_patient_tissue'] == index) & (cells.obs['tissue'] == 'Ascites'), 'total_counts'].sum()

cells.obs

pd.crosstab(cells.obs.cluster_from_seacells, cells.obs.tissue)

ascites.obs

sc.pl.umap(ascites, color=['leiden-0.31', 'total_counts_seacell'], frameon=False)

cells.obs.cluster_from_seacells = cells.obs.cluster_from_seacells.astype('str')

del ascites
gc.collect()

## Plotting the clusters from metacells on the original cells

sc.pl.umap(cells, color=['tissue', 'cluster_from_seacells'], frameon=False, save='_cancer_clusters_from_metacells.png')

## Checking for NA values

cells.obs.SEACell.isna().sum()
cells.obs.SEACell_patient_tissue.isna().sum()

## Saving the data

cells.write_h5ad(outDir + 'atlas_cancer_clusters_from_seacells.h5ad')

