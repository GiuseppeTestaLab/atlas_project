# Assignment of clusters identified in the metacells to the original cells

## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np

## Initialize directories

import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'atlas_annotated/'
fibroDir = rawPath + 'downstream/clustering/fibroblasts/cluster_assignments/'
endoDir = rawPath + 'downstream/clustering/endothelial/cluster_assignments/'
cancerDir = rawPath + 'downstream/clustering/cancer/cluster_assignments/'

sc.settings.figdir = figPath + 'atlas_annotated/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

## Load the data
cells = sc.read(initDir + 'atlas_embeddings_cell_labelled.h5ad')
fibro = sc.read(fibroDir + 'atlas_fibroblasts_ontologies_from_seacells.h5ad')
endo = sc.read(endoDir + 'atlas_endothelial_ontologies_from_seacells.h5ad')
cancer = sc.read(cancerDir + 'atlas_cancer_ontologies_from_seacells.h5ad')

cells
cells.obs

fibro
fibro.obs

endo
endo.obs

cancer
cancer.obs

## Creating a new column called major_celltypes
values = []

for index, row in cells.obs.iterrows():
    if row['cell_types'] == 'T_CD4_naive':
        values.append('T_CD4_cells')
    elif row['cell_types'] == 'T_CD4_CXCL13':
        values.append('T_CD4_cells')
    elif row['cell_types'] == 'T_CD4_reg':
        values.append('T_CD4_cells')
    elif row['cell_types'] == 'T_CD8_cytotoxic':
        values.append('T_CD8_cells')
    elif row['cell_types'] == 'T_CD8_CXCL13':
        values.append('T_CD8_cells')
    elif row['cell_types'] == 'NK_CD56':
        values.append('NK_cells')
    elif row['cell_types'] == 'NK_cytotoxic':
        values.append('NK_cells')
    else:
        values.append(row['cell_types'])

cells.obs['major_celltypes'] = values
cells.obs

## Appending the clusters identified in metacells to the cells they belong to, creating in the cells object a new column called cluster_from_seacells and using "max" column as category to set from whihc anndata object the clusters come from

for data, name in zip([fibro, endo, cancer], ['fibro', 'endo', 'cancer']):
    cells.obs.index.equals(data.obs.index)
    cells.obs = pd.concat([cells.obs, data.obs['major_ontologies'].rename(f'{name}_cluster')], axis='columns')

## I will now merge the three new columns into a single one, called cluster_from_seacells using the "max" column as category to set from whihc anndata object the clusters come from. So if a cell has adata.obs['max'] == 'FibroblastsMSK', then the value in adata.obs['cluster_from_seacells'] will be the value in adata.obs['fibro_cluster'].

cells.obs['ontologies_from_seacells'] = np.nan

for index, row in cells.obs.iterrows():
    if row['max'] == 'FibroblastsMSK':
        cells.obs.loc[index, 'ontologies_from_seacells'] = row['fibro_cluster']
    elif row['max'] == 'EndothelialMSK':
        cells.obs.loc[index, 'ontologies_from_seacells'] = row['endo_cluster']
    elif row['max'] == 'CancerMSK':
        cells.obs.loc[index, 'ontologies_from_seacells'] = row['cancer_cluster']
    elif row['max'] == 'HematopoieticMSK':
        cells.obs.loc[index, 'ontologies_from_seacells'] = row['major_celltypes']

cells.obs

cells.obs.drop(columns=['fibro_cluster', 'endo_cluster', 'cancer_cluster'], inplace=True)

## Plotting the results

sc.pl.umap(cells, color=['max'], frameon=False, save='_main_celltypes.png')
sc.pl.umap(cells, color=['ontologies_from_seacells'], frameon=False, save='_ontologies_from_seacells.png')
sc.pl.umap(cells, color=['tissue'], frameon=False, save='_tissue.png')
sc.pl.umap(cells, color=['treatment'], frameon=False, save='_treatment.png')
sc.pl.umap(cells, color=['anatomical_location'], frameon=False, save='_anatomical_location.png')
sc.pl.umap(cells, color=['paper_ID'], frameon=False, save='_patients.png')
sc.pl.umap(cells, color=['dataset'], frameon=False, save='_dataset.png')

## Save the results
cells.write_h5ad(initDir + 'atlas_embeddings_cell_labelled_with_ontologies.h5ad')





