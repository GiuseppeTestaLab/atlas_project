# Plotting atlas cells annotated

## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'

## Set figures parameters
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.set_figure_params(transparent=True)
sc.settings.figdir = "/group/testa/Project/OvarianAtlas/atlas_project/plots_def/atlas_annotated/"
palette_tissue = {'Metastasis':'#C20078', 'Ascites':'#00008B', 'Primary':'#FAC205'}
palette_max = {'HematopoieticMSK':'#9addfb', 'CancerMSK':'#F4BBC9', 'FibroblastsMSK':'#B3E6B5', 'EndothelialMSK':'#FC5A50'}

## Load the data
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled_with_ontologies.h5ad')

## Plot UMAPs
sc.pl.umap(adata, color=['max'], palette = palette_max, frameon=False, save='_main_celltypes.png')
sc.pl.umap(adata, color=['ontologies_from_seacells'], frameon=False, save='_ontologies_from_seacells.png')
sc.pl.umap(adata, color=['tissue'], palette = palette_tissue, frameon=False, save='_tissue.png')
sc.pl.umap(adata, color=['treatment'], frameon=False, save='_treatment.png')
sc.pl.umap(adata, color=['anatomical_location'], frameon=False, save='_anatomical_location.png')
sc.pl.umap(adata, color=['paper_ID'], frameon=False, save='_patients.png')
sc.pl.umap(adata, color=['dataset'], frameon=False, save='_dataset.png')