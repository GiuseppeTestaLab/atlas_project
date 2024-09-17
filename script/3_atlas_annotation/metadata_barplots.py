# Building metadata stacked barplots

#%%
## Imports
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read('../../utils/config.ini')

utilsPath = config.get('DEFAULT', 'utilsPath')
rawPath = config.get('DEFAULT', 'rawPath')
figPath = config.get('DEFAULT', 'figPath')

#%%
## Initialize directories
initDir = rawPath + 'atlas_annotated/'
figDir = figPath + 'atlas_annotated/'

#%%
## Load data
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled.h5ad')

#%%
## Plotting stacked bar plots 

### Cell type
df = pd.crosstab(adata.obs['dataset'], adata.obs['max'])
df

totals = df.sum(axis=1)

# Calculate percentages for each cell type
cell_type1_percent = df['CancerMSK'] / totals * 100
cell_type2_percent = df['EndothelialMSK'] / totals * 100
cell_type3_percent = df['FibroblastsMSK'] / totals * 100
cell_type4_percent = df['HematopoieticMSK'] / totals * 100

# Plotting
barWidth = 0.85
r = range(len(df))

# Set figure dimensions
plt.figure(figsize=(10, 5))

# Create stacked bars
plt.bar(r, cell_type1_percent, color='#F4BBC9', edgecolor='white', width=barWidth, label='Cancer cells')
plt.bar(r, cell_type2_percent, bottom=cell_type1_percent, color='#FC5A50', edgecolor='white', width=barWidth, label='Endothelial cells')
plt.bar(r, cell_type3_percent, bottom=cell_type1_percent+cell_type2_percent, color='#B3E6B5', edgecolor='white', width=barWidth, label='Fibroblasts')
plt.bar(r, cell_type4_percent, bottom=cell_type1_percent+cell_type2_percent+cell_type3_percent, color='#9addfb', edgecolor='white', width=barWidth, label='Immune cells')

# Custom x axis
plt.xticks(r, df.reset_index().dataset)
#plt.xlabel("Dataset")
plt.ylabel("Percentage")

# Show legend
plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

# Save plot
plt.savefig(figDir + 'celltypes_across_datasets.pdf', format='pdf', bbox_inches='tight')

# Show plot
plt.show()

### Tissue
df = pd.crosstab(adata.obs['dataset'], adata.obs['tissue'])
df

totals = df.sum(axis=1)

# Calculate percentages for each cell type
cell_type1_percent = df['Primary'] / totals * 100
cell_type2_percent = df['Ascites'] / totals * 100
cell_type3_percent = df['Metastasis'] / totals * 100

# Plotting
barWidth = 0.85
r = range(len(df))

# Set figure dimensions
plt.figure(figsize=(10, 5))

# Create stacked bars
plt.bar(r, cell_type1_percent, color='#FAC205', edgecolor='white', width=barWidth, label='Primary tissue')
plt.bar(r, cell_type2_percent, bottom=cell_type1_percent, color='#00008B', edgecolor='white', width=barWidth, label='Ascites')
plt.bar(r, cell_type3_percent, bottom=cell_type1_percent+cell_type2_percent, color='#C20078', edgecolor='white', width=barWidth, label='Metastatic tissue')


# Custom x axis
plt.xticks(r, df.reset_index().dataset)
#plt.xlabel("Dataset")
plt.ylabel("Percentage")

# Show legend
plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

# Save plot
plt.savefig(figDir + 'tissues_across_datasets.pdf', format='pdf', bbox_inches='tight')

# Show plot
plt.show()

