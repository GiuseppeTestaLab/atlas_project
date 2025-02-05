# Plotting information about cell states

#%%
## Import libraries
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

#%%
## Initialize directories
tissueDir = rawPath + 'downstream/clustering/fibroblasts/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/'

## Plotting parameters
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Avenir'
plt.rcParams['axes.edgecolor']='#333F4B'
plt.rcParams['axes.linewidth']=0.8
plt.rcParams['xtick.color']='#333F4B'
plt.rcParams['ytick.color']='#333F4B'

#%%
## Loading data
primary = sc.read(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

#%%
## Create tables primary
tab1 = pd.DataFrame(primary.obs.cell_states.value_counts())
tab1
tab1.drop(index=tab1.index[2], axis=0, inplace=True)
tab1

tab1.index.astype('object')

my_range=list(range(1,len(tab1.index)+1))
#%%
fig, ax = plt.subplots(figsize=(5,5))

# create for each expense type an horizontal line that starts at x = 0 with the length 
# represented by the specific expense percentage value.
plt.hlines(y=my_range, xmin=0, xmax=tab1['count'], color='#FAC205', alpha=0.2, linewidth=5)

# create for each expense type a dot at the level of the expense percentage value
plt.plot(tab1['count'], my_range, "o", markersize=5, color='#FAC205', alpha=0.6)

# set labels style
ax.set_xlabel('Number of metacells', fontsize=10, color = '#333F4B')
# ax.set_ylabel('Cell States', fontweight='black')

# # set axis
ax.tick_params(axis='both', which='major', labelsize=12)
plt.yticks(my_range, tab1.index)

# add an horizonal label for the y axis 
fig.text(-0.23, 0.96, 'Cell States', fontsize=10, fontweight='black', color = '#333F4B')

# # change the style of the axis spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ax.spines['left'].set_bounds((1, len(my_range)))
# ax.set_xlim(0,25)

ax.spines['left'].set_position(('outward', 8))
ax.spines['bottom'].set_position(('outward', 5))

plt.savefig(figDir + 'hist_primary_cellstates.pdf', format='pdf', bbox_inches='tight')

#%%
## Create tables ascites
tab1 = pd.DataFrame(ascites.obs.cell_states.value_counts())
tab1
tab1.drop(index=tab1.index[0], axis=0, inplace=True)
tab1

tab1.index.astype('object')

my_range=list(range(1,len(tab1.index)+1))
#%%
fig, ax = plt.subplots(figsize=(5,5))

# create for each expense type an horizontal line that starts at x = 0 with the length 
# represented by the specific expense percentage value.
plt.hlines(y=my_range, xmin=0, xmax=tab1['count'], color='#00008B', alpha=0.2, linewidth=5)

# create for each expense type a dot at the level of the expense percentage value
plt.plot(tab1['count'], my_range, "o", markersize=5, color='#00008B', alpha=0.6)

# set labels style
ax.set_xlabel('Number of metacells', fontsize=10, color = '#333F4B')
# ax.set_ylabel('Cell States', fontweight='black')

# # set axis
ax.tick_params(axis='both', which='major', labelsize=12)
plt.yticks(my_range, tab1.index)

# add an horizonal label for the y axis 
fig.text(-0.23, 0.96, 'Cell States', fontsize=10, fontweight='black', color = '#333F4B')

# # change the style of the axis spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ax.spines['left'].set_bounds((1, len(my_range)))
# ax.set_xlim(0,25)

ax.spines['left'].set_position(('outward', 8))
ax.spines['bottom'].set_position(('outward', 5))

plt.savefig(figDir + 'hist_ascites_cellstates.pdf', format='pdf', bbox_inches='tight')

#%%
## Create tables metastasis
tab1 = pd.DataFrame(metastasis.obs.cell_states.value_counts())
tab1
tab1.drop(index=tab1.index[2], axis=0, inplace=True)
tab1

tab1.index.astype('object')

my_range=list(range(1,len(tab1.index)+1))
#%%
fig, ax = plt.subplots(figsize=(5,5))

# create for each expense type an horizontal line that starts at x = 0 with the length 
# represented by the specific expense percentage value.
plt.hlines(y=my_range, xmin=0, xmax=tab1['count'], color='#C20078', alpha=0.2, linewidth=5)

# create for each expense type a dot at the level of the expense percentage value
plt.plot(tab1['count'], my_range, "o", markersize=5, color='#C20078', alpha=0.6)

# set labels style
ax.set_xlabel('Number of metacells', fontsize=10, color = '#333F4B')
# ax.set_ylabel('Cell States', fontweight='black')

# # set axis
ax.tick_params(axis='both', which='major', labelsize=12)
plt.yticks(my_range, tab1.index)

# add an horizonal label for the y axis 
fig.text(-0.23, 0.96, 'Cell States', fontsize=10, fontweight='black', color = '#333F4B')

# # change the style of the axis spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ax.spines['left'].set_bounds((1, len(my_range)))
# ax.set_xlim(0,25)

ax.spines['left'].set_position(('outward', 8))
ax.spines['bottom'].set_position(('outward', 5))

plt.savefig(figDir + 'hist_metastasis_cellstates.pdf', format='pdf', bbox_inches='tight')