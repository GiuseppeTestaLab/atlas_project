# 

#%%
## Import libraries
from mpl_chord_diagram import chord_diagram
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy
import pandas as pd

#%%
## Initialize directories
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'downstream/LR_interactions/cancer_immune/'
figDir = '/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_gc/'

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_metastasis.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

subset_values = ['T_CD8_ISG', 'T_CD8_cells', 'ILC', 'Plasma_cells', 'T_CD4_cells',
                 'NK_cells', 'Myeloid_cells', 'Dendritic_cells', 'Mast_cells',
                 'M1_macrophages', 'B_cells']

mask = matrix_df['source'].isin(subset_values)

source_immune = matrix_df[mask]

#%%
matrix = source_immune.pivot_table(index='target', columns='source', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
## Create color palette
pastel = ['#e5bedd', 
          '#f7d8c2', 
          '#fae087', 
          '#a1edce', 
          '#ffc8c2', 
          '#bdcef4',
          '#b5dcf9',
          '#a9e6e3',
          '#9fd997',
          '#c5cf84',
          '#d2c897']

#%%
grads = (True, False, False, False)               # gradient
gaps  = (0.03, 0, 0.03, 0)                        # gap value
sorts = ("size", "distance", None, "distance")    # sort type
cclrs = (None, None, "slategrey", None)           # chord colors
nrota = (False, False, True, True)                # name rotation
cmaps = ("Set3", "Set3", "magma", "Set3")         # colormap
fclrs = "grey"                                    # fontcolors
fontsize = 7
drctd = (False, False, False, True)               # directed
# pastel = pastel

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "example{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()
# %%

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0)                        # gap value
sorts = ("size", "distance")             # sort type
#cclrs = (None, "slategrey")              # chord colors
nrota = (False, True)                    # name rotation
cmaps = ("Set3", "Set3")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (False, True)                    # directed
pastel = pastel                          # chord_colors

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=pastel, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "example{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()

