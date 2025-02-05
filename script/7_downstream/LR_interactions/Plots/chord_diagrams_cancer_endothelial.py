# Chord diagrams plotted for cancer-endothelial interactions

#%%
## Import libraries
from mpl_chord_diagram import chord_diagram
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import matplotlib.colors as mcolors
import scipy
import pandas as pd

#%%
## Initialize directories
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'downstream/LR_interactions/cancer_endothelial/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/LR_interactions/cancer_endothelial/'

## Create Color Palette
#%%
# define top and bottom colormaps 
def combine_colormaps(cmap1, cmap2):
    # Get the colors from each colormap
    colors1 = cmap1(np.linspace(0, 1, cmap1.N))
    colors2 = cmap2(np.linspace(0, 1, cmap1.N))

    # Combine the colors
    colors = np.vstack((colors1, colors2))

    # Create and return a new colormap
    return mcolors.ListedColormap(colors)

# Get the tab20b and tab20c colormaps
tab20b = plt.get_cmap('tab20b')
tab20c = plt.get_cmap('tab20c')

# Combine the colormaps
personal_cmap = combine_colormaps(tab20b, tab20c)

# Metastasis
#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_metastasis_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_metastasis{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()


# Ascites

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_ascites_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_ascites{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()

# Primary

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_primary_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_primary{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()

# Naive

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_naive_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_naive{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()

# NACT

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_nact_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_nact{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()

# CHT

#%%
## Load the data
adata = sc.read_h5ad(initDir + 'liana_aggregaterank_cht_labels.h5ad')
adata
adata.obs

#%%
matrix_df = adata.uns['liana_res']
matrix_df

#%%
matrix = matrix_df.pivot_table(index='source', columns='target', values='magnitude_rank', aggfunc='sum', fill_value=0)
matrix

#%%
matrix_values = matrix.values.tolist()

labels = matrix.columns.tolist()
labels

names = labels
flux = matrix_values

#%%
grads = (True, False)                    # gradient
gaps  = (0.03, 0.03)                        # gap value
sorts = ("distance", "distance")             # sort type
cclrs = (None, "slategrey")              # chord colors
nrota = (True, True)                    # name rotation
cmaps = (personal_cmap, "magma")                # colormap
fclrs = "grey"                           # fontcolors
fontsize = 7                             # fontsize
drctd = (True, False)                    # directed

args = (grads, gaps, sorts, cclrs, nrota, cmaps, drctd)

for grd, gap, srt, cc, nr, cm, d in zip(*args):
    chord_diagram(flux, names, gap=gap, use_gradient=grd, sort=srt, directed=d,
                  cmap=cm, chord_colors=None, rotate_names=nr, fontcolor=fclrs, fontsize=fontsize)

    str_grd = "_gradient" if grd else ""

    plt.savefig(figDir + 
        "circos_cht{}_sort-{}{}.png".format(str_grd, srt,
                                                "_directed" if d else ""),
        dpi=600, transparent=True, bbox_inches='tight',
        pad_inches=0.02)

plt.show()
# %%
