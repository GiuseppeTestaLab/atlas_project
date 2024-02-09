#

## Imports
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import kaleido
import matplotlib
import plotly.express as px

sc.logging.print_versions()

## Inizializing folders
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/cancer/'

## Loading data
#%%
adata = sc.read(initDir + 'adata_primary_embeddings.h5ad')

#%%
leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    sc.tl.leiden(adata_pr,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

