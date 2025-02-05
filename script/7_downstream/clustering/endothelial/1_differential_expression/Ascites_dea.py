# Clustering and cell states annotation of endothelial seacells from ascites

## Imports
#%%
import os
import logging
import shutil
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams
from gprofiler import GProfiler
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)
from plotting_bubble import scale_data_5_75, plot_enrich 
from ontologies import annotate_ontolgies

#%%
sc.logging.print_versions()

## inizializing directories
#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'integration/metacells/endothelial/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/endothelial/'

## loading data
#%%
adata = sc.read(initDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad')
adata
adata.obs

adata.obs.drop(columns = ['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 
                          'sample_ID', 'cell_labels_ratio', 
                          'assignment', 'leiden-1.8'], inplace = True)

## Clustering
#%%
adata_as = adata[(adata.obs['tissue'] == 'Ascites')]
sc.tl.pca(adata_as)
sc.pp.neighbors(adata_as, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata_as)

sc.pl.umap(adata_as, color=["treatment"], frameon=False)
sc.pl.umap(adata_as, color=["phase"], frameon=False)

leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    sc.tl.leiden(adata_as,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

# for i in leidenTotal:
#    sc.pl.umap(adata_as, color=i, frameon=False)

## Differential expression analysis
#%%
dedf={}
for lei in leidenTotal:
    if len(adata_as.obs[lei].unique()) == 1:
        continue
    dedf[lei]={}
    sc.tl.rank_genes_groups(adata_as, groupby=lei, method='wilcoxon', key_added = "wilcoxon_"+lei)
    for cl in adata_as.obs[lei].unique():
        dedf[lei][cl] = sc.get.rank_genes_groups_df(adata_as, group=cl, key ='wilcoxon_'+lei)

## Assigning gene ontologies to clusters
#%%
directory_root = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/endothelial/ascites/"
log_file = directory_root + 'ascites.log'
adata = adata_as
adata_as = annotate_ontolgies(adata, directory_root, leidenTotal, dedf, log_file)

logging.shutdown()

## Savings
adata_as.write_h5ad(outDir + 'adata_ascites_embeddings.h5ad')