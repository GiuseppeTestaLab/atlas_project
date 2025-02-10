# Clustering and cell states annotation of fibroblasts seacells from primary tissue

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
initDir = rawPath + 'integration/metacells/fibroblasts/'
outDir = rawPath + 'downstream/clustering/fibroblasts/'

## loading data
#%%
adata = sc.read(initDir + 'seacells_hdg_patients_batch_corr_scgen_celltypes_embeddings_HDG.h5ad')
adata
adata.obs

adata.obs.drop(columns = ['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 
                          'sample_ID', 'cell_labels_ratio', 
                          'assignment', 'leiden-1.8'], inplace = True)

## Clustering
#%%
adata_pr = adata[(adata.obs['tissue'] == 'Primary')]
sc.tl.pca(adata_pr)
sc.pp.neighbors(adata_pr, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata_pr)

sc.pl.umap(adata_pr, color=["treatment"], frameon=False)
sc.pl.umap(adata_pr, color=["phase"], frameon=False)

leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    sc.tl.leiden(adata_pr,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

# for i in leidenTotal:
#    sc.pl.umap(adata_pr, color=i, frameon=False)

## Differential expression analysis
#%%
dedf={}
for lei in leidenTotal:
    if len(adata_pr.obs[lei].unique()) == 1:
        continue
    dedf[lei]={}
    sc.tl.rank_genes_groups(adata_pr, groupby=lei, method='wilcoxon', key_added = "wilcoxon_"+lei)
    for cl in adata_pr.obs[lei].unique():
        dedf[lei][cl] = sc.get.rank_genes_groups_df(adata_pr, group=cl, key ='wilcoxon_'+lei)

## Assigning gene ontologies to clusters
#%%
directory_root = rawPath + "downstream/clustering/fibroblasts/primary/"
log_file = directory_root + 'primary.log'
if not os.path.exists(directory_root):
    os.makedirs(directory_root)

adata = adata_pr
adata_pr = annotate_ontolgies(adata, directory_root, leidenTotal, dedf, log_file)

logging.shutdown()

## Savings
if not os.path.exists(outDir):
    os.makedirs(outDir)
adata_pr.write_h5ad(outDir + 'adata_primary_embeddings.h5ad')