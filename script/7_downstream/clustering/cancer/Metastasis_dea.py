# Clustering and cell states annotation of cancer seacells from metastatic tissue

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
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from plotting_bubble import scale_data_5_75, plot_enrich 
from ontologies import annotate_ontolgies

#%%
sc.logging.print_versions()

## inizializing directories
#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/cancer/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/cancer/'

## loading data
#%%
adata = sc.read(initDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings.h5ad')
adata
adata.obs

## Clustering
#%%
adata_mt = adata[(adata.obs['tissue'] == 'Ascites')]
sc.tl.pca(adata_mt, use_highly_variable = True)
sc.pp.neighbors(adata_mt, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata_mt)

sc.pl.umap(adata_mt, color=["treatment"], frameon=False)
sc.pl.umap(adata_mt, color=["phase"], frameon=False)

leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    sc.tl.leiden(adata_mt,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

# for i in leidenTotal:
#    sc.pl.umap(adata_mt, color=i, frameon=False)

## Differential expression analysis
#%%
dedf={}
for lei in leidenTotal:
    dedf[lei]={}
    sc.tl.rank_genes_groups(adata_mt, groupby=lei, method='wilcoxon', key_added = "wilcoxon_"+lei)
    for cl in adata_mt.obs[lei].unique():
        dedf[lei][cl] = sc.get.rank_genes_groups_df(adata_mt, group=cl, key ='wilcoxon_'+lei)

## Assigning gene ontologies to clusters
#%%
directory_root = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/7_downstream/clustering/cancer/metastasis/"
log_file = directory_root + 'metastasis.log'
adata = adata_mt
adata_mt = annotate_ontolgies(adata, directory_root, leidenTotal, dedf, log_file)

logging.shutdown()

## Savings
adata_mt.write_h5ad(outDir + 'adata_metastasis_embeddings.h5ad')