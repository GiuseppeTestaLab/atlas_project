# Reference mapping with scarches
#%%
import os

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from anndata import AnnData
from typing import Optional, Union

## Initialize folders
import seaborn as sns
import configparser
#%%
enrich = pd.read_csv("enrichment/summary_enrichment_ingestion_bg.csv")

#%%

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

major_cell_types = ["fibroblasts", "endothelial", "cancer"]
cell_type_names = {"fibroblasts": "FibroblastsMSK", "immune": "ImmuneMSK", "endothelial": "EndothelialMSK", "cancer": "CancerMSK"}
integrated_data_path = {"fibroblasts": "seacells_hdg_patients_batch_corr_scgen_celltypes_embeddings_HDG.h5ad",
                        "endothelial": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",
                        "cancer": "seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad",}

#%%
adata_ref = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
major_cell_type = "cancer"
tissue = "primary"
#%%
ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca_ingestion/'
path = f'{ooseDir}{tissue}_ingest.h5ad'
adata = sc.read_h5ad(path)
#%%
adata_query = adata[adata.obs.origin == "new"]
adata_ref = adata[adata.obs.origin == "ref"]
#%%
sc.pp.pca(adata_ref)
sc.pp.pca(adata_query)
sc.pl.pca(adata_ref, color=["02_tissue"])
sc.pl.pca_overview(adata_ref, color=["02_tissue"])
sc.pl.pca(adata_query, color=["02_tissue"])
sc.pl.pca_overview(adata_query, color=["02_tissue"])
#%%
sc.pp.neighbors(adata_ref, use_rep='X_pca')
sc.pp.neighbors(adata_query, use_rep='X_pca')
sc.tl.umap(adata_ref)
sc.tl.umap(adata_query)
sc.pl.umap(adata_ref, color=["02_tissue"])
sc.pl.umap(adata_query, color=["02_tissue"])
#%%
CCGenes = config.get("DEFAULT", "CCGenes")
cell_cycle_genes = [x.strip() for x in open(CCGenes)]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_ref.var_names]
sc.tl.score_genes_cell_cycle(adata_ref, s_genes=s_genes, g2m_genes=g2m_genes)
sc.tl.score_genes_cell_cycle(adata_query, s_genes=s_genes, g2m_genes=g2m_genes)
#%%
sc.pl.umap(adata_ref, color=["S_score", "G2M_score", "phase", "07_cell_states"], use_raw=False)
sc.pl.umap(adata_query, color=["S_score", "G2M_score", "phase", "07_cell_states"], use_raw=False)
#%%
adata_ref.obs["07_cell_states"].value_counts()
adata_query.obs["07_cell_states"].value_counts()
# %%
# %%
