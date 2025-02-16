
#library imports
#%%
import scanpy as sc
import pandas as pd
import sys
import os
import zipfile
import tempfile
import glob
import argparse
import numpy as np
import anndata as ad
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")
scriptsPath = config.get("DEFAULT", "scriptsPath")
rawPath = config.get("DEFAULT", "rawPath")

# inputs
parameters = pd.read_csv(
    scriptsPath + "8_out_of_sample_extension/01_Zheng2023/preprocess_params.tsv", sep=";"
)
#inputs

init_dir = os.path.join(rawPath, parameters.init_dir[0])
out_dir = os.path.join(rawPath, parameters.out_dir[0])
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum = float(parameters.target_sum[0])
final_dir = os.path.join(rawPath, parameters.final_dir[0])

#create anndata
#%%
adata = sc.read_mtx(init_dir + "Counts/counts.mtx")

adata = ad.AnnData(adata)
adata = adata.T

var = pd.read_csv(init_dir + "Counts/var_names.csv", index_col=0)
var.columns
var.columns[0]

#%%
var.rename(columns= {var.columns[0] : 'genes'}, inplace=True)
adata.var.index = var.iloc[:,0]
adata.var.index= adata.var_names

#Paper:  PMID: 
#Counts file: mtx

# Creating metadata
#%%
meta = pd.read_csv(init_dir + "Counts/metadata.csv", index_col=0)
adata.obs = meta 
adata.obs.drop(columns=['nCount_RNA', 'nFeature_RNA', 'Cellname', 'percent.mt', 'percent.ribo', 'percent.HSP',
       'S.Score', 'G2M.Score', 'Annotation', 'maintypes_2',
       'maintypes_3', 'UMAP_1', 'UMAP_2'], inplace=True)

#%%
columns = {}

columns['tissue'] = {'HGSOC1_BC':'Blood', 'HGSOC1_PT':'Primary', 'HGSOC1_MT':'Metastasis', 'HGSOC1_AS':'Ascites', 'HGSOC3_BC':'Blood',
       'HGSOC3_PT':'Primary', 'HGSOC3_MT':'Metastasis', 'HGSOC3_LN':'Lymph node', 'HGSOC3_AS':'Ascites', 'HGSOC2_BC':'Blood',
       'HGSOC2_PT':'Primary', 'HGSOC2_LN':'Lymph node', 'HGSOC2_AS':'Ascites', 'HGSOC7_PT':'Primary', 'UOC1_PT':'Primary',
       'HGSOC6_AS':'Ascites', 'HGSOC6_LN':'Lymph node', 'HGSOC6_MT':'Metastasis', 'HGSOC6_BC':'Blood', 'HGSOC6_PT':'Primary',
       'OCCC1_AS':'Ascites', 'HGSOC5_PT':'Primary', 'HGSOC5_AS':'Ascites', 'HGSOC4_LN':'Lymph node', 'HGSOC4_MT':'Metastasis',
       'HGSOC4_BC':'Blood', 'HGSOC4_PT':'Primary', 'HGSOC8_PT':'Primary', 'HGSOC8_AS':'Ascites', 'HGSOC9_PT':'Primary',
       'HGSOC9_AS':'Ascites', 'ECO1_BC':'Blood', 'ECO1_PT':'Primary', 'ECO1_MT':'Metastasis', 'ECO1_LN':'Lymph node', 'ECO1_AS':'Ascites',
       'HGSOC10_AS':'Ascites', 'HGSOC10_PT':'Primary', 'C1_PT':'Primary'}

columns['developmental_stage'] = {'HGSOC1_BC':'Advanced stage', 'HGSOC1_PT':'Advanced stage', 'HGSOC1_MT':'Advanced stage', 'HGSOC1_AS':'Advanced stage', 'HGSOC3_BC':'Advanced stage',
       'HGSOC3_PT':'Advanced stage', 'HGSOC3_MT':'Advanced stage', 'HGSOC3_LN':'Advanced stage', 'HGSOC3_AS':'Advanced stage', 'HGSOC2_BC':'Advanced stage',
       'HGSOC2_PT':'Advanced stage', 'HGSOC2_LN':'Advanced stage', 'HGSOC2_AS':'Advanced stage', 'HGSOC7_PT':'Advanced stage', 'UOC1_PT':'Advanced stage',
       'HGSOC6_AS':'Advanced stage', 'HGSOC6_LN':'Advanced stage', 'HGSOC6_MT':'Advanced stage', 'HGSOC6_BC':'Advanced stage', 'HGSOC6_PT':'Advanced stage',
       'OCCC1_AS':'Advanced stage', 'HGSOC5_PT':'Advanced stage', 'HGSOC5_AS':'Advanced stage', 'HGSOC4_LN':'Advanced stage', 'HGSOC4_MT':'Advanced stage',
       'HGSOC4_BC':'Advanced stage', 'HGSOC4_PT':'Advanced stage', 'HGSOC8_PT':'Advanced stage', 'HGSOC8_AS':'Advanced stage', 'HGSOC9_PT':'Advanced stage',
       'HGSOC9_AS':'Advanced stage', 'ECO1_BC':'Advanced stage', 'ECO1_PT':'Advanced stage', 'ECO1_MT':'Advanced stage', 'ECO1_LN':'Advanced stage', 'ECO1_AS':'Advanced stage',
       'HGSOC10_AS':'Advanced stage', 'HGSOC10_PT':'Advanced stage', 'C1_PT':'Advanced stage'}

columns['patient_id'] = {'HGSOC1_BC':'P1', 'HGSOC1_PT':'P1', 'HGSOC1_MT':'P1', 'HGSOC1_AS':'P1', 'HGSOC3_BC':'P3',
       'HGSOC3_PT':'P3', 'HGSOC3_MT':'P3', 'HGSOC3_LN':'P3', 'HGSOC3_AS':'P3', 'HGSOC2_BC':'P2',
       'HGSOC2_PT':'P2', 'HGSOC2_LN':'P2', 'HGSOC2_AS':'P2', 'HGSOC7_PT':'P7', 'UOC1_PT':'Undifferentiated',
       'HGSOC6_AS':'P6', 'HGSOC6_LN':'P6', 'HGSOC6_MT':'P6', 'HGSOC6_BC':'P6', 'HGSOC6_PT':'P6',
       'OCCC1_AS':'Clear cell carcinoma', 'HGSOC5_PT':'P5', 'HGSOC5_AS':'P5', 'HGSOC4_LN':'P4', 'HGSOC4_MT':'P4',
       'HGSOC4_BC':'P4', 'HGSOC4_PT':'P4', 'HGSOC8_PT':'P8', 'HGSOC8_AS':'P8', 'HGSOC9_PT':'P9',
       'HGSOC9_AS':'P9', 'ECO1_BC':'Endometroid_carcinoma', 'ECO1_PT':'Endometroid_carcinoma', 'ECO1_MT':'Endometroid_carcinoma', 'ECO1_LN':'Endometroid_carcinoma', 'ECO1_AS':'Endometroid_carcinoma',
       'HGSOC10_AS':'P10', 'HGSOC10_PT':'P10', 'C1_PT':'Carcinosarcoma'}

columns['treatment'] = {'HGSOC1_BC':'Naive', 'HGSOC1_PT':'Naive', 'HGSOC1_MT':'Naive', 'HGSOC1_AS':'Naive', 'HGSOC3_BC':'Naive',
       'HGSOC3_PT':'Naive', 'HGSOC3_MT':'Naive', 'HGSOC3_LN':'Naive', 'HGSOC3_AS':'Naive', 'HGSOC2_BC':'Naive',
       'HGSOC2_PT':'Naive', 'HGSOC2_LN':'Naive', 'HGSOC2_AS':'Naive', 'HGSOC7_PT':'Naive', 'UOC1_PT':'NACT',
       'HGSOC6_AS':'Naive', 'HGSOC6_LN':'Naive', 'HGSOC6_MT':'Naive', 'HGSOC6_BC':'Naive', 'HGSOC6_PT':'Naive',
       'OCCC1_AS':'Naive', 'HGSOC5_PT':'Naive', 'HGSOC5_AS':'Naive', 'HGSOC4_LN':'Naive', 'HGSOC4_MT':'Naive',
       'HGSOC4_BC':'Naive', 'HGSOC4_PT':'Naive', 'HGSOC8_PT':'Naive', 'HGSOC8_AS':'Naive', 'HGSOC9_PT':'Naive',
       'HGSOC9_AS':'Naive', 'ECO1_BC':'Naive', 'ECO1_PT':'Naive', 'ECO1_MT':'Naive', 'ECO1_LN':'Naive', 'ECO1_AS':'Naive',
       'HGSOC10_AS':'Naive', 'HGSOC10_PT':'Naive', 'C1_PT':'Naive'}

columns['recurrence'] = {'HGSOC1_BC':'Sensitive', 'HGSOC1_PT':'Sensitive', 'HGSOC1_MT':'Sensitive', 'HGSOC1_AS':'Sensitive', 'HGSOC3_BC':'Recurrence',
       'HGSOC3_PT':'Recurrence', 'HGSOC3_MT':'Recurrence', 'HGSOC3_LN':'Recurrence', 'HGSOC3_AS':'Recurrence', 'HGSOC2_BC':'Sensitive',
       'HGSOC2_PT':'Sensitive', 'HGSOC2_LN':'Sensitive', 'HGSOC2_AS':'Sensitive', 'HGSOC7_PT':'Recurrence', 'UOC1_PT':'Unknown',
       'HGSOC6_AS':'Recurrence', 'HGSOC6_LN':'Recurrence', 'HGSOC6_MT':'Recurrence', 'HGSOC6_BC':'Recurrence', 'HGSOC6_PT':'Recurrence',
       'OCCC1_AS':'Sensitive', 'HGSOC5_PT':'Sensitive', 'HGSOC5_AS':'Sensitive', 'HGSOC4_LN':'Sensitive', 'HGSOC4_MT':'Sensitive',
       'HGSOC4_BC':'Sensitive', 'HGSOC4_PT':'Sensitive', 'HGSOC8_PT':'Sensitive', 'HGSOC8_AS':'Sensitive', 'HGSOC9_PT':'Sensitive',
       'HGSOC9_AS':'Sensitive', 'ECO1_BC':'Recurrence', 'ECO1_PT':'Recurrence', 'ECO1_MT':'Recurrence', 'ECO1_LN':'Recurrence', 'ECO1_AS':'Recurrence',
       'HGSOC10_AS':'Sensitive', 'HGSOC10_PT':'Sensitive', 'C1_PT':'Sensitive'}

columns['tumor_stage'] = {'HGSOC1_BC':'IIIC', 'HGSOC1_PT':'IIIC', 'HGSOC1_MT':'IIIC', 'HGSOC1_AS':'IIIC', 'HGSOC3_BC':'IIIC',
       'HGSOC3_PT':'IIIC', 'HGSOC3_MT':'IIIC', 'HGSOC3_LN':'IIIC', 'HGSOC3_AS':'IIIC', 'HGSOC2_BC':'IIB',
       'HGSOC2_PT':'IIB', 'HGSOC2_LN':'IIB', 'HGSOC2_AS':'IIB', 'HGSOC7_PT':'IIIC', 'UOC1_PT':'IIIC',
       'HGSOC6_AS':'IIIC', 'HGSOC6_LN':'IIIC', 'HGSOC6_MT':'IIIC', 'HGSOC6_BC':'IIIC', 'HGSOC6_PT':'IIIC',
       'OCCC1_AS':'IIIC', 'HGSOC5_PT':'IIIC', 'HGSOC5_AS':'IIIC', 'HGSOC4_LN':'IIIC', 'HGSOC4_MT':'IIIC',
       'HGSOC4_BC':'IIIC', 'HGSOC4_PT':'IIIC', 'HGSOC8_PT':'IIIC', 'HGSOC8_AS':'IIIC', 'HGSOC9_PT':'IIIC',
       'HGSOC9_AS':'IIIC', 'ECO1_BC':'IV', 'ECO1_PT':'IV', 'ECO1_MT':'IV', 'ECO1_LN':'IV', 'ECO1_AS':'IV',
       'HGSOC10_AS':'IIIC', 'HGSOC10_PT':'IIIC', 'C1_PT':'IIIC'}

columns['paper_ID'] = {'HGSOC1_BC':'Zheng_1', 'HGSOC1_PT':'Zheng_1', 'HGSOC1_MT':'Zheng_1', 'HGSOC1_AS':'Zheng_1', 'HGSOC3_BC':'Zheng_3',
       'HGSOC3_PT':'Zheng_3', 'HGSOC3_MT':'Zheng_3', 'HGSOC3_LN':'Zheng_3', 'HGSOC3_AS':'Zheng_3', 'HGSOC2_BC':'Zheng_2',
       'HGSOC2_PT':'Zheng_2', 'HGSOC2_LN':'Zheng_2', 'HGSOC2_AS':'Zheng_2', 'HGSOC7_PT':'Zheng_7', 'UOC1_PT':'Zheng_Undifferentiated',
       'HGSOC6_AS':'Zheng_6', 'HGSOC6_LN':'Zheng_6', 'HGSOC6_MT':'Zheng_6', 'HGSOC6_BC':'Zheng_6', 'HGSOC6_PT':'Zheng_6',
       'OCCC1_AS':'Zheng_Clear cell carcinoma', 'HGSOC5_PT':'Zheng_5', 'HGSOC5_AS':'Zheng_5', 'HGSOC4_LN':'Zheng_4', 'HGSOC4_MT':'Zheng_4',
       'HGSOC4_BC':'Zheng_4', 'HGSOC4_PT':'Zheng_4', 'HGSOC8_PT':'Zheng_8', 'HGSOC8_AS':'Zheng_8', 'HGSOC9_PT':'Zheng_9',
       'HGSOC9_AS':'Zheng_9', 'ECO1_BC':'Zheng_Endometroid_carcinoma', 'ECO1_PT':'Zheng_Endometroid_carcinoma', 'ECO1_MT':'Zheng_Endometroid_carcinoma', 'ECO1_LN':'Zheng_Endometroid_carcinoma', 'ECO1_AS':'Zheng_Endometroid_carcinoma',
       'HGSOC10_AS':'Zheng_10', 'HGSOC10_PT':'Zheng_10', 'C1_PT':'Zheng_Carcinosarcoma'}

columns['anatomical_location'] = {'HGSOC1_BC':'Blood', 'HGSOC1_PT':'Ovary', 'HGSOC1_MT':'Omentum', 'HGSOC1_AS':'Ascites', 'HGSOC3_BC':'Blood',
       'HGSOC3_PT':'Ovary', 'HGSOC3_MT':'Omentum', 'HGSOC3_LN':'Lymph node', 'HGSOC3_AS':'Ascites', 'HGSOC2_BC':'Blood',
       'HGSOC2_PT':'Ovary', 'HGSOC2_LN':'Lymph node', 'HGSOC2_AS':'Ascites', 'HGSOC7_PT':'Ovary', 'UOC1_PT':'Ovary',
       'HGSOC6_AS':'Ascites', 'HGSOC6_LN':'Lymph node', 'HGSOC6_MT':'Omentum', 'HGSOC6_BC':'Blood', 'HGSOC6_PT':'Ovary',
       'OCCC1_AS':'Ascites', 'HGSOC5_PT':'Ovary', 'HGSOC5_AS':'Ascites', 'HGSOC4_LN':'Lymph node', 'HGSOC4_MT':'Omentum',
       'HGSOC4_BC':'Blood', 'HGSOC4_PT':'Ovary', 'HGSOC8_PT':'Ovary', 'HGSOC8_AS':'Ascites', 'HGSOC9_PT':'Ovary',
       'HGSOC9_AS':'Ascites', 'ECO1_BC':'Blood', 'ECO1_PT':'Ovary', 'ECO1_MT':'Omentum', 'ECO1_LN':'Lymph node', 'ECO1_AS':'Ascites',
       'HGSOC10_AS':'Ascites', 'HGSOC10_PT':'Ovary', 'C1_PT':'Ovary'}

for column in columns:
    adata.obs[column] = adata.obs.Samples.replace(columns[column])

adata.obs.rename(columns = {'Samples':'sample_name'}, inplace = True)
adata.obs.drop(columns=['orig.ident', 'Groups', 'Patients', 'Group_abb'], inplace=True)
#%%
adata.obs['dataset'] = 'Zheng'

adata.obs.index = adata.obs_names
   
#%%
#Write raw adata with metadata

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

adata.write_h5ad(out_dir + "zheng2023_wholedataset_rawcounts.h5ad")

#%%
##Preprocessing

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(out_dir + "zheng2023_wholedataset_filt_norm_nolog.h5ad")

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

adata.write(out_dir + "zheng2023_wholedataset_embeddings.h5ad")
#%%
#Remove blood, lymph nodes and all non-HGSOC cancer types
adata = sc.read(out_dir + "zheng2023_wholedataset_rawcounts.h5ad")
#%%
#adata.obs.drop(columns=['n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'], inplace=True)

mask = ((adata.obs['tissue'] != 'Blood') & (adata.obs['tissue'] != 'Lymph node'))
mask.sum()
(adata.obs['tissue'] == 'Blood').sum()
(adata.obs['tissue'] == 'Lymph node').sum()
adata = adata[mask]
adata.obs

mask2 = ((adata.obs['paper_ID'] != 'Zheng_Carcinosarcoma') & 
        (adata.obs['paper_ID'] != 'Zheng_Clear cell carcinoma') & 
        (adata.obs['paper_ID'] != 'Zheng_Endometroid_carcinoma') & 
        (adata.obs['paper_ID'] != 'Zheng_Undifferentiated'))
adata = adata[mask2]

adata.write(out_dir + "zheng2023_rawcounts.h5ad")
#%%
##Preprocessing

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(out_dir + "zheng2023_filt_norm_nolog.h5ad")

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

adata.write(out_dir + "zheng2023_embeddings.h5ad")

# %%
