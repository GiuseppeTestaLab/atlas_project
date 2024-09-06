
#library imports

import scanpy as sc
import pandas as pd
import sys
import os
import zipfile
import tempfile
import glob
import argparse
import numpy as np
import anndata
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read('../../../utils/config.ini')
scriptsPath = config.get('DEFAULT', 'scriptsPath')
#inputs
parameters = pd.read_csv(scriptsPath+ "1_original_counts/Xu2022/preprocess_params.csv", sep = ';')

init_dir = rowPath+parameters.init_dir[0]
out_dir = rowPath+parameters.out_dir[0]
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum= float(parameters.target_sum[0])

run = glob.glob(os.path.join(init_dir, "Cancer/*"))
                 
files = [os.path.join(x, "feature_bc_matrix") for x in run]
samplenames = ["Cancer1", "Cancer2", "Cancer3", "Cancer4", "Cancer5", "Cancer6", "Cancer7"]
filesdict =  dict(zip(samplenames, files))

#create anndata

#Paper: Xu J, Fang Y, Chen K, Li S et al. Single-cell RNA sequencing reveals the tissue architecture in human high-grade serous ovarian cancer. Clin Cancer Res 2022 Jun 8. PMID: 35675036

#Counts file: standard 10X files (matrix + features + barcodes)

adata=['']*len(files)
for i, item in enumerate(files):
    print(i)
    print(item)
    filename_data = os.path.join(item, 'matrix.mtx.gz')
    filename_genes = os.path.join(item, 'genes.tsv.gz')
    filename_barcodes = os.path.join(item, 'barcodes.tsv.gz')
    adata[i] = sc.read(filename_data, cache=True).transpose()
    adata[i].var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
    adata[i].obs_names = np.genfromtxt(filename_barcodes, dtype=str)
    adata[i].obs['sample_ID'] = samplenames[i]
    print("SAMPLENAME:", samplenames[i], sep = "\t")
    print(adata[i].X.shape)
    adata[i].var_names_make_unique()


adata = adata[0].concatenate(adata[1:], batch_categories=samplenames)

# Creating metadata

columns = {}

columns['sample_name'] = {"Cancer1": "GSM5599225", "Cancer2":"GSM5599226", "Cancer3":"GSM5599227", "Cancer4":"GSM5599228", "Cancer5":"GSM5599229", "Cancer6":"GSM5599230", "Cancer7":"GSM5599231"}
columns['tissue'] = {"Cancer1": "Primary", "Cancer2":"Primary", "Cancer3":"Primary", "Cancer4":"Primary", "Cancer5":"Primary", "Cancer6":"Primary", "Cancer7":"Primary"}
columns['developmental_stage'] = {"Cancer1": "Advanced stage", "Cancer2":"Advanced stage", "Cancer3":"Early stage", "Cancer4":"Early stage", "Cancer5":"Advanced stage", "Cancer6":"Advanced stage", "Cancer7":"Early stage"}
columns['treatment'] = {"Cancer1": "Naive", "Cancer2":"Naive", "Cancer3":"Naive", "Cancer4":"Naive", "Cancer5":"Naive", "Cancer6":"Naive", "Cancer7":"Naive"}
columns['recurrence'] = {"Cancer1": "Sensitive", "Cancer2":"Sensitive", "Cancer3":"Sensitive", "Cancer4":"Sensitive", "Cancer5":"Sensitive", "Cancer6":"Recurrence", "Cancer7":"Sensitive"}
columns['tumor_stage'] = {"Cancer1": "IIIB", "Cancer2":"IIB", "Cancer3":"IC", "Cancer4":"IC", "Cancer5":"IIB", "Cancer6":"IIIB", "Cancer7":"IC"}
columns['paper_ID'] = {"Cancer1": "Xu_1", "Cancer2":"Xu_2", "Cancer3":"Xu_3", "Cancer4":"Xu_4", "Cancer5":"Xu_5", "Cancer6":"Xu_6", "Cancer7":"Xu_7"}
columns['anatomical_location'] = {"Cancer1": "Ovary", "Cancer2":"Ovary", "Cancer3":"Ovary", "Cancer4":"Ovary", "Cancer5":"Ovary", "Cancer6":"Ovary", "Cancer7":"Ovary"}

for column in columns:
    adata.obs[column] = adata.obs.batch.replace(columns[column])
    

adata.obs.rename(columns = {'sample_ID':'ID'}, inplace = True)
del adata.obs['batch']
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 8000, "pct_counts_mt" : 20, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(out_dir + "xu2022_rawcounts.h5ad")

#Preprocessing

final_dir = parameters.final_dir[0]

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(final_dir + "Xu2022_filt_norm_nolog.h5ad")