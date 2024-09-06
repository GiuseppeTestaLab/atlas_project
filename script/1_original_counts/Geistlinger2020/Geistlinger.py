
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
parameters = pd.read_csv(scriptsPath+ "1_original_counts/Geistlinger2020/preprocess_params.csv", sep = ';')

init_dir = rowPath+parameters.init_dir[0]
out_dir = rowPath+parameters.out_dir[0]
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum= float(parameters.target_sum[0])
                
files = glob.glob(os.path.join(init_dir, "Geistlinger2020/Samples/*"))
samplenames = ["T59", "T76", "T77", "T89", "T90"]
filesdict =  dict(zip(samplenames, files))

#create anndata

#Paper: Geistlinger L et al. Multiomic Analysis of Subtype Evolution and Heterogeneity in High-Grade Serous Ovarian Carcinoma. Cancer Res 2020 Oct 15. PMID: 32747365

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


columns['sample_name'] = {"T59": "GSM4675273", "T76":"GSM4675274", "T77":"GSM4675275", "T89":"GSM4675276", "T90":"GSM4675277"}

columns['tissue'] = {"T59": "Metastasis", "T76":"Metastasis", "T77":"Metastasis", "T89":"Metastasis", "T90":"Metastasis"}

columns['developmental_stage'] = {"T59": "Advanced stage", "T76":"Advanced stage", "T77":"Advanced stage", "T89":"Advanced stage", "T90":"Advanced stage"}

columns['treatment'] = {"T59": "CHT", "T76":"CHT", "T77":"CHT", "T89":"CHT", "T90":"CHT"}

columns['recurrence'] = {"T59": "Recurrence", "T76":"Recurrence", "T77":"Recurrence", "T89":"Sensitive", "T90":"Sensitive"}

columns['tumor_stage'] = {"T59": "IV", "T76":"IIIC", "T77":"IVB", "T89":"IV", "T90":"IVA"}

columns['paper_ID'] = {"T59": "Geistlinger_59", "T76":"Geistlinger_76", "T77":"Geistlinger_77", "T89":"Geistlinger_89", "T90":"Geistlinger_90"}

columns['anatomical_location'] = {"T59": "Omentum", "T76":"Omentum", "T77":"Omentum", "T89":"Omentum", "T90":"Omentum"}

for column in columns:
    adata.obs[column] = adata.obs.batch.replace(columns[column])
    

adata.obs.rename(columns = {'sample_ID':'ID'}, inplace = True)
del adata.obs['batch']
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 5000, "pct_counts_mt" : 15, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(out_dir + "geistlinger2020_rawcounts.h5ad")

#Preprocessing

final_dir = parameters.final_dir[0]

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(final_dir + "Geistlinger2020_filt_norm_nolog.h5ad")

