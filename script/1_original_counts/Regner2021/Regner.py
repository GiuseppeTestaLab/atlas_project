
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
rawPath = config.get('DEFAULT', 'rawPath')

#inputs
parameters = pd.read_csv(scriptsPath+ "1_original_counts/Regner2021/preprocess_params.tsv", sep = ';')

init_dir = rawPath+parameters.init_dir[0]
out_dir = rawPath+parameters.out_dir[0]
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum= float(parameters.target_sum[0])


run = glob.glob(os.path.join(init_dir, "HGSOC_tum/*"))
                 
files = [os.path.join(x, "feature_bc_matrix") for x in run]
samplenames = ["Patient_8", "Patient_9"]
filesdict =  dict(zip(samplenames, files))

#create anndata

#Paper: Regner MJ, Wisniewska K, Garcia-Recio S, Thennavan A et al. A multi-omic single-cell landscape of human gynecologic malignancies. Mol Cell 2021 Dec 2;81(23):4924-4941.e10. PMID: 34739872
#Counts file: standard 10X files (matrix + features + barcodes)

adata=['']*len(files)
for i, item in enumerate(files):
    print(i)
    print(item)
    filename_data = os.path.join(item, 'matrix.mtx')
    filename_genes = os.path.join(item, 'features.tsv')
    filename_barcodes = os.path.join(item, 'barcodes.tsv')
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

columns['sample_name'] = {"Patient_8":"GSM5276940", "Patient_9":"GSM5276943"}
columns['tissue'] = {"Patient_8":"Primary", "Patient_9":"Primary"}
columns['developmental_stage'] = {"Patient_8":"Advanced stage", "Patient_9":"Advanced stage"}
columns['treatment'] = {"Patient_8":"Naive", "Patient_9":"Naive"}
columns['recurrence'] = {"Patient_8":"Unknown", "Patient_9":"Unknown"}
columns['tumor_stage'] = {"Patient_8":"IIB", "Patient_9":"IIIC"}
columns['paper_ID'] = {"Patient_8":"Regner_8", "Patient_9":"Regner_9"}
columns['anatomical_location'] = {"Patient_8":"Ovary", "Patient_9":"Ovary"}

for column in columns:
    adata.obs[column] = adata.obs.batch.replace(columns[column])
    

adata.obs.rename(columns = {'sample_ID':'ID'}, inplace = True)
del adata.obs['batch']
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 8000, "pct_counts_mt" : 20, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(out_dir + "regner2021_rawcounts.h5ad")

#Preprocessing

final_dir = parameters.final_dir[0]

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(final_dir + "Regner2021_filt_norm_nolog.h5ad")