
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

#inputs

args = sys.argv

initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Ren2022/" 
outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Ren2022/Adata/"
min_genes = int(args[1])
min_cells = int(args[2])
genes_by_counts = int(args[3])
pct_counts_mt = float(args[4])
target_sum= float(args[5])


run = glob.glob(os.path.join(initDir, "Counts/*"))
                 
files = [os.path.join(x, "feature_bc_matrix") for x in run]
samplenames = ["Patient_1_as", "Patient_1_tum"]
filesdict =  dict(zip(samplenames, files))

#create anndata

#Paper: Ren Y et al. Single-cell sequencing reveals effects of chemotherapy on the immune landscape and TCR/BCR clonal expansion in a relapsed ovarian cancer patients. Frontiers in Immunology 2022 Sep 28. PMID: 36248860
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

columns['ID'] = {"Patient_1_tum":"Patient_1", "Patient_1_as":"Patient_1"}
columns['sample_name'] = {"Patient_1_tum":"GSM6576367", "Patient_1_as":"GSM6576368"}
columns['tissue'] = {"Patient_1_tum":"Primary", "Patient_1_as":"Ascites"}
columns['developmental_stage'] = {"Patient_1_tum":"Advanced stage", "Patient_1_as":"Advanced stage"}
columns['treatment'] = {"Patient_1_tum":"CHT", "Patient_1_as":"CHT"}
columns['recurrence'] = {"Patient_1_tum":"Recurrence", "Patient_1_as":"Recurrence"}
columns['tumor_stage'] = {"Patient_1_tum":"IIIC", "Patient_1_as":"IIIC"}
columns['paper_ID'] = {"Patient_1_tum":"Ren_1", "Patient_1_as":"Ren_1"}
columns['anatomical_location'] = {"Patient_1_tum":"Ovary", "Patient_1_as":"Ascites"}

for column in columns:
    adata.obs[column] = adata.obs.batch.replace(columns[column])
    
del adata.obs['batch']
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 10000, "pct_counts_mt" : 30, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(outDir + "ren2022_rawcounts.h5ad")

#Preprocessing

finalDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Ren2022/"

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(finalDir + "ren2022_filt_norm_nolog.h5ad")
