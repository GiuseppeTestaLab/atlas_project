
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

#inputs
parameters = pd.read_csv(scriptsPath+ "1_original_counts/Olbrecht2021/preprocess_params.csv", sep = ';')

init_dir = rowPath+parameters.init_dir[0]
out_dir = rowPath+parameters.out_dir[0]
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum= float(parameters.target_sum[0])

                 
files = glob.glob(os.path.join(init_dir, "Olbrecht2021/10xCounts"))
samplenames = ["SOL1303", "SOL1304", "SOL1305", "SOL1306", "SOL1307", "SOL003", "SOL004", "SOL006", "SOL007", "SOL008", "SOL012", "SOL016"]
filesdict =  dict(zip(samplenames, files))

#create anndata

#Paper: Olbrecht S et al. High-grade serous tubo-ovarian cancer refined with single-cell RNA sequencing: specific cell subtypes influence survival and determine molecular subtype classification. Genome Medicine 2021 Jul 9. PMID: 34238352

#Counts file: standard 10X files (matrix + features + barcodes)

filename_data = os.path.join('/group/testa/Project/OvarianAtlas/Olbrecht2021/10xCounts', 'matrix.mtx')
filename_genes = os.path.join('/group/testa/Project/OvarianAtlas/Olbrecht2021/10xCounts', 'genes.tsv')
filename_barcodes = os.path.join('/group/testa/Project/OvarianAtlas/Olbrecht2021/10xCounts', 'barcodes.tsv')
adata = sc.read(filename_data, cache=True).transpose()
adata.var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
adata.obs_names = np.genfromtxt(filename_barcodes, dtype=str)

samples = []

for i in adata.obs_names:
    name = i.split("_")[1]
    samples.append(name)

adata.obs['sample_ID'] = samples

adata.var_names_make_unique()

# Creating metadata

columns = {}

columns['patient_id'] = {"SOL1303": "P1", "SOL1304":"P1", "SOL1305":"P1", "SOL1306":"P1", "SOL1307":"P1", "SOL003":"P2", "SOL004":"P3", "SOL006":"P4", "SOL007":"P4", "SOL008":"P5", "SOL012":"P6", "SOL016":"P7"}

columns['sample_name'] = {"SOL1303": "EGAD00001006898", "SOL1304":"EGAD00001006898", "SOL1305":"EGAD00001006898", "SOL1306":"EGAD00001006898", "SOL1307":"EGAD00001006898", "SOL003":"EGAD00001006898", "SOL004":"EGAD00001006898", "SOL006":"EGAD00001006898", "SOL007":"EGAD00001006898", "SOL008":"EGAD00001006898", "SOL012":"EGAD00001006898", "SOL016":"EGAD00001006898"}

columns['tissue'] = {"SOL1303": "Metastasis", "SOL1304":"Normal", "SOL1305":"Normal", "SOL1306":"Metastasis", "SOL1307":"Primary", "SOL003":"Metastasis", "SOL004":"Metastasis", "SOL006":"Primary", "SOL007":"Normal", "SOL008":"Metastasis", "SOL012":"Metastasis", "SOL016":"Metastasis"}

columns['developmental_stage'] = {"SOL1303": "Advanced stage", "SOL1304":"Advanced stage", "SOL1305":"Advanced stage", "SOL1306":"Advanced stage", "SOL1307":"Advanced stage", "SOL003":"Advanced stage", "SOL004":"Advanced stage", "SOL006":"Early stage", "SOL007":"Early stage", "SOL008":"Advanced stage", "SOL012":"Advanced stage", "SOL016":"Advanced stage"}

columns['treatment'] = {"SOL1303": "Naive", "SOL1304":"Naive", "SOL1305":"Naive", "SOL1306":"Naive", "SOL1307":"Naive", "SOL003":"Naive", "SOL004":"Naive", "SOL006":"Naive", "SOL007":"Naive", "SOL008":"Naive", "SOL012":"Naive", "SOL016":"Naive"}

columns['recurrence'] = {"SOL1303": "Sensitive", "SOL1304":"Sensitive", "SOL1305":"Sensitive", "SOL1306":"Sensitive", "SOL1307":"Sensitive", "SOL003":"Sensitive", "SOL004":"Recurrence", "SOL006":"Sensitive", "SOL007":"Sensitive", "SOL008":"Recurrence", "SOL012":"Sensitive", "SOL016":"Recurrence"}

columns['tumor_stage'] = {"SOL1303": "IIIC", "SOL1304":"IIIC", "SOL1305":"IIIC", "SOL1306":"IIIC", "SOL1307":"IIIC", "SOL003":"IVB", "SOL004":"IVB", "SOL006":"IC", "SOL007":"IC", "SOL008":"IVB", "SOL012":"IIIC", "SOL016":"IVB"}

columns['paper_ID'] = {"SOL1303": "Olbrecht_1", "SOL1304":"Olbrecht_1", "SOL1305":"Olbrecht_1", "SOL1306":"Olbrecht_1", "SOL1307":"Olbrecht_1", "SOL003":"Olbrecht_2", "SOL004":"Olbrecht_3", "SOL006":"Olbrecht_4", "SOL007":"Olbrecht_4", "SOL008":"Olbrecht_5", "SOL012":"Olbrecht_6", "SOL016":"Olbrecht_7"}

columns['anatomical_location'] = {"SOL1303": "Omentum", "SOL1304":"Omentum", "SOL1305":"Peritoneum", "SOL1306":"Peritoneum", "SOL1307":"Ovary", "SOL003":"Peritoneum", "SOL004":"Peritoneum", "SOL006":"Ovary", "SOL007":"Ovary", "SOL008":"Peritoneum", "SOL012":"Peritoneum", "SOL016":"Peritoneum"}

for column in columns:
    adata.obs[column] = adata.obs.sample_ID.replace(columns[column])
    

adata.obs.rename(columns = {'sample_ID':'ID'}, inplace = True)
adata = adata[adata.obs.tissue != "Normal"]
adata = adata[adata.obs.ID != "SOL004"]
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 6000, "pct_counts_mt" : 20, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(out_dir + "olbrecht2021_rawcounts.h5ad")

#Preprocessing

final_dir = parameters.final_dir[0]

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(final_dir + "olbrecht2021_filt_norm_nolog.h5ad")