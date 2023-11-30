
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
import anndata as ad

#inputs

args = sys.argv

initDir = "/group/testa/Project/OvarianAtlas/Qian2020/atlas_project/raw_data/original_counts/" 
outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Adata/"
min_genes = int(args[1])
min_cells = int(args[2])
genes_by_counts = int(args[3])
pct_counts_mt = float(args[4])
target_sum= float(args[5])

listpath = ["/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/BT1305.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/BT1306.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/BT1307.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/scrSOL001.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/scrSOL003.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/scrSOL004.counts.csv",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/Counts/scrSOL006.counts.csv"]

#create anndata

counter = 1

counteradata = 0

adata = ['']*len(listpath)

for i in listpath:
    print(counteradata + 1) 
    print(i)
    counts = pd.read_csv(i, index_col=0)
    counts = counts.T
    adata[counteradata] = ad.AnnData(X = counts)
    print(adata[counteradata].X.shape)
    adata[counteradata].var_names_make_unique()
    counteradata = counteradata + 1
    counter = counter + 1

adata = adata[0].concatenate(adata[1:], index_unique=None)

#Paper: Qian J et al. A pan-cancer blueprint of the heterogeneous tumor microenvironment revealed by single-cell profiling. Cell Res 2020 Sep;30(9):745-762. PMID: 32561858
#Counts file: csv

# Creating metadata

samples = []

for a in adata.obs_names:
    name = a.split('_')[0]
    samples.append(name)
adata.obs['sample_name'] = samples
adata.obs['batch'] = adata.obs['sample_name']

columns = {}

columns['ID'] = {"BT1305":"P11_BT1305", 
                "BT1306":"P11_BT1306",
                "BT1307":"P11_BT1307",
                "scrSOL001":"P12_SOL001",
                "scrSOL003":"P13_SOL003",
                "scrSOL004":"P14_SOL004",
                "scrSOL006":"P15_SOL006"}

columns['tissue'] = {"BT1305":"Metastasis", 
                    "BT1306":"Metastasis",
                    "BT1307":"Primary",
                    "scrSOL001":"Metastasis",
                    "scrSOL003":"Metastasis",
                    "scrSOL004":"Metastasis",
                    "scrSOL006":"Primary"}

columns['developmental_stage'] = {"BT1305":"Advanced stage", 
                                    "BT1306":"Advanced stage",
                                    "BT1307":"Advanced stage",
                                    "scrSOL001":"Advanced stage",
                                    "scrSOL003":"Advanced stage",
                                    "scrSOL004":"Advanced stage",
                                    "scrSOL006":"Early stage"}

columns['treatment'] = {"BT1305":"Naive", 
                    "BT1306":"Naive",
                    "BT1307":"Naive",
                    "scrSOL001":"Naive",
                    "scrSOL003":"Naive",
                    "scrSOL004":"Naive",
                    "scrSOL006":"Naive"}

columns['recurrence'] = {"BT1305":"Unknown", 
                    "BT1306":"Unknown",
                    "BT1307":"Unknown",
                    "scrSOL001":"Unknown",
                    "scrSOL003":"Unknown",
                    "scrSOL004":"Unknown",
                    "scrSOL006":"Unknown"}

columns['tumor_stage'] = {"BT1305":"IIIC", 
                    "BT1306":"IIIC",
                    "BT1307":"IIIC",
                    "scrSOL001":"IVB",
                    "scrSOL003":"IVB",
                    "scrSOL004":"IVB_mixedCCC",
                    "scrSOL006":"IA"}

columns['paper_ID'] = {"BT1305":"Qian_11", 
                    "BT1306":"Qian_11",
                    "BT1307":"Qian_11",
                    "scrSOL001":"Qian_12",
                    "scrSOL003":"Qian_13",
                    "scrSOL004":"Qian_14",
                    "scrSOL006":"Qian_15"}

columns['anatomical_location'] = {"BT1305":"Peritoneum", 
                                    "BT1306":"Peritoneum",
                                    "BT1307":"Ovary",
                                    "scrSOL001":"Peritoneum",
                                    "scrSOL003":"Peritoneum",
                                    "scrSOL004":"Peritoneum",
                                    "scrSOL006":"Ovary"}

for column in columns:
    adata.obs[column] = adata.obs.batch.replace(columns[column])
    

del adata.obs['batch']
adata = adata[adata.obs.sample_name != "scrSOL004"]
adata = adata[(adata.obs.ID == "P12_SOL001") |
          (adata.obs.ID == "P11_BT1305")]
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 4000, "pct_counts_mt" : 12, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(outDir + "qian2020_rawcounts.h5ad")

#Preprocessing

finalDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Qian2020/"

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(finalDir + "qian2020_filt_norm_nolog.h5ad")
