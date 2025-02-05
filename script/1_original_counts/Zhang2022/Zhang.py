# library imports

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
    scriptsPath + "1_original_counts/Zhang2022/preprocess_params.tsv", sep=";"
)

init_dir = os.path.join(rawPath, parameters.init_dir[0])
out_dir = os.path.join(rawPath, parameters.out_dir[0])
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum = float(parameters.target_sum[0])

# create anndata

# Paper: Zhang K, Erkan EP, Jamalzadeh S, Dai J et al. Longitudinal single-cell RNA-seq analysis reveals stress-promoted chemoresistance in metastatic ovarian cancer. Sci Adv 2022 Feb 25;8(8):eabm1831. PMID: 35196078
# Counts file: .tsv

counts = pd.read_table(
    init_dir + "GSE165897_UMIcounts_HGSOC.tsv", sep="\t", index_col=0
)
adata = ad.AnnData(X=counts.transpose())
metadata = pd.read_csv(init_dir + "GSE165897_cellInfo_HGSOC.tsv", sep="\t")
metadata = metadata.set_index("cell")
adata.obs = metadata
adata.obs.rename(columns={"sample": "ID"}, inplace=True)
# adata.obs['batch'] = adata.obs['ID']
adata.obs = adata.obs.drop(columns=["nCount_RNA", "nFeature_RNA", "percent.mt"])

# Creating metadata

columns = {}

columns["sample_name"] = {
    "EOC372_primary_Peritoneum": "GSM5057586",
    "EOC372_interval_Peritoneum": "GSM5057587",
    "EOC443_interval_Omentum": "GSM5057591",
    "EOC443_primary_Omentum": "GSM5057590",
    "EOC540_primary_Omentum": "GSM5057592",
    "EOC540_interval_Omentum": "GSM5057593",
    "EOC3_interval_Omentum": "GSM5057589",
    "EOC3_primary_Peritoneum": "GSM5057588",
    "EOC87_interval_Omentum": "GSM5057597",
    "EOC87_primary_Peritoneum": "GSM5057596",
    "EOC136_interval_Omentum": "GSM5057579",
    "EOC136_primary_Mesentery": "GSM5057578",
    "EOC1005_primary_Peritoneum": "GSM5057576",
    "EOC1005_interval_Tumor": "GSM5057577",
    "EOC733_primary_Peritoneum": "GSM5057594",
    "EOC733_interval_Omentum": "GSM5057595",
    "EOC153_interval_Omentum": "GSM5057581",
    "EOC153_primary_Omentum": "GSM5057580",
    "EOC349_interval_Omentum": "GSM5057585",
    "EOC349_primary_Peritoneum": "GSM5057584",
    "EOC227_interval_Omentum": "GSM5057583",
    "EOC227_primary_Omentum": "GSM5057582",
}

columns["tissue"] = {
    "EOC372_primary_Peritoneum": "Metastasis",
    "EOC372_interval_Peritoneum": "Metastasis",
    "EOC443_interval_Omentum": "Metastasis",
    "EOC443_primary_Omentum": "Metastasis",
    "EOC540_primary_Omentum": "Metastasis",
    "EOC540_interval_Omentum": "Metastasis",
    "EOC3_interval_Omentum": "Metastasis",
    "EOC3_primary_Peritoneum": "Metastasis",
    "EOC87_interval_Omentum": "Metastasis",
    "EOC87_primary_Peritoneum": "Metastasis",
    "EOC136_interval_Omentum": "Metastasis",
    "EOC136_primary_Mesentery": "Metastasis",
    "EOC1005_primary_Peritoneum": "Metastasis",
    "EOC1005_interval_Tumor": "Metastasis",
    "EOC733_primary_Peritoneum": "Metastasis",
    "EOC733_interval_Omentum": "Metastasis",
    "EOC153_interval_Omentum": "Metastasis",
    "EOC153_primary_Omentum": "Metastasis",
    "EOC349_interval_Omentum": "Metastasis",
    "EOC349_primary_Peritoneum": "Metastasis",
    "EOC227_interval_Omentum": "Metastasis",
    "EOC227_primary_Omentum": "Metastasis",
}

columns["developmental_stage"] = {
    "EOC372_primary_Peritoneum": "Advanced stage",
    "EOC372_interval_Peritoneum": "Advanced stage",
    "EOC443_interval_Omentum": "Advanced stage",
    "EOC443_primary_Omentum": "Advanced stage",
    "EOC540_primary_Omentum": "Advanced stage",
    "EOC540_interval_Omentum": "Advanced stage",
    "EOC3_interval_Omentum": "Advanced stage",
    "EOC3_primary_Peritoneum": "Advanced stage",
    "EOC87_interval_Omentum": "Advanced stage",
    "EOC87_primary_Peritoneum": "Advanced stage",
    "EOC136_interval_Omentum": "Advanced stage",
    "EOC136_primary_Mesentery": "Advanced stage",
    "EOC1005_primary_Peritoneum": "Advanced stage",
    "EOC1005_interval_Tumor": "Advanced stage",
    "EOC733_primary_Peritoneum": "Advanced stage",
    "EOC733_interval_Omentum": "Advanced stage",
    "EOC153_interval_Omentum": "Advanced stage",
    "EOC153_primary_Omentum": "Advanced stage",
    "EOC349_interval_Omentum": "Advanced stage",
    "EOC349_primary_Peritoneum": "Advanced stage",
    "EOC227_interval_Omentum": "Advanced stage",
    "EOC227_primary_Omentum": "Advanced stage",
}

columns["recurrence"] = {
    "EOC372_primary_Peritoneum": "Sensitive",
    "EOC372_interval_Peritoneum": "Sensitive",
    "EOC443_interval_Omentum": "Sensitive",
    "EOC443_primary_Omentum": "Sensitive",
    "EOC540_primary_Omentum": "Sensitive",
    "EOC540_interval_Omentum": "Sensitive",
    "EOC3_interval_Omentum": "Sensitive",
    "EOC3_primary_Peritoneum": "Sensitive",
    "EOC87_interval_Omentum": "Sensitive",
    "EOC87_primary_Peritoneum": "Sensitive",
    "EOC136_interval_Omentum": "Sensitive",
    "EOC136_primary_Mesentery": "Sensitive",
    "EOC1005_primary_Peritoneum": "Sensitive",
    "EOC1005_interval_Tumor": "Sensitive",
    "EOC733_primary_Peritoneum": "Sensitive",
    "EOC733_interval_Omentum": "Sensitive",
    "EOC153_interval_Omentum": "Sensitive",
    "EOC153_primary_Omentum": "Sensitive",
    "EOC349_interval_Omentum": "Sensitive",
    "EOC349_primary_Peritoneum": "Sensitive",
    "EOC227_interval_Omentum": "Sensitive",
    "EOC227_primary_Omentum": "Sensitive",
}

columns["tumor_stage"] = {
    "EOC372_primary_Peritoneum": "IIIC",
    "EOC372_interval_Peritoneum": "IIIC",
    "EOC443_interval_Omentum": "IVA",
    "EOC443_primary_Omentum": "IVA",
    "EOC540_primary_Omentum": "IIIC",
    "EOC540_interval_Omentum": "IIIC",
    "EOC3_interval_Omentum": "IVA",
    "EOC3_primary_Peritoneum": "IVA",
    "EOC87_interval_Omentum": "IIIC",
    "EOC87_primary_Peritoneum": "IIIC",
    "EOC136_interval_Omentum": "IVA",
    "EOC136_primary_Mesentery": "IVA",
    "EOC1005_primary_Peritoneum": "IVA",
    "EOC1005_interval_Tumor": "IVA",
    "EOC733_primary_Peritoneum": "IVA",
    "EOC733_interval_Omentum": "IVA",
    "EOC153_interval_Omentum": "IVA",
    "EOC153_primary_Omentum": "IVA",
    "EOC349_interval_Omentum": "IVB",
    "EOC349_primary_Peritoneum": "IVB",
    "EOC227_interval_Omentum": "IVA",
    "EOC227_primary_Omentum": "IVA",
}

columns["paper_ID"] = {
    "EOC372_primary_Peritoneum": "Zhang_372",
    "EOC372_interval_Peritoneum": "Zhang_372",
    "EOC443_interval_Omentum": "Zhang_443",
    "EOC443_primary_Omentum": "Zhang_443",
    "EOC540_primary_Omentum": "Zhang_540",
    "EOC540_interval_Omentum": "Zhang_540",
    "EOC3_interval_Omentum": "Zhang_3",
    "EOC3_primary_Peritoneum": "Zhang_3",
    "EOC87_interval_Omentum": "Zhang_87",
    "EOC87_primary_Peritoneum": "Zhang_87",
    "EOC136_interval_Omentum": "Zhang_136",
    "EOC136_primary_Mesentery": "Zhang_136",
    "EOC1005_primary_Peritoneum": "Zhang_1005",
    "EOC1005_interval_Tumor": "Zhang_1005",
    "EOC733_primary_Peritoneum": "Zhang_733",
    "EOC733_interval_Omentum": "Zhang_733",
    "EOC153_interval_Omentum": "Zhang_153",
    "EOC153_primary_Omentum": "Zhang_153",
    "EOC349_interval_Omentum": "Zhang_349",
    "EOC349_primary_Peritoneum": "Zhang_349",
    "EOC227_interval_Omentum": "Zhang_227",
    "EOC227_primary_Omentum": "Zhang_227",
}


for column in columns:
    adata.obs[column] = adata.obs.ID.replace(columns[column])

adata.obs["dataset"] = adata.obs.paper_ID.str.split("_").str[0]

treat = []
for i in adata.obs.treatment_phase:
    if i == "treatment-naive":
        treat.append("Naive")
    elif i == "post-NACT":
        treat.append("NACT")
    else:
        treat.append(i)
adata.obs["treatment"] = treat

adata.obs = adata.obs.drop(columns=["treatment_phase"])

adata.obs = adata.obs[
    [
        "ID",
        "patient_id",
        "sample_name",
        "tissue",
        "developmental_stage",
        "treatment",
        "recurrence",
        "tumor_stage",
        "paper_ID",
        "anatomical_location",
        "dataset",
        "cell_type",
        "cell_subtype",
    ]
]

adata.uns["preprocessing"] = {
    "min_genes": 200,
    "min_cells": 3,
    "genes_by_counts": 5000,
    "pct_counts_mt": 8,
    "target_sum": 1e4,
}

# Write raw adata with metadata
adata.write_h5ad(out_dir + "zhang2022_rawcounts.h5ad")

# Preprocessing

final_dir = os.path.join(rawPath, parameters.final_dir[0])

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var["mt"] = adata.var_names.str.startswith(
    "MT-"
)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

if not os.path.exists(final_dir):
    os.makedirs(final_dir)

adata.write(final_dir + "Zhang2022_filt_norm_nolog.h5ad")
