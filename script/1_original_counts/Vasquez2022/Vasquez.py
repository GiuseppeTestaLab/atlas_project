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

#inputs
parameters = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/1_original_counts/Vasquez2022/preprocess_params.csv', sep = ';')

init_dir = parameters.init_dir[0]
out_dir = parameters.out_dir[0]
min_genes = int(parameters.min_genes[0])
min_cells = int(parameters.min_cells[0])
genes_by_counts = int(parameters.genes_by_counts[0])
pct_counts_mt = float(parameters.pct_counts_mt[0])
target_sum= float(parameters.target_sum[0])

#create anndata

#Paper: Vasquez-Garcia et al, Ovarian cancer mutational processes drive site-specific immune evasion. Nature 2022 Dec 14. PMID: 
#Counts file: .h5
#%%
adata = sc.read(init_dir + 'GSE180661_matrix.h5ad')
#%%
adata.obs['patient_id'] = adata.obs.index.str.split("_").str[0]
adata.obs['treatment'] = 'Naive'
adata.obs['developmental_stage'] = 'Advanced stage'
adata.obs['recurrence'] = 'Unknown'

#%%
paper_ID =[]
for i in adata.obs.patient_id:
    paper_ID.append("Vasquez_" + i.split('-')[2])

adata.obs["paper_ID"] = paper_ID

adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]
 #%%
tissue = []

for i in adata.obs.index:
    tissue.append(i.rsplit('_', 1)[0].rsplit('_')[-1])

adata.obs["tissue_type"] = tissue

#%%
columns = {}

columns['tissue'] = {'ADNEXA':'Primary', 'TUBE':'Primary', 'OVARY':'Primary', 'OMENTUM':'Metastasis', 
                     'BLADDER':'Metastasis', 'QUADRANT':'Metastasis', 'CECUM':'Metastasis', 'TUMOUR':'Metastasis',
                     'SURFACE':'Metastasis', 'PERITONEUM':'Metastasis', 'BOWEL':'Metastasis',  'MESENTARY':'Metastasis',
                     'DIAPHRAGM':'Metastasis', 'IMPLANT':'Metastasis', 'NODULE':'Metastasis', 'NODE':'Metastasis', 
                     'GUTTER':'Metastasis', 'WALL':'Metastasis', 'ASCITES':'Ascites', 'PELVIS':'Metastasis'} 

columns['anatomical_location'] = {'ADNEXA':'Ovary', 'TUBE':'Ovary', 'OVARY':'Ovary', 'OMENTUM':'Omentum', 
                     'BLADDER':'Urinary_bladder', 'QUADRANT':'Upper_quadrant', 'CECUM':'Other', 'TUMOUR':'Other',
                     'SURFACE':'Other', 'PERITONEUM':'Peritoneum', 'BOWEL':'Bowel',  'MESENTARY':'Mesentery',
                     'DIAPHRAGM':'Other', 'IMPLANT':'Other', 'NODULE':'Other', 'NODE':'Other', 
                     'GUTTER':'Other', 'WALL':'Other', 'ASCITES':'Ascites', 'PELVIS':'Other'}

for column in columns:
    adata.obs[column] = adata.obs.tissue_type.replace(columns[column])

adata.obs = adata.obs.drop(columns=['tissue_type'])
 #%%
columns = {}

columns['tumor_stage'] = {'SPECTRUM-OV-002':'IVB', 'SPECTRUM-OV-003':'IIIC', 'SPECTRUM-OV-007':'IIIC',
       'SPECTRUM-OV-008':'IIIC', 'SPECTRUM-OV-009':'IVB', 'SPECTRUM-OV-014':'IIIC',
       'SPECTRUM-OV-022':'IIIB', 'SPECTRUM-OV-024':'IIIC', 'SPECTRUM-OV-025':'IIIC',
       'SPECTRUM-OV-026':'IVB', 'SPECTRUM-OV-031':'IVB', 'SPECTRUM-OV-036':'IIIC',
       'SPECTRUM-OV-037':'IVB', 'SPECTRUM-OV-041':'IIIC', 'SPECTRUM-OV-042':'IIIC',
       'SPECTRUM-OV-045':'IVB', 'SPECTRUM-OV-049':'IVB', 'SPECTRUM-OV-050':'IIIC',
       'SPECTRUM-OV-051':'IIIC', 'SPECTRUM-OV-052':'IVB', 'SPECTRUM-OV-053':'IIIC',
       'SPECTRUM-OV-054':'IVA', 'SPECTRUM-OV-065':'IIIC', 'SPECTRUM-OV-067':'IVB',
       'SPECTRUM-OV-068':'IIIC', 'SPECTRUM-OV-070':'IIIC', 'SPECTRUM-OV-071':'IVB',
       'SPECTRUM-OV-075':'IVB', 'SPECTRUM-OV-077':'IVB', 'SPECTRUM-OV-080':'IVB',
       'SPECTRUM-OV-081':'IVB', 'SPECTRUM-OV-082':'IVB', 'SPECTRUM-OV-083':'IVB',
       'SPECTRUM-OV-090':'IIIC', 'SPECTRUM-OV-105':'IVB', 'SPECTRUM-OV-107':'IVA',
       'SPECTRUM-OV-110':'IVB', 'SPECTRUM-OV-112':'IVB', 'SPECTRUM-OV-115':'IVB',
       'SPECTRUM-OV-116':'IIIC', 'SPECTRUM-OV-118':'IVB'}

for column in columns:
    adata.obs[column] = adata.obs.patient_id.replace(columns[column])

#%%
adata.obs = adata.obs[['patient_id', 'tissue',  'developmental_stage', 'treatment', 
                           'recurrence', 'tumor_stage', 'paper_ID', 'anatomical_location', 'dataset']]

# %%
adata.write_h5ad(out_dir + 'vasquez2022_rawcounts.h5ad')

#Preprocessing

final_dir = parameters.final_dir[0]

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(final_dir + "vasquez2022_filt_norm_nolog.h5ad")
