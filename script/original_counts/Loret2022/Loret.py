
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

initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/"
outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Adata/"
min_genes = int(args[1])
min_cells = int(args[2])
genes_by_counts = int(args[3])
pct_counts_mt = float(args[4])
target_sum= float(args[5])


listpath = ["/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049611_2_N_A_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049612_3_N_PER_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049613_4_N_OM_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049614_5_N_BL_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049615_6_T_OT_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049616_7_T_PER_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049617_8_T_OM_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049618_9_T_A_PT1_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049619_10_N_OT_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049620_11_N_A_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049621_12_N_OM_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049622_13_T_OT_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049623_14_T_OM_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049624_15_T_PER_PT2_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049625_16_N_OT_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049626_17_N_A_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049627_18_N_PER_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049628_19_T_OT_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049629_20_T_A_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049630_21_T_OM_PT3_filtered_gene_bc_matrices_h5.h5",
           "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/Counts/GSM6049631_22_T_PER_PT3_filtered_gene_bc_matrices_h5.h5"]
                 

#create anndata

counter = 1

counteradata = 0

adata = ['']*len(listpath)

for i in listpath:
    print(counteradata + 1) 
    print(i)
    adata[counteradata] = sc.read_10x_h5(i)
    adata[counteradata].obs["ID"] = "sample_" + str(counter)
    adata[counteradata].obs["cell"] = adata[counteradata].obs.index + "-" + adata[counteradata].obs["ID"]
    adata[counteradata].obs.index = adata[counteradata].obs["cell"]
    del adata[counteradata].obs["cell"]
    print(adata[counteradata].X.shape)
    adata[counteradata].var_names_make_unique()
    counteradata = counteradata + 1
    counter = counter + 1

adata = adata[0].concatenate(adata[1:], index_unique=None)


#Paper: Loret N, Vandamme N, De Coninck J, Taminau J et al. Distinct transcriptional programs in ascitic and solid cancer cells induce different responses to chemotherapy in high-grade serous ovarian cancer. Mol Cancer Res 2022 Jun 24. PMID: 35749080

#Counts file: h5 files

# Creating metadata

columns = {}

columns['patient_id'] = {"sample_1":"P1", 
                          "sample_2":"P1", 
                          "sample_3":"P1", 
                          "sample_4":"P1", 
                          "sample_5":"P1", 
                          "sample_6":"P1", 
                          "sample_7":"P1", 
                          "sample_8":"P1", 
                          "sample_9":"P1", 
                          "sample_10":"P2", 
                          "sample_11":"P2", 
                          "sample_12":"P2", 
                          "sample_13":"P2", 
                          "sample_14":"P2", 
                          "sample_15":"P2", 
                          "sample_16":"P3", 
                          "sample_17":"P3", 
                          "sample_18":"P3", 
                          "sample_19":"P3", 
                          "sample_20":"P3", 
                          "sample_21":"P3", 
                          "sample_22":"P3"}

columns['sample_name'] = {"sample_1": "GSM6049610", 
                          "sample_2":"GSM6049611", 
                          "sample_3":"GSM6049612", 
                          "sample_4":"GSM6049613", 
                          "sample_5":"GSM6049614", 
                          "sample_6":"GSM6049615", 
                          "sample_7":"GSM6049616", 
                          "sample_8":"GSM6049617", 
                          "sample_9":"GSM6049618", 
                          "sample_10":"GSM6049619", 
                          "sample_11":"GSM6049620", 
                          "sample_12":"GSM6049621", 
                          "sample_13":"GSM6049622", 
                          "sample_14":"GSM6049623", 
                          "sample_15":"GSM6049624", 
                          "sample_16":"GSM6049625", 
                          "sample_17":"GSM6049626", 
                          "sample_18":"GSM6049627", 
                          "sample_19":"GSM6049628", 
                          "sample_20":"GSM6049629", 
                          "sample_21":"GSM6049630", 
                          "sample_22":"GSM6049631"}

columns['tissue'] = {"sample_1":"Primary", 
                          "sample_2":"Ascites", 
                          "sample_3":"Metastasis", 
                          "sample_4":"Metastasis", 
                          "sample_5":"Metastasis", 
                          "sample_6":"Primary", 
                          "sample_7":"Metastasis", 
                          "sample_8":"Metastasis", 
                          "sample_9":"Ascites", 
                          "sample_10":"Primary", 
                          "sample_11":"Ascites", 
                          "sample_12":"Metastasis", 
                          "sample_13":"Primary", 
                          "sample_14":"Metastasis", 
                          "sample_15":"Metastasis", 
                          "sample_16":"Primary", 
                          "sample_17":"Ascites", 
                          "sample_18":"Metastasis", 
                          "sample_19":"Primary", 
                          "sample_20":"Ascites", 
                          "sample_21":"Metastasis", 
                          "sample_22":"Metastasis"}

columns['developmental_stage'] = {"sample_1":"Advanced stage",                           
                                  "sample_2":"Advanced stage", 
                                  "sample_3":"Advanced stage", 
                                  "sample_4":"Advanced stage", 
                                  "sample_5":"Advanced stage", 
                                  "sample_6":"Advanced stage", 
                                  "sample_7":"Advanced stage", 
                                  "sample_8":"Advanced stage", 
                                  "sample_9":"Advanced stage", 
                                  "sample_10":"Advanced stage", 
                                  "sample_11":"Advanced stage", 
                                  "sample_12":"Advanced stage", 
                                  "sample_13":"Advanced stage", 
                                  "sample_14":"Advanced stage", 
                                  "sample_15":"Advanced stage", 
                                  "sample_16":"Advanced stage", 
                                  "sample_17":"Advanced stage", 
                                  "sample_18":"Advanced stage", 
                                  "sample_19":"Advanced stage", 
                                  "sample_20":"Advanced stage", 
                                  "sample_21":"Advanced stage", 
                                  "sample_22":"Advanced stage"}

columns['treatment'] = {"sample_1":"Naive",                           
                        "sample_2":"Naive", 
                        "sample_3":"Naive", 
                        "sample_4":"Naive", 
                        "sample_5":"Naive", 
                        "sample_6":"NACT", 
                        "sample_7":"NACT", 
                        "sample_8":"NACT", 
                        "sample_9":"NACT", 
                        "sample_10":"Naive", 
                        "sample_11":"Naive", 
                        "sample_12":"Naive", 
                        "sample_13":"NACT", 
                        "sample_14":"NACT", 
                        "sample_15":"NACT", 
                        "sample_16":"Naive", 
                        "sample_17":"Naive", 
                        "sample_18":"Naive", 
                        "sample_19":"NACT", 
                        "sample_20":"NACT", 
                        "sample_21":"NACT", 
                        "sample_22":"NACT"}

columns['recurrence'] = {"sample_1":"Recurrence",                           
                        "sample_2":"Recurrence", 
                        "sample_3":"Recurrence", 
                        "sample_4":"Recurrence", 
                        "sample_5":"Recurrence", 
                        "sample_6": "Recurrence", 
                        "sample_7":"Recurrence", 
                        "sample_8":"Recurrence", 
                        "sample_9":"Recurrence", 
                        "sample_10":"Recurrence", 
                        "sample_11":"Recurrence", 
                        "sample_12":"Recurrence", 
                        "sample_13":"Recurrence", 
                        "sample_14":"Recurrence", 
                        "sample_15":"Recurrence", 
                        "sample_16":"Recurrence", 
                        "sample_17":"Recurrence", 
                        "sample_18":"Recurrence", 
                        "sample_19":"Recurrence", 
                        "sample_20":"Recurrence", 
                        "sample_21":"Recurrence", 
                        "sample_22":"Recurrence"}

columns['tumor_stage'] = {"sample_1":"IIIC", 
                          "sample_2":"IIIC", 
                          "sample_3":"IIIC", 
                          "sample_4":"IIIC", 
                          "sample_5":"IIIC", 
                          "sample_6":"IIIC", 
                          "sample_7":"IIIC", 
                          "sample_8":"IIIC", 
                          "sample_9":"IIIC", 
                          "sample_10":"IIIC", 
                          "sample_11":"IIIC", 
                          "sample_12":"IIIC", 
                          "sample_13":"IIIC", 
                          "sample_14":"IIIC", 
                          "sample_15":"IIIC", 
                          "sample_16":"IIIB", 
                          "sample_17":"IIIB", 
                          "sample_18":"IIIB", 
                          "sample_19":"IIIB", 
                          "sample_20":"IIIB", 
                          "sample_21":"IIIB", 
                          "sample_22":"IIIB"}

columns['paper_ID'] = {"sample_1":"Loret_1", 
                          "sample_2":"Loret_1", 
                          "sample_3":"Loret_1", 
                          "sample_4":"Loret_1", 
                          "sample_5":"Loret_1", 
                          "sample_6":"Loret_1", 
                          "sample_7":"Loret_1", 
                          "sample_8":"Loret_1", 
                          "sample_9":"Loret_1", 
                          "sample_10":"Loret_2", 
                          "sample_11":"Loret_2", 
                          "sample_12":"Loret_2", 
                          "sample_13":"Loret_2", 
                          "sample_14":"Loret_2", 
                          "sample_15":"Loret_2", 
                          "sample_16":"Loret_3", 
                          "sample_17":"Loret_3", 
                          "sample_18":"Loret_3", 
                          "sample_19":"Loret_3", 
                          "sample_20":"Loret_3", 
                          "sample_21":"Loret_3", 
                          "sample_22":"Loret_3"}

columns['anatomical_location'] = {"sample_1":"Ovary",                           
                                  "sample_2":"Ascites", 
                                  "sample_3":"Peritoneum", 
                                  "sample_4":"Omentum", 
                                  "sample_5":"Urinary_bladder", 
                                  "sample_6":"Ovary", 
                                  "sample_7":"Peritoneum", 
                                  "sample_8":"Omentum", 
                                  "sample_9":"Ascites", 
                                  "sample_10":"Ovary", 
                                  "sample_11":"Ascites", 
                                  "sample_12":"Omentum", 
                                  "sample_13":"Ovary", 
                                  "sample_14":"Omentum", 
                                  "sample_15":"Peritoneum", 
                                  "sample_16":"Ovary", 
                                  "sample_17":"Ascites", 
                                  "sample_18":"Peritoneum", 
                                  "sample_19":"Ovary", 
                                  "sample_20":"Ascites", 
                                  "sample_21":"Omentum", 
                                  "sample_22":"Peritoneum"}

for column in columns:
    adata.obs[column] = adata.obs.ID.replace(columns[column])
    

adata.obs.rename(columns = {'sample_ID':'ID'}, inplace = True)
del adata.obs['batch']
del adata.var["gene_ids"]
adata.uns["preprocessing"] = ({"min_genes" : 200, "min_cells" : 3,  "genes_by_counts" : 5000, "pct_counts_mt" : 20, "target_sum" : 1e4})

#Write raw adata with metadata
adata.write_h5ad(outDir + "loret2022_rawcounts.h5ad")

#Preprocessing

finalDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Loret2022/"

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < genes_by_counts, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
sc.pp.normalize_total(adata, target_sum=target_sum)

adata.write(finalDir + "loret2022_filt_norm_nolog.h5ad")