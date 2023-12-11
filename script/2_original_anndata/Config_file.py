# Config file

## HELP:
# To run this script, the argument word takes one of these elements of this list:
# ['Zhang2022',
#  'Ren2022',
#  'Xu2022',
#  'Qian2020',
#  'Regner2021',
#  'Geistlinger2020'
#  'Olbrecht2021',
#  'Loret2022',
#  'Vasquez2022']

# Then, on HPC run the following:
# sbatch  --job-name Zhang2022 --output %x_%jZhang2022.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Zhang2022
# sbatch  --job-name Ren2022 --output %x_%jRen2022.log /home/marta.sallese/ov_cancer_atlasatlas_project/script/original_anndata/parallel_submission.sh Ren2022
# sbatch  --job-name Xu2022 --output %x_%jXu2022.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Xu2022
# sbatch  --job-name Regner2021 --output %x_%jRegner2021.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Regner2021
# sbatch  --job-name Olbrecht2021 --output %x_%jOlbrecht2021.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Olbrecht2021
# sbatch  --job-name Geistlinger2020 --output %x_%jGeistlinger2020.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Geistlinger2020
# sbatch  --job-name Loret2022 --output %x_%jLoret2022.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Loret2022
# sbatch  --job-name Vasquez2022 --output %x_%jVasquez2022.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Vasquez2022
# sbatch  --job-name Qian2020 --output %x_%jQian2020.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Qian2020
# sbatch  --job-name Loret2022 --output %x_%jLoret2022.log /home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/parallel_submission.sh Loret2022

#%%
import os
import re
import argparse
import sys
import subprocess

python_bin = '/home/marta.sallese/ov_cancer_atlas/miniconda3/bin/python3'
args = sys.argv

word = args[1] #see HELP

#%%
def preprocess(adataDir, outPath) :
    import scanpy as sc
    adata = sc.read(adataDir)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    adata.write_h5ad(outPath)

#%%
dirs = [ x for x in list(os.listdir("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/")) if "." not in x]
dir_name = list(filter(lambda x: word in x, dirs))[0]
file_name = dir_name.lower()
script_name = re.sub("[0-9]*","",dir_name)
print(dir_name)
print(file_name)
print(script_name)
adataDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/" + dir_name + "/Adata/" + file_name + "_filt_norm_nolog.h5ad" 
outPath = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/" + dir_name + "/Adata/" + file_name + "_preproccessed.h5ad"

if word == "Zhang2022":
    ## Zhang
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Zhang2022/" 
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Zhang2022/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 5000
    pct_counts_mt = 8
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Zhang2022/Zhang.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

#%%
if word == "Vasquez2022":
    ## Vasquez
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Vasquez2022/"
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Vasquez2022/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 5000
    pct_counts_mt = 20
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Vasquez2022/Vasquez.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Xu2022":
    ## Xu
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Xu2022/"
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Xu2022/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 5000
    pct_counts_mt = 20
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Xu_2022/Xu.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)
    
if word == "Geistlinger2020":
    ## Geistlinger
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Geistlinger2020/"
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Geistlinger2020/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 5000
    pct_counts_mt = 15
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Geistlinger2020/Geistlinger.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Loret2022":
    ## Loret
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Loret2022/"
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Loret2022/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 5000
    pct_counts_mt = 20
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Loret2022/Loret.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Olbrecht2021":
    ## Olbrecht
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Olbrecht2021/"
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Olbrecht2021/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 6000
    pct_counts_mt = 20
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/OLbrecht2021/Olbrecht.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Qian2020":
    ## Qian
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Qian2020/" 
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Qian2020/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 4000
    pct_counts_mt = 12
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Qian2020/Qian.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Regner2021":
    ## Regner
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Regner2021/" 
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Regner2021/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 8000
    pct_counts_mt = 20
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Regner2021/Regner.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

if word == "Ren2022":
    ## Ren
    initDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_counts/Ren2022/" 
    outDir = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/Ren2022/Adata/"
    min_genes = 200
    min_cells = 3
    genes_by_counts = 10000
    pct_counts_mt = 30
    target_sum=1e4
    script_name = "/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_counts/Ren2022/Ren.py"
    subprocess.check_output("{0} {1} {2} {3} {4} {5} {6}".format(python_bin, script_name, min_genes, min_cells, genes_by_counts, pct_counts_mt, target_sum), shell=True)
    preprocess(adataDir, outPath)

