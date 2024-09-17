import configparser
import scanpy as sc
import pandas as pd

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

# Get datasets and initial path from the configuration file
datasets = config.get("DATASETS", "datasets").split(", ")
rawPath = config.get("DEFAULT", "rawPath")

# %%
var_names = {}
common_var_names = []
for dataset in datasets:
    path = rawPath + "original_anndata/" + dataset
    adata = sc.read(path + "/" + dataset + "_filt_norm_nolog.h5ad")
    common_var_names = common_var_names & adata.var_names

pd.DataFrame(index=common_var_names).to_csv(
    rawPath + "original_anndata/common_varnames_datasets.csv"
)
