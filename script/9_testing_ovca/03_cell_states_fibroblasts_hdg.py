
## Imports
#%%
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams
import seaborn as sns
from gprofiler import GProfiler
import sys


import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")
major_cell_types = ["fibroblasts", "endothelial", "cancer"]

major_cell_type = "endothelial"
initDir = rawPath + f'metacells_step0/{major_cell_type}/'
outDir = rawPath + f'integration/metacells/{major_cell_type}_testing_ovca/'
ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca/'
genes_path = scriptsPath + f'4_hdg/Tables/atlas_hdg_dispersion_patients_{major_cell_type}.csv'

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")


tissues = ["metastasis", "ascites", "primary"]

adata = ooseDir + f'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates_{major_cell_type}_raw.h5ad'
adata = sc.read_h5ad(adata)

# ad_raw_ref = sc.read_h5ad("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration_backup/integration/metacells/fibroblasts/seacells_hdg_patients_batch_corr_scgen_celltypes_HDG.h5ad")
# ad_raw_query = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/raw_data/integration/metacells/fibroblasts/seacells_hdg_patients_batch_corr_scgen_celltypes_HDG.h5ad")

# New Strategy
adata.obs["02_tissue"] = adata.obs["02_tissue"].str.lower()

adata_ref = adata[~adata.obs_names.str.startswith("new")]
# adata_ref.raw = ad_raw_ref.raw.to_adata()
adata_query = adata[adata.obs_names.str.startswith("new")]
# ad_raw_query.obs_names = ["new_" + name for name in ad_raw_query.obs_names]
# adata_query.raw = ad_raw_query.raw.to_adata()
def adata_by_tissue(adata):
    adata_by_tissue = {}
    for tissue in adata.obs["02_tissue"].unique():
        if sum(adata.obs["02_tissue"] == tissue) > 10:
            adata_by_tissue[tissue] = adata[adata.obs["02_tissue"] == tissue]
    return adata_by_tissue

adata_ref_by_tissue = adata_by_tissue(adata_ref)
adata_query_by_tissue = adata_by_tissue(adata_query)

#%%

adata_ref = adata_ref_by_tissue["metastasis"]
adata_query = adata_query_by_tissue["metastasis"]
#%%
sc.pl.umap(adata_query, color=["predicted_cell_states"])
sc.pl.umap(adata_ref, color=["cell_states"])

#%%
sc.tl.rank_genes_groups(adata_ref, groupby="cell_states", method="wilcoxon", use_raw=True)
sc.tl.rank_genes_groups(adata_query, groupby="predicted_cell_states", method="wilcoxon",use_raw=True)
ref_deg = sc.get.rank_genes_groups_df(adata_ref, group=None)
#%%

sc.tl.rank_genes_groups(adata_ref, groupby="cell_states", method="wilcoxon", use_raw=False)
sc.tl.rank_genes_groups(adata_query, groupby="predicted_cell_states", method="wilcoxon", use_raw=False)
ref_deg = sc.get.rank_genes_groups_df(adata_ref, group=None)
query_deg = sc.get.rank_genes_groups_df(adata_query, group=None)
#%%


def run_gprof(query, background):
    gp = GProfiler(return_dataframe=True)
    enrichment_results = gp.profile(
        organism='hsapiens', 
        query=query,
        no_evidences=False, 
        background=background,
        sources=['GO:CC', 'GO:BP', 'GO:MF', 'REAC', 'KEGG'])
    return enrichment_results

def wrap_gprof(ranks, cluster, background):
    if(ranks[cluster].empty):
        return pd.DataFrame(columns=["source", "native", "name", "p_value", "significant", "description", "term_size", "query_size", "intersection_size", "effective_domain_size", "precision", "recall", "query", "parents", "intersections", "evidences"])
    return run_gprof(ranks[cluster]["names"].to_list(), background)


#%%
def filter_degs(df):
    return df[
        ((df['logfoldchanges'] > 1) | (df['logfoldchanges'] < -1)) &
        (df['pvals_adj'] < 0.05)
    ]
ref_deg = filter_degs(ref_deg)
query_deg = filter_degs(query_deg)
#%%
enrichment = {}
both = set(ref_deg.group) & set(query_deg.group)
for cluster in both:
    query = run_gprof(query_deg[query_deg.group == cluster].names.to_list() ,adata_query.var_names.to_list())
    ref = run_gprof(ref_deg[ref_deg.group == cluster].names.to_list() ,adata_ref.var_names.to_list())
    enrichment[cluster] = query, ref


#%%
import pandas as pd

# Initialize a list to store the results
results = []

# Iterate over the enrichment dictionary
for cluster, (query_df, ref_df) in enrichment.items():
    # Get the sets of names from both dataframes
    query_names = set(query_df['name'])
    ref_names = set(ref_df['name'])
    
    # Calculate the number of overlapping names
    overlap = len(query_names & ref_names)
    
    # Calculate the total number of names in each dataframe
    total_query = len(query_names)
    total_ref = len(ref_names)
    if(total_query == 0 or total_ref == 0):
        perc = 0
    else:
        perc = overlap/min(total_query, total_ref)*100
    
    cluster_tissue = cluster.split("_")[0]
    cluster_cell_state = "_".join(cluster.split("_")[1:])

    no_cells_query = sum(adata_query_by_tissue[cluster_tissue].obs.predicted_cell_states == cluster_cell_state)
    no_cells_ref = sum(adata_ref_by_tissue[cluster_tissue].obs.cell_states == cluster_cell_state)

    # Append the results to the list
    results.append({
        'cluster': cluster,
        'overlap': overlap,
        'total_query': total_query,
        'total_ref': total_ref,
        'no_cells_query': no_cells_query,
        'no_cells_ref': no_cells_ref,
        'percentage': perc
    })

# Convert the results list to a dataframe
overlap_df = pd.DataFrame(results)

# Display the dataframe
overlap_df.to_csv(f"enrichment/{major_cell_type}_overlap_hdg_full.csv")

# %%
