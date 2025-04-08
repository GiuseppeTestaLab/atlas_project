
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
for major_cell_type in major_cell_types:
    initDir = rawPath + f'metacells_step0/{major_cell_type}/'
    outDir = rawPath + f'integration/metacells/{major_cell_type}_testing/'
    ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing/'
    genes_path = scriptsPath + f'4_hdg/Tables/atlas_hdg_dispersion_patients_{major_cell_type}.csv'

    utilsPath = config.get("DEFAULT", "utilsPath")
    rawPath = config.get("DEFAULT", "rawPath")
    scriptsPath = config.get("DEFAULT", "scriptsPath")
    
    adata = ooseDir + 'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates.h5ad'
    adata = sc.read_h5ad(adata)

    # ad_raw_ref = sc.read_h5ad(f"/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration_backup/integration/metacells/{major_cell_type}/seacells_hdg_patients_batch_corr_scgen_celltypes_HDG.h5ad")
    # ad_raw_query = sc.read_h5ad(f"/group/testa/Project/OvarianAtlasTestStep0/raw_data/integration/metacells/{major_cell_type}/seacells_hdg_patients_batch_corr_scgen_celltypes_HDG.h5ad")


    tissues = ["metastasis", "ascites", "primary"]

    
    adata_ref = adata[~adata.obs_names.str.startswith("new")]
    # adata_ref.raw = ad_raw_ref.raw.to_adata()
    adata_query = adata[adata.obs_names.str.startswith("new")]
    # ad_raw_query.obs_names = ["new_" + name for name in ad_raw_query.obs_names]
    # adata_query.raw = ad_raw_query.raw.to_adata()

    # New Strategy
    def adata_by_tissue(adata):
        adata_by_tissue = {}
        for tissue in adata.obs["tissue"].unique():
            if sum(adata.obs["tissue"] == tissue) > 10:
                adata_by_tissue[tissue] = adata[adata.obs["tissue"] == tissue]
        return adata_by_tissue

    adata_ref_by_tissue = adata_by_tissue(adata_ref)
    adata_query_by_tissue = adata_by_tissue(adata_query)


    
    both = adata_ref_by_tissue.keys() & adata_query_by_tissue.keys()

    for tissue in both:
        counts = adata_ref_by_tissue[tissue].obs.cell_states.value_counts()
        adata_ref_cell_states = [cell_states for cell_states in counts.index if counts[cell_states] > 10]
        adata_ref_by_tissue[tissue] = adata_ref_by_tissue[tissue][adata_ref_by_tissue[tissue].obs.cell_states.isin(adata_ref_cell_states)]

        counts = adata_query_by_tissue[tissue].obs.predicted_cell_states.value_counts()
        adata_query_cell_states = [cell_states for cell_states in counts.index if counts[cell_states] > 10]
        adata_query_by_tissue[tissue] = adata_query_by_tissue[tissue][adata_query_by_tissue[tissue].obs.predicted_cell_states.isin(adata_query_cell_states)]

    for tissue in both:
        sc.tl.rank_genes_groups(adata_ref_by_tissue[tissue], groupby="cell_states", method="wilcoxon", use_raw=False)
        sc.tl.rank_genes_groups(adata_query_by_tissue[tissue], groupby="predicted_cell_states", method="wilcoxon", use_raw=False)
    ## 2nd strategy

    
    def extract_degs(adata):
        ranks = adata.uns["rank_genes_groups"]
        pvals = pd.DataFrame(ranks["pvals_adj"])
        names = pd.DataFrame(ranks["names"])
        scores = pd.DataFrame(ranks["scores"])
        change = pd.DataFrame(ranks["logfoldchanges"])
        dfs = {}
        for c in change.columns:
            dfs[c] = pd.concat([names[c], pvals[c], scores[c], change[c]], axis=1)
            dfs[c].columns = ["names", "pvals_adj", "scores", "logfoldchanges"]
            dfs[c] = dfs[c][
                ((dfs[c]['logfoldchanges'] > 1) | (dfs[c]['logfoldchanges'] < -1)) &
                (dfs[c]['logfoldchanges']  < 100) &
                (dfs[c]['logfoldchanges'] > -100) &            
                (dfs[c]['pvals_adj'] < 0.05)
            ]
        return dfs
    
    ranks_query = {tissue: extract_degs(adata_query_by_tissue[tissue]) for tissue in both}
    ranks_ref = {tissue: extract_degs(adata_ref_by_tissue[tissue]) for tissue in both}
    from collections.abc import MutableMapping

    def flatten(dictionary, parent_key='', separator='_'):
        items = []
        for key, value in dictionary.items():
            new_key = parent_key + separator + key if parent_key else key
            if isinstance(value, MutableMapping):
                items.extend(flatten(value, new_key, separator=separator).items())
            else:
                items.append((new_key, value))
        return dict(items)
    ranks_query = flatten(ranks_query)
    ranks_ref = flatten(ranks_ref)
    def write_degs(degs, path):
        for cluster in degs:
            degs[cluster].to_csv(f"{path}_{cluster}_de_sea_type.csv")
    os.makedirs(f"enrichment/degs/{major_cell_type}", exist_ok=True)
    write_degs(ranks_query, f"enrichment/degs/{major_cell_type}/degs")
    write_degs(ranks_ref, f"enrichment/degs/{major_cell_type}/degs")

    
    both = ranks_query.keys() & ranks_ref.keys()
    #both = [b for b in both if not b.endswith("_up") and not b.endswith("_down")]
    
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

    enrichment = {}
    for cluster in both:
        query = wrap_gprof(ranks_query ,cluster, adata_query.var_names.to_list())
        ref = wrap_gprof(ranks_ref, cluster, adata_ref.var_names.to_list())
        enrichment[cluster] = query, ref
    
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
    overlap_df.to_csv(f"enrichment/{major_cell_type}_overlap_sea_type.csv")
