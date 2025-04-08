
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
import json

import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")
major_cell_types = ["fibroblasts", "endothelial", "cancer"]
ovca = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
ovca_genes = ovca.var_names.to_list()
del ovca
ranks = pd.read_csv(f"enrichment/ranks_by_tissue.csv")
#%%
gp = GProfiler(return_dataframe=True)
def run_gprof(query, background):
    enrichment_results = gp.profile(
        organism='hsapiens', 
        query=query,
        no_evidences=False, 
        background=background,
        sources=['GO:CC', 'GO:BP', 'GO:MF', 'REAC', 'KEGG'])
    return enrichment_results

def wrap_gprof(ranks, background):
    if(ranks.empty):
        return pd.DataFrame(columns=["source", "native", "name", "p_value", "significant", "description", "term_size", "query_size", "intersection_size", "effective_domain_size", "precision", "recall", "query", "parents", "intersections", "evidences"])
    return run_gprof(ranks.to_list(), background)

#%%
g_prof = pd.DataFrame()

for enrichment_group in ranks.enrichment_group.unique():
    ranks_ref = ranks[(ranks.is_ref) & 
                        (ranks.enrichment_group == enrichment_group)]
    ranks_query = ranks[(~ranks.is_ref) &
                        (ranks.enrichment_group == enrichment_group)]
    g_ref = wrap_gprof(ranks_ref["names"], ovca_genes)
    g_ref["is_ref"] = True
    g_ref["enrichment_group"] = enrichment_group
    g_query = wrap_gprof(ranks_query["names"], ovca_genes)
    g_query["is_ref"] = False
    g_query["enrichment_group"] = enrichment_group
    g_prof = pd.concat([g_prof, g_ref, g_query], ignore_index=True, axis=0)

#%%
g_prof.to_csv("enrichment/enrichment.csv")
# %%
summary_data = []

for enrichment_group in g_prof.enrichment_group.unique():
    # Get unique term values for each is_ref value
    ref_terms = set(g_prof[(g_prof.enrichment_group == enrichment_group) & 
                           (g_prof.is_ref == True)]['native'].unique())
    
    query_terms = set(g_prof[(g_prof.enrichment_group == enrichment_group) & 
                             (g_prof.is_ref == False)]['native'].unique())
    
    # Calculate intersection
    shared_terms = ref_terms.intersection(query_terms)
    
    split = enrichment_group.split("_")
    summary_data.append({
        'enrichment_group': enrichment_group,
        'major_celltype':split[0],
        'tissue':split[1],
        'cell_state':'_'.join(split[2:]),
        'shared_term_count': len(shared_terms),
        'ref_only_term_count': len(ref_terms - query_terms),
        'query_only_term_count': len(query_terms - ref_terms),
        'ref_total_term_count': len(ref_terms),
        'query_total_term_count': len(query_terms),
        'union_terms': len(ref_terms.union(query_terms)),
        'overlap_percentage': len(shared_terms) / (len(ref_terms.union(query_terms))) * 100 if len(ref_terms.union(query_terms)) > 0 else 0,
        'is_terms_zero': ref_terms.union(query_terms) == 0,
        'is_terms_few': len(ref_terms.union(query_terms)) < 5,
        'ref_cells':,
        'query_cells':
    })

summary_df = pd.DataFrame(summary_data)
summary_df = summary_df.sort_values('overlap_percentage', ascending=False)


# %%
