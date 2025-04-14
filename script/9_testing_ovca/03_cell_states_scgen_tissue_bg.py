
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
#%%
# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")
major_cell_types = ["fibroblasts", "endothelial", "cancer"]
ranks = pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'major_celltype', 'tissue', 'is_ref', 'enrichment_group'])
# ovca = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
# ovca_genes = ovca.var_names.to_list()
# ovca_genes_raw = ovca.raw.to_adata().var_names.to_list()
# del ovca
#%%
tissues = ["metastasis", "ascites", "primary"]
adata_ref_by = {}
adata_query_by = {}
ovca_genes_raw = []
for major_cell_type in major_cell_types:

    initDir = rawPath + f'metacells_step0/{major_cell_type}/'
    outDir = rawPath + f'integration/metacells/{major_cell_type}_testing_ovca_corrected/'
    ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca_corrected/'

    utilsPath = config.get("DEFAULT", "utilsPath")
    rawPath = config.get("DEFAULT", "rawPath")
    scriptsPath = config.get("DEFAULT", "scriptsPath")

    #adata = ooseDir + f'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates_{major_cell_type}_raw.h5ad'

    tissues = ["metastasis", "ascites", "primary"]
    adata_ref_by_tissue = {}
    adata_query_by_tissue = {}

    for tissue in tissues:
        adata = ooseDir + f'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates_{major_cell_type}_{tissue}corrected_raw.h5ad'
        adata = sc.read_h5ad(adata)
        adata.obs["02_tissue"] = adata.obs["02_tissue"].str.lower()
        if ovca_genes_raw == []:
            ovca_genes_raw = adata.raw.to_adata().var_names.to_list()

        adata_ref = adata[~adata.obs_names.str.startswith("new")]
        adata_ref.obs["cell_states"] = adata_ref.obs["07_cell_states"]
        adata_query = adata[adata.obs_names.str.startswith("new")]
        adata_query.obs["cell_states"] = adata_query.obs["predicted_cell_states"]

        # adata_ref.raw = ad_raw_ref.raw.to_adata()
        adata_ref_by_tissue[tissue] = adata_ref[adata_ref.obs["02_tissue"] == tissue]
        adata_query_by_tissue[tissue] = adata_query[adata_query.obs["02_tissue"] == tissue]
        

    both = adata_ref_by_tissue.keys() & adata_query_by_tissue.keys()
    
    for tissue in both:
        counts = adata_ref_by_tissue[tissue].obs.cell_states.value_counts()
        adata_ref_cell_states = [cell_states for cell_states in counts.index if counts[cell_states] > 10]
        adata_ref_by_tissue[tissue] = adata_ref_by_tissue[tissue][adata_ref_by_tissue[tissue].obs.cell_states.isin(adata_ref_cell_states)]

        counts = adata_query_by_tissue[tissue].obs.cell_states.value_counts()
        adata_query_cell_states = [cell_states for cell_states in counts.index if counts[cell_states] > 10]
        adata_query_by_tissue[tissue] = adata_query_by_tissue[tissue][adata_query_by_tissue[tissue].obs.cell_states.isin(adata_query_cell_states)]

    for tissue in both:
        sc.tl.rank_genes_groups(adata_ref_by_tissue[tissue], groupby="cell_states", method="wilcoxon", use_raw=False)         # TODO vedere se cambia se filtro dopo
        df_ref = pd.concat([sc.get.rank_genes_groups_df(adata_ref_by_tissue[tissue], group=None, pval_cutoff=0.05, log2fc_max=100, log2fc_min=1),
                            sc.get.rank_genes_groups_df(adata_ref_by_tissue[tissue], group=None, pval_cutoff=0.05, log2fc_max=-1, log2fc_min=-100)])
        df_ref["tissue"] = tissue
        df_ref["major_celltype"] = major_cell_type
        df_ref["is_ref"] = True
        df_ref["enrichment_group"] = [major_cell_type + "_" + tissue + "_" + g for g in df_ref["group"]]
        sc.tl.rank_genes_groups(adata_query_by_tissue[tissue], groupby="cell_states", method="wilcoxon", use_raw=False)
        df_query = pd.concat([sc.get.rank_genes_groups_df(adata_query_by_tissue[tissue], group=None, pval_cutoff=0.05, log2fc_max=100, log2fc_min=1),
                              sc.get.rank_genes_groups_df(adata_query_by_tissue[tissue], group=None, pval_cutoff=0.05, log2fc_max=-1, log2fc_min=-100)])
        df_query["tissue"] = tissue
        df_query["major_celltype"] = major_cell_type
        df_query["is_ref"] = False
        df_query["enrichment_group"] = [major_cell_type + "_" + tissue + "_" + g for g in df_query["group"]] #da rivedere check che non si sputtani
        ranks = pd.concat([ranks ,df_ref, df_query], axis=0)
    adata_ref_by[major_cell_type] = adata_ref_by_tissue
    adata_query_by[major_cell_type] = adata_query_by_tissue

ranks.reset_index(drop=True, inplace=True)
ranks["is_ref"] = ranks["is_ref"].astype(bool)
ranks.to_csv("enrichment/ranks_scgen_tissue_bg.csv", index=False)

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
    g_ref = wrap_gprof(ranks_ref["names"], ovca_genes_raw)
    g_ref["is_ref"] = True
    g_ref["enrichment_group"] = enrichment_group
    g_query = wrap_gprof(ranks_query["names"], ovca_genes_raw)
    g_query["is_ref"] = False
    g_query["enrichment_group"] = enrichment_group
    g_prof = pd.concat([g_prof, g_ref, g_query], ignore_index=True, axis=0)

#%%
g_prof.to_csv("enrichment/enrichment_scgen_tissue_bg.csv")

#%%
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
    major_cell_type = split[0]
    tissue = split[1]
    cell_state = '_'.join(split[2:])

    summary_data.append({
        'enrichment_group': enrichment_group,
        'major_celltype':major_cell_type,
        'tissue':tissue,
        'cell_state':cell_state,
        'shared_term_count': len(shared_terms),
        'ref_only_term_count': len(ref_terms - query_terms), 
        'query_only_term_count': len(query_terms - ref_terms),
        'ref_total_term_count': len(ref_terms),
        'query_total_term_count': len(query_terms),
        'union_terms': len(ref_terms.union(query_terms)),
        'overlap_percentage': len(shared_terms) / len(ref_terms) * 100 if len(ref_terms) > 0 else 0,
        'is_terms_zero': ref_terms.union(query_terms) == 0,
        'is_terms_few': len(ref_terms.union(query_terms)) < 5,
        'ref_cells':sum(adata_ref_by[major_cell_type][tissue].obs["cell_states"] == cell_state),
        'query_cells':sum(adata_query_by[major_cell_type][tissue].obs["cell_states"] == cell_state),
    })

summary_df = pd.DataFrame(summary_data)
summary_df = summary_df.sort_values('overlap_percentage', ascending=False)
#%%
summary_df.to_csv("enrichment/summary_enrichment_scgen_tissue_bg.csv", index=False)

#%%
nonZero = summary_df[(summary_df["ref_cells"] > 0) & (summary_df["query_cells"] > 0)]
sns.boxenplot(data=nonZero, y="overlap_percentage", x="major_celltype", hue="tissue")
plt.ylim(0, 100)  # Set y-axis range from 0 to 100

# %%
