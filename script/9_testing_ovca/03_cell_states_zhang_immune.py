
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
from sklearn.neighbors import KNeighborsClassifier
#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")


integrationDir = rawPath + 'integration/metacells/cancer/'

path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/out_of_sample_extension_backup/out_of_sample_extension/immune/integrated_query_zheng_seacells_scarches_tissuetreat_predicted_cellsubtypes.h5ad"
adata = sc.read_h5ad(path)
#%%
atlas_adata = adata[adata.obs["reference_map"] == "reference"]
new_dataset_adata = adata[adata.obs["reference_map"] != "reference"]

#%%
X_train = atlas_adata.obsm["latent_corrected"]
y_train = atlas_adata.obs['cell_subtypes'].to_numpy()

X_test = new_dataset_adata.obsm["latent_corrected"]

#%%
# Train kNN classifier
knn = KNeighborsClassifier(n_neighbors=5, metric="cosine")  # Use cosine similarity for high-dimensional data
knn.fit(X_train, y_train)

#%%
# Predict labels for new dataset
new_dataset_adata.obs["predicted_cell_subtypes"] = knn.predict(X_test)
new_dataset_adata.obs["predicted_cell_subtypes"] = new_dataset_adata.obs["predicted_cell_subtypes"].astype("category")
#%%
all_categories = pd.api.types.union_categoricals([
    adata.obs['cell_subtypes'],
    new_dataset_adata.obs["predicted_cell_subtypes"]
    ]).categories

# Align the categories of 'ascites_state' and 'metastasis_state' to match the combined categories
adata.obs['cell_subtypes'] = adata.obs['cell_subtypes'].cat.set_categories(all_categories)
new_dataset_adata.obs["predicted_cell_subtypes"] = new_dataset_adata.obs["predicted_cell_subtypes"].cat.set_categories(all_categories)
adata.obs.loc[new_dataset_adata.obs_names, "cell_subtypes"] = new_dataset_adata.obs["predicted_cell_subtypes"]

#%%

adata_ref = adata[~(adata.obs["dataset"] == "Zheng")]
adata_query = adata[adata.obs["dataset"] == "Zheng"]
adata_query.obs["predicted_cell_subtypes"] = adata_query.obs["cell_subtypes"]
#%% New Strategy
def adata_by_tissue(adata):
    adata_by_tissue = {}
    for tissue in adata.obs["tissue"].unique():
        if sum(adata.obs["tissue"] == tissue) > 10:
            adata_by_tissue[tissue] = adata[adata.obs["tissue"] == tissue]
    return adata_by_tissue

adata_ref_by_tissue = adata_by_tissue(adata_ref)
adata_query_by_tissue = adata_by_tissue(adata_query)


#%%
both = adata_ref_by_tissue.keys() & adata_query_by_tissue.keys()

for tissue in both:
    counts = adata_ref_by_tissue[tissue].obs.cell_subtypes.value_counts()
    adata_ref_cell_subtypes = [cell_subtypes for cell_subtypes in counts.index if counts[cell_subtypes] > 10]
    adata_ref_by_tissue[tissue] = adata_ref_by_tissue[tissue][adata_ref_by_tissue[tissue].obs.cell_subtypes.isin(adata_ref_cell_subtypes)]

    counts = adata_query_by_tissue[tissue].obs.predicted_cell_subtypes.value_counts()
    adata_query_cell_subtypes = [cell_subtypes for cell_subtypes in counts.index if counts[cell_subtypes] > 10]
    adata_query_by_tissue[tissue] = adata_query_by_tissue[tissue][adata_query_by_tissue[tissue].obs.predicted_cell_subtypes.isin(adata_query_cell_subtypes)]

for tissue in both:
    sc.tl.rank_genes_groups(adata_ref_by_tissue[tissue], groupby="cell_subtypes", method="wilcoxon")

#%%
ranks_query = {}
for tissue in both:
    ranks_query[tissue] = {}
    for cell_state in adata_query_by_tissue[tissue].obs.cell_subtypes.cat.categories:
        adata_ref_test = adata_ref_by_tissue[tissue][adata_ref_by_tissue[tissue].obs.cell_subtypes != cell_state]
        adata_test = sc.concat([adata_ref_test, adata_query_by_tissue[tissue][adata_query_by_tissue[tissue].obs.cell_subtypes == cell_state]])
        sc.tl.rank_genes_groups(adata_test, groupby="cell_subtypes", method="wilcoxon", groups=[cell_state])
        ranks = adata_test.uns["rank_genes_groups"]
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
        ranks_query[tissue][cell_state] = dfs[cell_state]

## 2nd strategy

#%%
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
# %%
#ranks_query = {tissue: extract_degs(adata_query_by_tissue[tissue]) for tissue in both}
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

#%%
both = ranks_query.keys() & ranks_ref.keys()
#both = [b for b in both if not b.endswith("_up") and not b.endswith("_down")]
# %%
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
# %%
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

    no_cells_query = sum(adata_query_by_tissue[cluster_tissue].obs.predicted_cell_subtypes == cluster_cell_state)
    no_cells_ref = sum(adata_ref_by_tissue[cluster_tissue].obs.cell_subtypes == cluster_cell_state)

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
overlap_df.to_csv("zheng_immune_overlap_filter.csv")
# %%
