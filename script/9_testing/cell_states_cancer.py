
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

#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'metacells_step0/endothelial/'
outDir = rawPath + 'integration/metacells/endothelial_testing/'
ooseDir = rawPath + 'out_of_sample_extension/endothelial_testing/'
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_endothelial.csv'

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
#%%
adata = ooseDir + 'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates.h5ad'
adata = sc.read_h5ad(adata)

#%%
adata_ref = adata[~adata.obs_names.str.startswith("new")]
adata_query = adata[adata.obs_names.str.startswith("new")]

counts = adata_ref.obs.cell_type.value_counts()
adata_ref_cell_type = [cell_type for cell_type in counts.index if counts[cell_type] > 10]
adata_ref = adata_ref[adata_ref.obs.cell_type.isin(adata_ref_cell_type)]

counts = adata_query.obs.predicted_cell_types.value_counts()
adata_query_cell_type = [cell_type for cell_type in counts.index if counts[cell_type] > 10]
adata_query = adata_query[adata_query.obs.predicted_cell_types.isin(adata_query_cell_type)]
#%% New Strategy
sc.tl.rank_genes_groups(adata_ref, groupby="cell_type", method="wilcoxon")
sc.tl.rank_genes_groups(adata_query, groupby="predicted_cell_types", method="wilcoxon")

## 2nd strategy
#%%
def get_degs(adata, groupby):  
    for cluster in adata.obs[groupby].unique():

        deg_df = sc.get.rank_genes_groups_df(adata, group=cluster)  # Replace 'wilcoxon' with your actual key if needed

        # Apply your thresholds
        deg_df = deg_df.dropna(subset=['names'])
        deg_df = deg_df[
            (deg_df['logfoldchanges'] > 1) & 
            (deg_df['logfoldchanges'] < 100) & 
            (deg_df['scores'] > 5) & 
            (deg_df['pvals_adj'] < 0.05)
        ]

        # Use the filtered gene names for enrichment
        dc_cluster_genes = deg_df['names'].tolist()

        if len(dc_cluster_genes) > 0:
            gp = GProfiler(return_dataframe=True)
            enrichment_results = gp.profile(
            organism='hsapiens', 
            query=dc_cluster_genes,
            no_evidences=False, 
            background=adata.var_names.to_list(),
            sources=['GO:CC', 'GO:BP', 'GO:MF', 'REAC', 'KEGG']
        )

            df = enrichment_results
            df["-log10(p-value)"] = -np.log10(df["p_value"])

            plt.figure(figsize=(10, 6))
            sns.scatterplot(
                data=df.head(20),
                x="precision",
                y="name",  
                size="intersection_size",
                hue="-log10(p-value)",
                sizes=(20, 300),
                palette="coolwarm",
                edgecolor="black"
            )
            plt.xlabel("Precision")
            plt.ylabel("Enriched Term")
            plt.title(f"Bubble Plot of Enriched Pathways - Cluster {cluster}")
            plt.legend(title="-log10(p-value)", bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.grid(True)
            plt.show()
        else:
            print(f"No significant genes found for cluster {cluster}. Skipping enrichment.")
    return (deg_df, df)

#%%
adata.obs

# %%
adata.obs['zheng_ontologies'] = np.nan

col = {"Immunoreactive_cells":"Immunoreactive_cells", 
        "Cellular_metabolism-extracellular_signaling":"Immunoreactive_cells", 
        "Unknown_primary":"Immunoreactive_cells", 
        "Organelles_organization-metabolism":"RNA_metabolism", 
        "Cycling_cells":"Cycling_cells", 
        "Response_to_extracellular_signals":"Response_to_extracellular_signals", 
        "Unknown_ascites":"Unknown_ascites",
        "Response_to_stress":"Unknown", 
        "Cellular_metabolism":"Cellular_metabolism", 
        "ECM_shaping_cells":"ECM_shaping_cells", 
        "Unknown_metastasis":"Unknown_metastasis", 
        "Organelles_organization-cell_cycle":"Organelles_organization-cell_cycle", 
        "Ciliated_cancer_cells":"Ciliated_cancer_cells",
        "Extracellular_signaling-immune_cells":"Unknown",
        "Extracellular_signaling":"Unknown",
        "Organelles_organization-cell_movement":"Organelles_organization-cell_movement",
        "RNA_metabolism":"RNA_metabolism",
        "INF_mediated_signaling":"INF_mediated_signaling"}

adata.obs['zheng_ontologies'] = adata.obs.cell_states.replace(col)

#%%
adata.obs
set(adata.obs['zheng_ontologies'])

#%%
### compute percentage of matching ontolgies between cell_states and zheng_ontologies
matches = adata.obs['cell_states'].astype(str) == adata.obs['zheng_ontologies'].astype(str)
percentage_matching = matches.mean() * 100
print(f"Percentage of matching ontologies: {percentage_matching:.2f}%")

#%%
adata.write(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat_predicted_and_real_cellstates.h5ad')
# %%
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
        dfs[c + "_up"] = dfs[c][
            (dfs[c]['logfoldchanges'] > 2) & 
            (dfs[c]['scores'] > 5) & 
            (dfs[c]['pvals_adj'] < 0.05)
        ]       
        dfs[c + "_down"] = dfs[c][
            (dfs[c]['logfoldchanges'] < -0.1) & 
            (dfs[c]['scores'] > 5) & 
            (dfs[c]['pvals_adj'] < 0.05)
        ]

    return dfs
# %%
ranks_query = extract_degs(adata_query)
ranks_ref = extract_degs(adata_ref)
both = ranks_query.keys() & ranks_ref.keys()
both = [b for b in both if not b.endswith("_up") and not b.endswith("_down")]
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
    query_up = wrap_gprof(ranks_query ,cluster + "_up", adata_query.var_names.to_list())
    ref_up = wrap_gprof(ranks_ref, cluster + "_up", adata_ref.var_names.to_list())
    query_down = wrap_gprof(ranks_query, cluster + "_down", adata_ref.var_names.to_list())
    ref_down = wrap_gprof(ranks_ref, cluster + "_down", adata_query.var_names.to_list())
    enrichment[cluster + "_up"] = query_up, ref_up
    enrichment[cluster + "_down"] = query_down, ref_down

# %%
