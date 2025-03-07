
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
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from plotting_bubble import scale_data_5_75, plot_enrich 
from ontologies import annotate_ontolgies

#%%
ooseDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/out_of_sample_extension_backup/out_of_sample_extension/cancer/'

#%%
adata = sc.read(ooseDir + 'integrated_query_zheng_seacells_scarches_tissuetreat_predicted_cellstates.h5ad')
adata

#%%
adata.obs

#%%
#check nan values in the cell_states column
adata.obs.cell_states.isna().sum()

#%%
#drop nan values in the cell_states column
adata = adata[~adata.obs.cell_states.isna()]

#%% New Strategy
sc.tl.rank_genes_groups(adata, groupby="cell_states", method="wilcoxon")

## 2nd strategy
#%%
for cluster in adata.obs.cell_states.unique():

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
