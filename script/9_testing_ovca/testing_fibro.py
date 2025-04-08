#%% 
import os
import scanpy as sc
import pandas as pd
import numpy as np
from gprofiler import GProfiler
import sys
# %%
rank_marta = pd.read_csv("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream_backup/downstream/clustering/fibroblasts/metastasis/leiden-0.41/rank_gene_groups_df_0.csv")
# %%
mask = ((rank_marta['logfoldchanges'] > 1) | (rank_marta['logfoldchanges'] < -1)) & (rank_marta['logfoldchanges']  < 100) & (rank_marta['logfoldchanges'] > -100) & (rank_marta['pvals_adj'] < 0.05)
# %%
adata = sc.read_h5ad("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration_backup/integration/metacells/fibroblasts/seacells_hdg_patients_batch_corr_scgen_celltypes_HDG.h5ad")
# %%
adata_mt = adata[(adata.obs['tissue'] == 'Metastasis')]
sc.tl.pca(adata_mt)
sc.pp.neighbors(adata_mt, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata_mt)

# %%
sc.tl.leiden(adata_mt,resolution = 0.41,key_added="leiden-{}".format(round(0.41,2)))
# %%
sc.tl.rank_genes_groups(adata_mt, groupby="leiden-0.41", method='wilcoxon', key_added = "wilcoxon_0.41", use_raw=False)

# %%
ranks = {}
for cl in adata_mt.obs["leiden-0.41"].unique():
    ranks[cl] = sc.get.rank_genes_groups_df(adata_mt, group=cl, key ='wilcoxon_0.41')

# %%
mask_ori = ((ranks["0"]['logfoldchanges'] > 1) | (ranks["0"]['logfoldchanges'] < -1)) & (ranks["0"]['logfoldchanges']  < 100) & (ranks["0"]['logfoldchanges'] > -100) & (ranks["0"]['pvals_adj'] < 0.05)

# %%
cluster_results_ori_path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream_backup/downstream/clustering/fibroblasts/adata_metastasis_embeddings.h5ad"
cluster_results_ori = sc.read_h5ad(cluster_results_ori_path)
# %%
(cluster_results_ori.obs["leiden-0.41"] == "0").sum()
# %%
(adata_mt.obs["leiden-0.41"] == "0").sum()
# %%
adata_mt.obs["leiden-ori"] = cluster_results_ori.obs["leiden-0.41"]
# %%
(adata_mt.obs["leiden-ori"] == "0").sum()
# %%
sc.tl.rank_genes_groups(adata_mt, groupby="leiden-ori", method='wilcoxon', key_added = "wilcoxon_ori")

# %%
ranks = {}
for cl in adata_mt.obs["leiden-ori"].unique():
    ranks[cl] = sc.get.rank_genes_groups_df(adata_mt, group=cl, key ='wilcoxon_ori')

# %%
mask_ori = ((ranks["0"]['logfoldchanges'] > 1) | (ranks["0"]['logfoldchanges'] < -1)) & (ranks["0"]['logfoldchanges']  < 100) & (ranks["0"]['logfoldchanges'] > -100) & (ranks["0"]['pvals_adj'] < 0.05)
mask_ori.sum()
# %%
# Get indices where leiden_ori and leiden-0.41 are the same
matching_indices = adata_mt.obs.index[adata_mt.obs["leiden-ori"] == adata_mt.obs["leiden-0.41"]]

# Count how many cells have matching clusters
print(f"Number of cells with matching leiden clusters: {len(matching_indices)} out of {adata_mt.shape[0]} total cells")
print(f"Percentage of matching cells: {len(matching_indices) / adata_mt.shape[0] * 100:.2f}%")

# Look at distribution of matching cells by cluster
matching_by_cluster = adata_mt.obs.loc[matching_indices].groupby("leiden-ori").size()
print("\nDistribution of matching cells by cluster:")
print(matching_by_cluster)
# %%

leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    print(i)
    if i > 0.4:
        break
    # sc.tl.leiden(adata_mt,resolution = i,key_added="leiden-{}".format(round(i,2)))
    # leidenTotal.append("leiden-{}".format(round(i,2)))
# %%
sc.tl.leiden(adata_mt,resolution = i,key_added="leiden-ori-{}".format(round(i,2)))
#%%
sc.tl.rank_genes_groups(adata_mt, groupby="leiden-ori-0.41", method='wilcoxon', key_added = "wilcoxon_ori")

# %%
matching_indices = adata_mt.obs.index[adata_mt.obs["leiden-0.41"] == adata_mt.obs["leiden-ori-0.41"]]

# Count how many cells have matching clusters
print(f"Number of cells with matching leiden clusters: {len(matching_indices)} out of {adata_mt.shape[0]} total cells")
print(f"Percentage of matching cells: {len(matching_indices) / adata_mt.shape[0] * 100:.2f}%")

# Count total cells in each leiden-ori cluster
total_by_cluster = adata_mt.obs.groupby("leiden-ori").size()

# Look at distribution of matching cells by cluster
matching_by_cluster = adata_mt.obs.loc[matching_indices].groupby("leiden-ori").size()

# Display both total and matching cells by cluster
print("\nDistribution by cluster:")
for cluster in sorted(total_by_cluster.index):
    total = total_by_cluster[cluster]
    matching = matching_by_cluster.get(cluster, 0)
    print(f"Cluster {cluster}: {matching} matching cells out of {total} total cells ({matching/total*100:.2f}%)")

# %%
sc.tl.leiden(cluster_results_ori,resolution = 0.41,key_added="leiden-new-{}".format(round(0.41,2)))

# %%
# %%
matching_indices = cluster_results_ori.obs.index[cluster_results_ori.obs["leiden-new-0.41"] == cluster_results_ori.obs["leiden-0.41"]]

# Count how many cells have matching clusters
print(f"Number of cells with matching leiden clusters: {len(matching_indices)} out of {cluster_results_ori.shape[0]} total cells")
print(f"Percentage of matching cells: {len(matching_indices) / cluster_results_ori.shape[0] * 100:.2f}%")

# Count total cells in each leiden-ori cluster
total_by_cluster = cluster_results_ori.obs.groupby("leiden-0.41").size()

# Look at distribution of matching cells by cluster
matching_by_cluster = cluster_results_ori.obs.loc[matching_indices].groupby("leiden-0.41").size()

# Display both total and matching cells by cluster
print("\nDistribution by cluster:")
for cluster in sorted(total_by_cluster.index):
    total = total_by_cluster[cluster]
    matching = matching_by_cluster.get(cluster, 0)
    print(f"Cluster {cluster}: {matching} matching cells out of {total} total cells ({matching/total*100:.2f}%)")

# %%
