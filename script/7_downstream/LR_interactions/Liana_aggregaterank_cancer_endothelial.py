# Running Liana on cancer and endothelial cells divided by treatment type and tissue type

#%%
import scanpy as sc
# Only needed for visualization:
from scanpy.pl import umap
import liana as li
import pandas as pd
import plotnine as p9

## Inizialize directories
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/LR_interactions/cancer_endothelial/'
sc.settings.figdir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/LR_interactions/cancer_endothelial/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled_with_ontologies.h5ad')
adata
adata.obs
adata = adata[adata.obs['ontologies_from_seacells'] != 'nan']
adata = adata[(adata.obs['max'] == 'CancerMSK') | (adata.obs['max'] == 'EndothelialMSK')]

### Create a new column based on conditions
adata.obs['ontologies_from_seacells'] = adata.obs['ontologies_from_seacells'].astype(str)
adata.obs['liana_categories'] = adata.obs['ontologies_from_seacells']
adata.obs.loc[adata.obs['max'] == "CancerMSK", 'liana_categories'] = 'C_' + adata.obs['ontologies_from_seacells']
adata.obs.loc[adata.obs['max'] == "EndothelialMSK", 'liana_categories'] = 'E_' + adata.obs['ontologies_from_seacells']
adata.obs['liana_categories'] = adata.obs['liana_categories'].astype('category')


#%%
adata_cht = adata[(adata.obs['treatment'] == 'CHT')]
adata_naive = adata[(adata.obs['treatment'] == 'Naive')]
adata_nact = adata[(adata.obs['treatment'] == 'NACT')]

adata_primary = adata[(adata.obs['tissue'] == 'Primary')]
adata_metastasis = adata[(adata.obs['tissue'] == 'Metastasis')]
adata_ascites = adata[(adata.obs['tissue'] == 'Ascites')]

#%%
## Methods
li.mt.show_methods()

#%%
## import liana's rank_aggregate
from liana.mt import rank_aggregate
## resources
li.resource.show_resources()

#%%
## import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

#%%
# Run rank_aggregate
li.mt.rank_aggregate(adata_cht, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_naive, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_nact, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)

li.mt.rank_aggregate(adata_primary, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_metastasis, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_ascites, groupby='liana_categories', expr_prop=0.1, verbose=True, use_raw=False)

#%%
adata_cht.uns['liana_res'].head()
adata_naive.uns['liana_res'].head()
adata_nact.uns['liana_res'].head()

rank_aggregate.describe()

## Save the results
#%%
adata_cht.write_h5ad(outDir + 'liana_aggregaterank_cht_labels.h5ad')
adata_naive.write_h5ad(outDir + 'liana_aggregaterank_naive_labels.h5ad')
adata_nact.write_h5ad(outDir + 'liana_aggregaterank_nact_labels.h5ad')

adata_primary.write_h5ad(outDir + 'liana_aggregaterank_primary_labels.h5ad')
adata_metastasis.write_h5ad(outDir + 'liana_aggregaterank_metastasis_labels.h5ad')
adata_ascites.write_h5ad(outDir + 'liana_aggregaterank_ascites_labels.h5ad')

