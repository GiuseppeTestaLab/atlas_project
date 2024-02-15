# Appending cancer clusters derived from metacells to the atlas of cells

#%%
## Import libraries
import scanpy as sc
import pandas as pd
import numpy as np

## Initialize directories
initDir = '/group/testa/Project/OvariaAtlas/atlas_project/raw_data/atlas_annotated/'
clustersDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/cancer/cluster_assignments/'

#%%
## Load the data
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled.h5ad')
clusters = sc.read(clustersDir + 'atlas_cancer_clusters_from_seacells.h5ad')
#%%
adata
adata.obs
clusters
clusters.obs

#%%
##
adata.obs = adata.obs.join(clusters.obs['cluster_from_seacells'])
adata.obs

adata.obs['cell_type'] = np.where(adata.obs['cell_type'] == "CancerMSK", adata.obs['cluster_from_seacells'], adata.obs['cell_type'])
adata.obs

adata.obs.drop(columns=['cluster_from_seacells'], inplace=True)
adata.obs

adata.obs.cell_type.unique()

#%%
## Creating a new columns with major cell types

values = []

for index, row in adata.obs.iterrows():
    if row['cell_type'] == 'T_CD4_naive':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD4_CXCL13':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD4_reg':
        values.append('T_CD4_cells')
    elif row['cell_type'] == 'T_CD8_cytotoxic':
        values.append('T_CD8_cells')
    elif row['cell_type'] == 'T_CD8_CXCL13':
        values.append('T_CD8_cells')
    elif row['cell_type'] == 'NK_CD56':
        values.append('NK_cells')
    elif row['cell_type'] == 'NK_cytotoxic':
        values.append('NK_cells')
    else:
        values.append(row['cell_type'])

adata.obs['major_celltypes'] = values
adata.obs

#%%
adata.obs.major_celltypes.isna().sum()
adata.obs.cell_type.isna().sum()

#%%
clusters.obs["cluster_from_seacells"][clusters.obs["cluster_from_seacells"] == "nan"] = np.nan
#%%
clusters.obs.isna().sum()
#%%
adata
#%%
## Save the data

adata.write_h5ad(clustersDir + 'atlas_embeddings_cell_labelled_cancer_clusters_from_mc.h5ad')


## Extra code maybe useful for CCC analysis
#%%
# adata.obs.dropna(subset=['cell_type', 'major_celltypes'], inplace=True)
# adata.obs.major_celltypes.isna().sum()
# adata.obs.cell_type.isna().sum()

# adata.obs['cell_type_treatment'] = adata.obs['treatment'].astype('str') + '_' + adata.obs['major_celltypes'].astype('str')
# adata.obs

# adata_primary = adata[(adata.obs['tissue'] == 'Primary')]
# adata_primary.obs['cell_type'] = adata_primary.obs['cell_type'].astype('category')
# adata_ascites = adata[(adata.obs['tissue'] == 'Ascites')]
# adata_ascites.obs['cell_type'] = adata_ascites.obs['cell_type'].astype('category')
# adata_metastasis = adata[(adata.obs['tissue'] == 'Metastasis')]
# adata_metastasis.obs['cell_type'] = adata_metastasis.obs['cell_type'].astype('category')

# adata_primary = adata[(adata.obs['tissue'] == 'Primary')]
# adata_primary.obs['major_celltypes'] = adata_primary.obs['major_celltypes'].astype('category')
# adata_ascites = adata[(adata.obs['tissue'] == 'Ascites')]
# adata_ascites.obs['major_celltypes'] = adata_ascites.obs['major_celltypes'].astype('category')
# adata_metastasis = adata[(adata.obs['tissue'] == 'Metastasis')]
# adata_metastasis.obs['major_celltypes'] = adata_metastasis.obs['major_celltypes'].astype('category')