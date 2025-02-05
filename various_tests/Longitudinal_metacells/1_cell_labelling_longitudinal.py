# longitudinal total cell labelling strategy

#%%
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)
from cell_labeller import assign_scores, actual_labeller, create_cancer_adata, immune_labeller, fibroblast_labeller

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/Longitudinal/'

## Cell labelling
#%%
adata = sc.read(initDir + "RnaAdataLongitudinalwithSubtypes_meta.h5ad")

### Creating annotation columns with Vasquez-Garcia et al. signatures genes
#%%
adata = assign_scores(adata)

### Cell labelling strategy
adata = actual_labeller(adata)

### Removing columns that are not needed

adata.obs.drop(columns = ['cell_labels_ratio', 'assignment', 'leiden-1.8'], inplace = True)

#%%
### Creating cell labels for cell subtypes

adata_immune = immune_labeller(adata)

immune_cells = pd.DataFrame(adata_immune.obs['cell_types'], index=adata_immune.obs.index)
immune_cells

adata.obs['cell_types'] = adata.obs['max'].astype("str")

adata.obs['cell_types'][adata.obs.index.isin(immune_cells.index)] = immune_cells['cell_types']

adata_fibroblast = fibroblast_labeller(adata)

fibro_cells = pd.DataFrame(adata_fibroblast.obs['cell_types'], index=adata_fibroblast.obs.index)
fibro_cells

adata.obs['cell_types'][adata.obs.index.isin(fibro_cells.index)] = fibro_cells['cell_types']

adata.obs.cell_types.astype('category')

adata.write(initDir + "RnaAdataLongitudinalwithSubtypes_meta_cell_labelled.h5ad")

#%%
## creating adata by cell type

adata = sc.read(initDir + "RnaAdataLongitudinalwithSubtypes_meta_cell_labelled.h5ad")
#### cancer
adata_cancer = sc.AnnData(X=adata[adata.obs['max'] == 'CancerMSK'].X, 
                            var=adata[adata.obs['max'] == "CancerMSK"].var, 
                            obs = adata[adata.obs['max'] == 'CancerMSK'].obs)

adata_cancer.write_h5ad(initDir + 'longitudinal_cancer_filt_norm_nolog.h5ad')

#### immune
adata_immune = sc.AnnData(X=adata[adata.obs['max'] == 'HematopoieticMSK'].X, 
                            var=adata[adata.obs['max'] == "HematopoieticMSK"].var, 
                            obs = adata[adata.obs['max'] == 'HematopoieticMSK'].obs)

adata_immune.write_h5ad(initDir + 'longitudinal_immune_filt_norm_nolog.h5ad')

#### fibroblasts
adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write_h5ad(initDir + 'longitudinal_fibroblast_filt_norm_nolog.h5ad')

#### endothelial
adata_endo = sc.AnnData(X=adata[adata.obs['max'] == 'EndothelialMSK'].X, 
                            var=adata[adata.obs['max'] == "EndothelialMSK"].var, 
                            obs = adata[adata.obs['max'] == 'EndothelialMSK'].obs)

adata_endo.write_h5ad(initDir + 'longitudinal_endothelial_filt_norm_nolog.h5ad')

# %%
