# atlas total cell labelling strategy

#%%
import numpy as np
import pandas as pd
import scanpy as sc
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from cell_labeller import assign_scores, actual_labeller, create_cancer_adata

## Cell labelling

adata = sc.read("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/original_anndata/atlas_embeddings.h5ad")
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

### Creating annotation columns with Vasquez-Garcia et al. signatures genes

adata = assign_scores(adata)

### Cell labelling strategy

adata = actual_labeller(adata)

adata.write('/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/atlas_embeddings_cell_labelled.h5ad')

