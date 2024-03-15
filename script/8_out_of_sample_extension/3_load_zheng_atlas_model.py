

#%%
## imports
import scanpy as sc
import scgen
import numpy as np
import pandas as pd
import anndata

#%%
## load data
adata = sc.read_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Metacells/cancer_seacells_hdg_patients.h5ad')
genes = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'
adata.obs['tissue-treatment'] = adata.obs['tissue'].astype('str') + '_' + adata.obs['treatment'].astype('str')

sc.pp.log1p(adata)
hdg = pd.read_csv(genes, index_col=0)
adata.var['highly_variable'] = hdg.highly_variable
adata.var.highly_variable = adata.var.highly_variable.fillna(False)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)

#%%
model = scgen.SCGEN.load("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_HDG.pt", adata=adata)
# ValueError: Number of vars in adata_target not the same as source. Expected: 5191 Received: 5192

#%%
## I need to add the missing gene to the adata_target
## Identify the missing gene from genes and add it to adata_target
genes = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
missing_gene = genes[~genes.index.isin(adata.var_names)].index
missing_gene

#%%
missing_gene = 'ZBTB20-AS2'

## Create a new anndata with scanpy function with .var and .var_names equals to the ones of the first adata plus the missing gene and set its expression value to 0 in the .X matrix

# Create a new adata object by copying the original
new_adata = adata.copy()

# Create an array of zeros for the new gene expression
new_gene_expression = np.zeros((new_adata.shape[0], 1))  # assuming the axis 0 represents cells

# Add the new gene expression to the existing data
new_adata = new_adata.concatenate(anndata.AnnData(X=new_gene_expression), join='outer')

# Update the variable names to include the new gene
new_adata.var_names = list(adata.var_names) + ['ZBTB20-AS2']

# Add 'highly_variable' column to .var attribute and set True for the new gene
new_adata.var['highly_variable'] = True
new_adata.var.loc['ZBTB20-AS2', 'highly_variable'] = True
new_adata.var.drop(columns=['highly_variable-0'], inplace=True)

#%%
model = scgen.SCGEN.load("/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_HDG.pt", adata=new_adata)

### TypeError: '<' not supported between instances of 'float' and 'str'

#%%
## I need to convert the adata_target to a float
new_adata.X = new_adata.X.astype(float)


# %%
