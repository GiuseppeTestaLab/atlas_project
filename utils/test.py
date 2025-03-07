#%%
import scanpy as sc
import os
# %%
#%%
ori_path = "/group/testa/Project/OvarianAtlasTestStep0/raw_data/integration/cells/cancer/cells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad"
adata = sc.read_h5ad(ori_path)
# %%
ori_path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells_backup/metacells/endothelial/seacells_assignment_hdg_patients.h5ad"
adata_ma = sc.read_h5ad(ori_path)
# %%
step0_path = "/group/testa/Project/OvarianAtlasTestStep0/raw_data/integration/cells/cancer/cells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad"
step0_adata = sc.read_h5ad(step0_path)

#%%
new_path = "/group/testa/Project/OvarianAtlasTestStep0/raw_data/metacells/endothelial"
rels = {"seacells_assignment_hdg_patients_marta_tables.h5ad", "seacells_assignment_hdg_patients_seed_2.h5ad"}
seed_adatas = [sc.read_h5ad(os.path.join(new_path, rel)) for rel in rels]

#%%
# adata_ma.obs["SEACell2"] = step0_path.obs["SEACell"]
# adata_ma.obs['SEACell_patient_tissue'] = adata_ma.obs['SEACell'].astype('str') + '_' + adata_ma.obs['paper_ID'].astype('str') + '_' + adata_ma.obs['tissue'].astype('str')
# adata_ma.obs['SEACell2_patient_tissue'] = adata_ma.obs['SEACell2'].astype('str') + '_' + adata_ma.obs['paper_ID'].astype('str') + '_' + adata_ma.obs['tissue'].astype('str')

def copy_seacell(adata, adata2, name):
    adata.obs[name] = adata2.obs["SEACell"]
    adata.obs[name + '_patient_tissue'] = adata.obs[name].astype('str') + '_' + adata.obs['paper_ID'].astype('str') + '_' + adata.obs['tissue'].astype('str')
#%%
adata_ma.obs['SEACell_patient_tissue'] = adata_ma.obs['SEACell'].astype('str') + '_' + adata_ma.obs['paper_ID'].astype('str') + '_' + adata_ma.obs['tissue'].astype('str')
# copy_seacell(adata_ma, step0_adata, "SEACell2")
copy_seacell(adata_ma, seed_adatas[0], "SEACell_seed1")
copy_seacell(adata_ma, seed_adatas[1], "SEACell_seed2")
#copy_seacell(adata_ma, seed_adatas[2], "SEACell_seed3")

# %%
import numpy as np
from sklearn import metrics

def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
    # return purity
    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix) 

#%%
import pandas as pd
def celltype_frac(x, col_name):
    val_counts = x[col_name].value_counts()
    return val_counts.values[0] / val_counts.values.sum()


def compute_celltype_purity(ad, col_name):
    """
    Compute the purity (prevalence of most abundant value) of the specified col_name from ad.obs within each metacell.
    @param: ad - AnnData object with SEACell assignment and col_name in ad.obs dataframe
    @param: col_name - (str) column name within ad.obs representing celltype groupings for each cell.
    """

    celltype_fraction = ad.obs.groupby('SEACell_patient_tissue').apply(lambda x: celltype_frac(x, col_name))
    celltype = ad.obs.groupby('SEACell_patient_tissue').apply(lambda x: x[col_name].value_counts().index[0])

    return pd.concat([celltype, celltype_fraction], axis=1).rename(columns={0: col_name, 1: f'{col_name}_purity'})

#%%

purity = compute_celltype_purity(adata_ma, "SEACell2_patient_tissue")


# %%
adata_ma_na = adata_ma[~adata_ma.obs.SEACell.isna(), :]
purity_na = compute_celltype_purity(adata_ma_na, "SEACell2_patient_tissue")

# %%
purity_seed1 = compute_celltype_purity(adata_ma, "SEACell_seed1_patient_tissue")
purity_seed2 = compute_celltype_purity(adata_ma, "SEACell_seed2_patient_tissue")
#purity_seed3 = compute_celltype_purity(adata_ma, "SEACell_seed3_patient_tissue")
# %%
def compute_celltype_purity(ad, col_name, starting):
    """
    Compute the purity (prevalence of most abundant value) of the specified col_name from ad.obs within each metacell.
    @param: ad - AnnData object with SEACell assignment and col_name in ad.obs dataframe
    @param: col_name - (str) column name within ad.obs representing celltype groupings for each cell.
    """

    celltype_fraction = ad.obs.groupby(starting).apply(lambda x: celltype_frac(x, col_name))
    celltype = ad.obs.groupby(starting).apply(lambda x: x[col_name].value_counts().index[0])

    return pd.concat([celltype, celltype_fraction], axis=1).rename(columns={0: col_name, 1: f'{col_name}_purity'})

# %%
purity_seed1 = compute_celltype_purity(adata_ma, "SEACell_seed1_patient_tissue",  "SEACell_seed1_patient_tissue")
purity_seed2 = compute_celltype_purity(adata_ma, "SEACell_seed2_patient_tissue",  "SEACell_seed1_patient_tissue")
#purity_seed3 = compute_celltype_purity(adata_ma, "SEACell_seed3_patient_tissue",  "SEACell_seed1_patient_tissue")

# %%
# %%
ori_path = "/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells_backup/metacells/cancer/seacells_assignment_hdg_patients.h5ad"
adata_ma = sc.read_h5ad(ori_path)

# %%
copy_seacell(seed_adatas[0], seed_adatas[1], "SEACell2")

# %%
seed_adatas[0].obs['SEACell_patient_tissue'] = seed_adatas[0].obs['SEACell'].astype('str') + '_' + seed_adatas[0].obs['paper_ID'].astype('str') + '_' + seed_adatas[0].obs['tissue'].astype('str')

purity = compute_celltype_purity(seed_adatas[0], "SEACell2_patient_tissue")
# %%
import pandas as pd
import configparser
# Read configuration file
config = configparser.ConfigParser()
config.read("config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

genes_marta = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_cancer.csv'
genes_miei = scriptsPath + '4_hdg/Tables_step0/atlas_hdg_dispersion_patients_cancer.csv'
#%%
gm = pd.read_csv(genes_marta, index_col=0)
gv = pd.read_csv(genes_miei, index_col=0)
#gv.reindex(gm.index)
# %%
