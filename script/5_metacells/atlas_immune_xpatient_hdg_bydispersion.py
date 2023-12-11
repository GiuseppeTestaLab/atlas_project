# Atlas metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'
destDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/'

## Load data
#%%
adata= sc.read(initDir + "atlas_immune_filt_norm_nolog.h5ad")

data = []
for i in adata.obs.paper_ID:
    data.append(i.split('_')[0])

adata.obs.dataset = data

adata.obs.tissue[adata.obs.tissue == 'PELVIS'] = 'Metastasis'
adata.obs.tissue = pd.Categorical(adata.obs.tissue, categories = set(adata.obs.tissue))

adata.obs.anatomical_location[adata.obs.anatomical_location == 'PELVIS'] = 'Peritoneum'
adata.obs.anatomical_location = pd.Categorical(adata.obs.anatomical_location, categories = set(adata.obs.anatomical_location))


## Preprocessing
#%%
sc.pp.neighbors(adata)
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad
sc.pp.log1p(adata)
hdg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/hdg/immune/atlas_hdg_dispersion_patients_immune.csv',  index_col=0)
hdg[hdg.highly_variable]
adata.var['highly_variable']=hdg.highly_variable
adata.var.highly_variable = adata.var.highly_variable.fillna(False)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.tl.umap(adata)
adata.write_h5ad(initDir + 'atlas_immune_embeddings.h5ad')

## Metacells generation per patient
#%%
sc.settings.verbosity = 0
models = {}
mc_obs = pd.DataFrame()
i=0
for patient in adata.obs.paper_ID.unique():
    ad_patient=adata[adata.obs.paper_ID==patient].copy()
    n_SEACells = np.round(ad_patient.shape[0]/75, 0)
    build_kernel_on = 'X_pca'
    n_waypoint_eigs = 10 
    model = SEACells.core.SEACells(ad_patient, 
                                   build_kernel_on=build_kernel_on, 
                                   n_SEACells=n_SEACells, 
                                   n_waypoint_eigs=n_waypoint_eigs, 
                                   convergence_epsilon = 1e-5,verbose=False)
    model.construct_kernel_matrix()
    M = model.kernel_matrix 
    try:
        model.initialize_archetypes()
    except IndexError as e:
        print('patient {} the number {} has too few cells'.format(patient,i))
        ad = adata[adata.obs.paper_ID != patient]
        continue
    model.fit(min_iter=10, max_iter=200)
    models[patient]=model
    labels,weights = model.get_soft_assignments()
    SEACell_soft_ad = SEACells.core.summarize_by_soft_SEACell(ad_patient, model.A_, celltype_label='paper_ID',summarize_layer='raw', minimum_weight=0.05)
    if mc_obs.shape[0] == 0:
        # mc_obs=SEACell_soft_ad
        mc_obs=ad_patient.obs
    else:
        mc_obs=pd.concat([mc_obs,ad_patient.obs])
        # mc_obs=mc_obs.concatenate(SEACell_soft_ad)
    i+=1

    print('done! patient {} the number {}'.format(patient,i))
    #if i>2:
    #    break
sc.settings.verbosity = 3  

import pickle

with open('models.pkl','wb') as file:
    pickle.dump(models,file)


with open('models.pkl','rb') as file:

    modelsLoaded=pickle.load(file)

adata.obs['SEACell'] = mc_obs.SEACell

adata.write_h5ad(destDir + 'atlas_seacells_assignment_hdg_patients_immune.h5ad')

## Creating metacell matrix
#%%
# adata = sc.read(destDir + 'atlas_seacells_assignment_hdg_patients_immune.h5ad')

adata.obs
adata

adata.obs['SEACell_patient_tissue'] = adata.obs['SEACell'].astype('str') + '_' + adata.obs['paper_ID'].astype('str') + '_' + adata.obs['tissue'].astype('str')
adata = adata[~adata.obs.SEACell.isna(), :]
adata_obs = adata.obs.copy()

adata_obs = adata_obs.set_index('SEACell_patient_tissue')
adata_obs 

adata_obs = adata_obs.loc[~adata_obs.index.duplicated(keep='first'), :]

#%%
from SEACells.core import summarize_by_SEACell

ad = summarize_by_SEACell(adata, SEACells_label='SEACell_patient_tissue')
ad.layers['raw'] = ad.X

#%%
ad.obs = adata_obs
ad

ad.obs = ad.obs.drop(columns=['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 'sample_ID', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'])

#%%
ad.obs['# Single Cells'] = adata.obs.groupby('SEACell_patient_tissue').count().tissue
ad.obs

#%%
ad.write(destDir + 'atlas_immune_seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells

#%%
adata = sc.read(destDir + 'atlas_immune_seacells_hdg_patients.h5ad')

#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = "/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_def/metacells/"

#%%
cell_cycle_genes = [x.strip() for x in open('/home/marta.sallese/ov_cancer_atlas/regev_lab_cell_cycle_genes.txt')]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#%%
# Plasma cells markers
sc.tl.score_genes(adata, ['IGKC','IGHG1','CD79A','IGHG2','IGLC2','IGLC3','IGHG3','IGHG4','JCHAIN','MZB1','XBP1'], 
score_name = "Plasma_cells", use_raw=False)

# T cells markers
sc.tl.score_genes(adata, ['CD2','CD3D','TRAC','GZMA','NKG7','CD3E','CD3G','CD4','TCF7','CD8A','PRF1','GZMB','CCL5','CCL4','IL32','CD52'], 
score_name = "T_cells", use_raw=False)

# Mast cells markers
sc.tl.score_genes(adata, ['KIT','CPA3','CTSG','MS4A2','TPSAB1','TPSB2','HPGD','HPGDS','GATA2'], 
score_name = "Mast_cells", use_raw=False)

# B cells markers
sc.tl.score_genes(adata, ['MS4A1', 'CD79A', 'CD19', 'BANK1', 'IGKC', 'IGHM'], score_name = "B_cells", use_raw=False)

# Myeloid cells markers
sc.tl.score_genes(adata, 
['CD14','FCER1G','FCGR3A','LYZ','CTSS','CD33','CD68','CD163','ITGAX','ITGAM','CD4','MRC1',
'VSIG4','SPP1','APOE','C1QA','C1QB','C1QC','APOC1','FTL','S100A9','TYROBP','AIF1','CD74','PSAP','CTSB'], 
score_name = "Myeloid_cells", use_raw=False)

# Dendritic cells markers
sc.tl.score_genes(adata, 
['IL3RA','IRF7','IRF8','GZMB','CD4','CLEC4C','JCHAIN',
'PTGDS','PLAC8','PLD4','TCF4','BCL11A','GPR183','CCDC50','LILRA4','TSPAN13','CLIC3','MPEG1'], 
score_name = "Dendritic_cells", use_raw=False)

#%%
adata.obs
adata
adata.obs['cell_types'] = adata.obs[['Plasma_cells', 'T_cells','Mast_cells', 'B_cells', 'Myeloid_cells', 'Dendritic_cells' ]].idxmax(axis=1)

#%%
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#%%
hdg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/atlas_project/script/hdg/immune/atlas_hdg_dispersion_patients_immune.csv',  index_col=0)
adata.var['highly_variable']=hdg.highly_variable

adata.var.highly_variable = adata.var.highly_variable.fillna(False)
adata.var

#%%
sc.tl.pca(adata, use_highly_variable=True)

#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='immune_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='immune_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='immune_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='immune_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='immune_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='immune_seacells_HDG_anatomy.png')
sc.pl.umap(adata, color=["cell_types"], frameon=False, save='immune_seacells_HDG_cell_types.png')

#%%
adata.write(destDir + 'atlas_immune_seacells_hdg_patients_embeddings.h5ad')

