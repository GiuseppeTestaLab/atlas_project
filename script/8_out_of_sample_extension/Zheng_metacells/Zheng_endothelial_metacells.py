# Zheng endothelial metacells generation

## Import libraries
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

## Load data
#%%
adata= sc.read("/group/testa/Project/OvarianAtlas/Zheng2023/Adata/zheng_endothelial_filt_norm_nolog.h5ad")


## Preprocessing
#%%
sc.pp.neighbors(adata)
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad
sc.pp.log1p(adata)
hdg = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/atlas_hdg_dispersion_patients_endothelial.csv',  index_col=0)
hdg[hdg.highly_variable]
adata.var['highly_variable']=hdg.highly_variable
adata.var.highly_variable = adata.var.highly_variable.fillna(False)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.tl.umap(adata)
#adata.write_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Adata/zheng_endothelial_embeddings.h5ad')

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

adata.write_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Metacells/endothelial_seacells_assignment_hdg_patients.h5ad')