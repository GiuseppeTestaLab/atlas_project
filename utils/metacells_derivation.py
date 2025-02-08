
import scanpy as sc
import numpy as np
import pandas as pd
import SEACells
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

CCGenes = config.get("DEFAULT", "CCGenes")


## Preprocessing

def preprocess(adata, genes):
    sc.pp.neighbors(adata)
    raw_ad = sc.AnnData(adata.X)
    raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
    adata.raw = raw_ad
    sc.pp.log1p(adata)
    print("reading HDG")
    hdg = pd.read_csv(genes, index_col=0)
    adata.var['highly_variable'] = hdg.highly_variable
    adata.var.highly_variable = adata.var.highly_variable.fillna(False)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.tl.umap(adata)

    return(adata)

def assign_metacells(adata):
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
                                    convergence_epsilon = 1e-5,verbose=False,use_gpu=False)
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
        
    sc.settings.verbosity = 3  

    import pickle

    with open('models.pkl','wb') as file:
        pickle.dump(models,file)

    with open('models.pkl','rb') as file:

        modelsLoaded=pickle.load(file)

    adata.obs['SEACell'] = mc_obs.SEACell
    return(adata)

def create_mc_matrix(adata):
    adata.obs['SEACell_patient_tissue'] = adata.obs['SEACell'].astype('str') + '_' + adata.obs['paper_ID'].astype('str') + '_' + adata.obs['tissue'].astype('str')
    adata = adata[~adata.obs.SEACell.isna(), :]
    adata_obs = adata.obs.copy()

    adata_obs = adata_obs.set_index('SEACell_patient_tissue')
    adata_obs = adata_obs.loc[~adata_obs.index.duplicated(keep='first'), :]
    
    from SEACells.core import summarize_by_SEACell
    ad = summarize_by_SEACell(adata, SEACells_label='SEACell_patient_tissue')
    ad.layers['raw'] = ad.X
    ad.obs = adata_obs
    # ad.obs = ad.obs.drop(columns=['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 'sample_ID', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'])
    ad.obs['# Single Cells'] = adata.obs.groupby('SEACell_patient_tissue').count().tissue
    return(ad)

def preprocess_mc(adata, genes):
    cell_cycle_genes = [x.strip() for x in open(CCGenes)]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    hdg = pd.read_csv(genes,  index_col=0)
    adata.var['highly_variable']=hdg.highly_variable
    adata.var.highly_variable = adata.var.highly_variable.fillna(False)
    adata.var

    sc.tl.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata)
    return(adata)
