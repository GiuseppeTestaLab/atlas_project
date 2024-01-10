
import scanpy as sc
import numpy as np

def assign_scores(adata):
    
    sc.tl.score_genes(adata, 
    ['WFDC2','CD24','CLDN3','KRT7','KRT8','KRT17','KRT18','KRT19','EPCAM','WT1','CLDN4','MSLN','FOLR1','MUC1'], 
    score_name = "CancerMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['CLEC14A','PECAM1','VWF','CAV1','EMCN','CDH5','MCAM','IL3RA','IGFBP7','COL4A1','COL4A2','COL15A1','SPARCL1','A2M','HSPG2','PLVAP','AQP1','ENG','RAMP2','GNG11','EGFL7','CLDN5','INSR'], 
    score_name = "EndothelialMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['COL1A1','COL3A1','WT1','ACTA2','CAV1','COL1A2','DCN','SPARC','COL6A1','CCDC80','LUM','COL6A2','COL6A3','CALD1','RARRES2','MGP','CTHRC1','AEBP1','POSTN','COL5A2','FBLN1','TAGLN','C1S','C1R','NNMT','MMP2','IGFBP5','TIMP1','FN1','IGFBP7','C3','COL5A1','LGALS1'], 
    score_name = "FibroblastsMSK", use_raw=True)

    sc.tl.score_genes(adata, 
    ['PTPRC'], 
    score_name = "HematopoieticMSK", use_raw=True)

    return adata

# Cell labelling strategy

def actual_labeller(adata):
    def positive_ratio(x):
        x[x<0]=0
        if (sum(x)-max(x))==0:
            return max(x)
        else:
            return max(x)/(sum(x)-max(x))

    adata.obs['cell_labels_ratio']= adata.obs[['EndothelialMSK', 'FibroblastsMSK', 'CancerMSK', 'HematopoieticMSK']].apply(lambda x:positive_ratio(x), axis=1)

    np.log(adata.obs.loc[adata.obs.cell_labels_ratio > 0.1, 'cell_labels_ratio']).hist(bins=30)

    adata.obs['max'] = adata.obs[['EndothelialMSK', 'FibroblastsMSK', 'CancerMSK','HematopoieticMSK']].idxmax(axis=1)

    percentage = 0.2
    counter = 0
    for dataset in adata.obs.dataset.unique():
        print("Processing the dataset: {}".format(dataset))
        adata_slice = sc.AnnData(X=adata[adata.obs.dataset == dataset].X, 
                                var=adata[adata.obs.dataset == dataset].var, 
                                obs = adata[adata.obs.dataset == dataset].obs) 
        sc.tl.pca(adata_slice, svd_solver='arpack')
        sc.pp.neighbors(adata_slice, n_neighbors=15, n_pcs=50)
        sc.tl.leiden(adata_slice,resolution = 1.8, key_added="leiden-{}".format('1.8'))
        sc.tl.umap(adata_slice)
        for leiden in adata_slice.obs["leiden-1.8"].unique():
            n_cells_xpopul = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].value_counts()[0]
            n_cells_total = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].shape[0]
            assignment = adata_slice.obs.loc[adata_slice.obs['leiden-1.8'] == leiden, 'max'].value_counts().index[0]
            if (n_cells_xpopul / n_cells_total) > percentage:
                adata_slice.obs["assignment"] = assignment
        #putting the assignment in the initial adata
        if counter == 0:
            adata.obs['assignment'] = adata_slice.obs['assignment']
            adata.obs['leiden-1.8'] = adata_slice.obs['leiden-1.8']
        else :
            adata.obs.loc[adata_slice.obs.index, 'assignment'] = adata_slice.obs['assignment']
            adata.obs['leiden-1.8'] = adata.obs['leiden-1.8'].astype('str')
            adata.obs.loc[adata_slice.obs.index, 'leiden-1.8'] = adata_slice.obs['leiden-1.8'].astype('str')
        counter = counter + 1

    adata = adata.raw.to_adata()

    return adata

def create_cancer_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'CancerMSK'].X, 
                            var=adata[adata.obs['max'] == "CancerMSK"].var, 
                            obs = adata[adata.obs['max'] == 'CancerMSK'].obs)
    
    return adata

def create_immune_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'HematopoieticMSK'].X, 
                            var=adata[adata.obs['max'] == "HematopoieticMSK"].var, 
                            obs = adata[adata.obs['max'] == 'HematopoieticMSK'].obs)
    
    return adata

def create_fibroblast_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastMSK'].obs)
    
    return adata

def create_endothelial_adata(adata):
    adata = sc.AnnData(X=adata[adata.obs['max'] == 'EndothelialMSK'].X, 
                            var=adata[adata.obs['max'] == "EndothelialMSK"].var, 
                            obs = adata[adata.obs['max'] == 'EndothelialMSK'].obs)
    
    return adata