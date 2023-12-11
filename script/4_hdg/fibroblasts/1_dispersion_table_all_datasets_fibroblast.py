#%%
import scanpy as sc
import pandas as pd
import numpy as np
import glob

#%%
initialize directories

initDir = "/group/testa/Project/OvarianAtlas/"
datasetName = ['regner2021', 'geistlinger2020', 'qian2020', 'ren2022', 'loret2022', 'olbrecht2021', 'xu2022_CCR', 'zhang2022', 'vasquez2022']
common_var_names = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/common_varnames_datasets.csv', index_col=0) #this file is derived from common_var_names.py

for dat in datasetName:
    path = (initDir + dat.capitalize() + '/Adata/{}_filt_norm_nolog.h5ad').format(dat)
    if dat == 'xu2022_CCR':
       path = initDir + 'Xu2022_CCR' + '/Adata/' + "xu2022_filt_norm_nolog.h5ad"
    print(path)
    adata = sc.read(path)
    adata = adata[:, common_var_names.index]
    dataset = path.split('/')[-1].split('_')[0]
    adata.write(("/group/testa/Project/OvarianAtlas/Adata_common_genes_filt_norm_nolog/{}_common_genes_filt_norm_nolog.h5ad").format(dataset))

commonDir = "/group/testa/Project/OvarianAtlas/Adata_common_genes_filt_norm_nolog/"

#%%
## Geistlinger
adata = sc.read(commonDir + "geistlinger2020_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Geistlinger2020/geistlinger2020_adata_fibroblasts.h5ad')

#%%
## Regner
adata = sc.read(commonDir + "regner2021_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Regner2021/regner2021_adata_fibroblasts.h5ad')

#%%
## Loret
adata = sc.read(commonDir + "loret2022_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Loret2022/loret2022_adata_fibroblasts.h5ad')

#%%
## Olbrecht
adata = sc.read(commonDir + "olbrecht2021_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Olbrecht2021/olbrecht2021_adata_fibroblasts.h5ad')

## Qian
adata = sc.read(commonDir + "qian2020_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Qian2020/qian2020_adata_fibroblasts.h5ad')

#%%
## Ren
adata = sc.read(commonDir + "ren2022_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Ren2022/ren2022_adata_fibroblasts.h5ad')

#%%
## Xu
adata = sc.read(commonDir + "xu2022_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Xu2022/xu2022_adata_fibroblasts.h5ad')

## Zhang
adata = sc.read(commonDir + "zhang2022_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Zhang2022/zhang2022_adata_fibroblasts.h5ad')

#%%
## Vasquez
adata = sc.read(commonDir + "vasquez2022_common_genes_filt_norm_nolog.h5ad")
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs['dataset'] = adata.obs.paper_ID.str.split("_").str[0]

## Creating annotation columns with Vasquez-Garcia et al. signatures genes

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

## Cell labelling strategy

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

adata_fibro = sc.AnnData(X=adata[adata.obs['max'] == 'FibroblastsMSK'].X, 
                            var=adata[adata.obs['max'] == "FibroblastsMSK"].var, 
                            obs = adata[adata.obs['max'] == 'FibroblastsMSK'].obs)

adata_fibro.write(commonDir + 'Vasquez2022/vasquez2022_adata_fibroblasts.h5ad')


#%%

dataDir = "/group/testa/Project/OvarianAtlas/Adata_common_genes_filt_norm_nolog/"
dataName = ['regner2021', 'geistlinger2020', 'qian2020', 'ren2022', 'loret2022', 'olbrecht2021', 'xu2022', 'zhang2022', 'vasquez2022']

common_var_names = pd.read_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/common_varnames_datasets.csv', index_col=0)

dispersion_table = pd.DataFrame(index = common_var_names.index)
hvg_table = pd.DataFrame(index = common_var_names.index)

for j in dataName:
    path = (dataDir + j.capitalize() + "/{}_adata_fibroblasts.h5ad").format(j)
    print(path)
    adata_fibro = sc.read(path)
    dispersion_gene_xpatient = {}
    highly_variable_genes_per_patient = {}
    for patient_code in adata_fibro.obs["paper_ID"].unique():
        patient_anndata = adata_fibro[adata_fibro.obs["paper_ID"] == patient_code].copy()
        sc.pp.highly_variable_genes(patient_anndata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        highly_variable_genes_per_patient[patient_code] = patient_anndata.var.highly_variable
        dispersion_gene_xpatient[patient_code] = patient_anndata.var.dispersions_norm
    for i in dispersion_gene_xpatient:
        print(i)
        print(dispersion_gene_xpatient[i])
        dispersion_table[i] = dispersion_gene_xpatient[i]
        hvg_table[i] = highly_variable_genes_per_patient[i]

dispersion_table.to_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/dispersion_table_fibroblasts.csv')
hvg_table.to_csv('/home/marta.sallese/ov_cancer_atlas/Atlas_scripts/HDG_new/Tables/hvg_table_fibroblasts.csv')

# %%
