# atlas total cell labelling strategy

#%%
import numpy as np
import pandas as pd
import scanpy as sc

# #%%
# adata = sc.read('/group/testa/Project/OvarianAtlas/Geistlinger2020/Adata/geistlinger2020_rawcounts.h5ad') 

# data = ['/group/testa/Project/OvarianAtlas/Loret2022/Adata/loret2022_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Olbrecht2021/Adata/olbrecht2021_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Xu2022_CCR/Adata/xu2022_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Zhang2022/Adata/zhang2022_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Regner2021/Adata/regner2021_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Ren2022/Adata/ren2022_rawcounts.h5ad',
#         '/group/testa/Project/OvarianAtlas/Qian2020/Adata/qian2020_rawcounts.h5ad',
#        '/group/testa/Project/OvarianAtlas/Vasquez2022/Adata/vasquez2022_rawcounts.h5ad']
# #%%
# for i in data:
#     print('Reading adata' + i)
#     adata = adata.concatenate(sc.read(i), index_unique=None)
# #%%
# adata
# #%%
# adata.obs = adata.obs.drop(columns=['batch'])
# #%%
# adata.obs
# #%%
# set(adata.obs.tissue)
# #%%
# data = []
# for i in adata.obs.paper_ID:
#     data.append(i.split('_')[0])

# adata.obs.dataset = data
# #%%
# adata.obs.tissue[adata.obs.tissue == 'PELVIS'] = 'Metastasis'
# adata.obs.tissue = pd.Categorical(adata.obs.tissue, categories = set(adata.obs.tissue))
# #%%
# adata.obs.anatomical_location[adata.obs.anatomical_location == 'PELVIS'] = 'Other'
# adata.obs.anatomical_location = pd.Categorical(adata.obs.anatomical_location, categories = set(adata.obs.anatomical_location))
# #%%
# adata.write_h5ad('/group/testa/Project/OvarianAtlas/Adata_atlas_def/atlas_rawcounts_def.h5ad')
# #%%
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# #%%
# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# #%%
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)
# #%%
# sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save = 'pct_counts_mt.png')
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = 'n_genes_by_counts.png')
# #%%
# adata = adata[adata.obs.n_genes_by_counts < 7000, :]
# adata = adata[adata.obs.pct_counts_mt < 30, :]
# #%%
# sc.pp.normalize_total(adata, target_sum=1e4)
# #%%
# adata.write("/group/testa/Project/OvarianAtlas/Adata_atlas_def/atlas_filt_norm_nolog_def.h5ad")
# #%%
# sc.pp.log1p(adata)
# #%%
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# #%%
# sc.pl.highly_variable_genes(adata)
# #%%
# adata.raw = adata
# #%%
# adata = adata[:, adata.var.highly_variable]
# #%%
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# #%%
# sc.pp.scale(adata, max_value=10)
# #%%
# sc.tl.pca(adata, svd_solver='arpack')
# #%%
# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
# #%%
# sc.tl.umap(adata)
# #%%
# adata.write("/group/testa/Project/OvarianAtlas/Adata_atlas_def/atlas_preprocessed_def.h5ad")

## Cell labelling

adata = sc.read("/group/testa/Project/OvarianAtlas/Adata_atlas_def/atlas_preprocessed_def.h5ad")

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

adata.write('/group/testa/Project/OvarianAtlas/Adata_atlas_def/atlas_preprocessed_cell_labelled_def.h5ad')

