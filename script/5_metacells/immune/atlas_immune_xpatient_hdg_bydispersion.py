# Atlas metacells generation

## Import libraries
#%%
import scanpy as sc
import sys
import configparser
import os
# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

sys.path.insert(1, utilsPath)
from metacells_derivation import preprocess, assign_metacells, create_mc_matrix, preprocess_mc # type: ignore

#%%
initDir = rawPath + 'atlas_annotated/'
destDir = rawPath + 'metacells/immune/'
if not os.path.exists(destDir):
    os.makedirs(destDir)
scriptsPath = config.get("DEFAULT", "scriptsPath")

## Load data
#%%
adata= sc.read(initDir + "atlas_immune_filt_norm_nolog.h5ad")
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_immune.csv'

## Preprocessing
#%%
adata = preprocess(adata, genes)
adata.write_h5ad(initDir + 'atlas_immune_embeddings.h5ad')

## Metacells generation per patient
#%%
adata = assign_metacells(adata)
adata.write_h5ad(destDir + 'seacells_assignment_hdg_patients.h5ad')

# Creating metacell matrix
#%%
adata = sc.read(destDir + 'seacells_assignment_hdg_patients.h5ad')
ad = create_mc_matrix(adata)

ad.obs.drop(columns = ['ID', 'sample_name', 'patient_id', 'cell_type', 'cell_subtype', 'sample_ID',
       'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt',
       'pct_counts_mt', 'CancerMSK', 'EndothelialMSK', 'FibroblastsMSK',
       'HematopoieticMSK', 'cell_labels_ratio', 'max', 'assignment',
       'leiden-1.8'])

ad.write(destDir + 'seacells_hdg_patients.h5ad')

## Compute embeddings and plot metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "metacells/immune/"

adata = sc.read(destDir + 'seacells_hdg_patients.h5ad')

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

adata = preprocess_mc(adata, genes)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='immune_seacells_HDG_treatm.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='immune_seacells_HDG_tissue.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='immune_seacells_HDG_dataset.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='immune_seacells_HDG_patients.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='immune_seacells_HDG_cellcycle.png')
sc.pl.umap(adata, color=["anatomical_location"], frameon=False, save='immune_seacells_HDG_anatomy.png')
sc.pl.umap(adata, color=["cell_types"], frameon=False, save='immune_seacells_HDG_cell_types.png')

#%%
adata.write(destDir + 'seacells_hdg_patients_embeddings.h5ad')

