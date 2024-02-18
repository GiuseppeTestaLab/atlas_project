# Adata atlas fibroblasts integration

#%%
import scanpy as sc
import scgen
import sys
sys.path.insert(1, '/home/marta.sallese/ov_cancer_atlas/atlas_project/utils')
from integration import preprocess_scgen

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/metacells/fibroblasts/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/integration/cells/fibroblasts/'
#%%
ad = sc.read(initDir + 'seacells_assignment_hdg_patients.h5ad')
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
# inflammatory CAFs (iCAF)
sc.tl.score_genes(ad, ['IL6', 'CXCL12', 'VTN', 'CADM3', 'PLIN2', 'SERPINB2', 'KRT19', 'DES', 'CALB2', 'WT1', 'KRT7'], 
score_name = "iCAF", use_raw=False)

# matrix CAFs (mCAF)
sc.tl.score_genes(ad, ['ACTA2', 'COL11A1', 'MFAP5', 'SFRP2', 'ISLR', 'COL10A1', 'VIM', 'COL3A', 
                       'MMP11', 'PTHLH', 'FGF1', 'WNT7B', 'WNT2', 'TGFB3', 'THRC1', 'POSTN', 'VCAN', 'ZEB1'], 
score_name = "mCAF", use_raw=False)

# vascular cancer-associated fibroblasts CAFs (vCAF)
sc.tl.score_genes(ad, ['COX4I2', 'HIGD1B', 'PTP4A3', 'MCAM', 'PPP1R14A'], 
score_name = "vCAF", use_raw=False)

# STAR gene expressing CAFs (starCAF)
sc.tl.score_genes(ad, ['STAR', 'IGFBP5', 'TSPAN8', 'C7', 'ALDH1A1', 'LGR5'], 
score_name = "starCAF", use_raw=False)

#%%
ad.obs['cell_types'] = ad.obs[['iCAF', 'mCAF','vCAF', 'starCAF']].idxmax(axis=1)
ad.obs

#%%
batch = 'paper_ID'
genes = 2500
ad = preprocess_scgen(ad, batch, genes)

#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="cell_types")

#%%
model = scgen.SCGEN(ad)
model.save(outDir + "saved_models/model_scgen_batch_removal_celltypes.pt", overwrite=True)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'scgen_batch_corr_celltypes.h5ad')

