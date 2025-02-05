# Integration (preserving tissues and treatment) of cancer seacells generated in the HDG space with prior HDG space selection

#%%
import scanpy as sc
import scgen
import pandas as pd
import numpy as np
from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField

#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'metacells/cancer/'
outDir = rawPath + 'integration/metacells/cancer/'
genes = scriptsPath + 'cancer.csv'

#%%
ad = sc.read_h5ad(initDir + "seacells_hdg_patients.h5ad")
ad.obs['tissue-treatment'] = ad.obs['tissue'].astype('str') + '_' + ad.obs['treatment'].astype('str')

#%%
sc.pp.log1p(ad)
hdg = pd.read_csv(genes, index_col=0)
ad.var['highly_variable'] = hdg.highly_variable
ad.var.highly_variable = ad.var.highly_variable.fillna(False)
ad.raw = ad
ad = ad[:, ad.var.highly_variable]
sc.tl.pca(ad, use_highly_variable=True)
sc.pp.neighbors(ad, use_rep='X_pca')
sc.tl.umap(ad)

## Creating an adata manager for ad in order to have to possibility to extend categories when creating the model
## This part works if scvi env is used
#%%
ad_fields = [LayerField(registry_key="x", layer=None, is_count_data=True),
             CategoricalObsField(registry_key="batch", attr_key="paper_ID")
            ]

#%%
adata_manager = AnnDataManager(fields=ad_fields)
adata_manager.register_fields(ad)

adata_manager = adata_manager.transfer_fields(
    adata_target=adata_new, extend_categories=True
)

adata_manager.view_registry()

print(
    adata_manager.registry.keys()
)

data_registry = adata_manager.data_registry
data_registry

print(data_registry["batch"])
print(data_registry.batch.attr_key)

batch_state_registry = adata_manager.get_state_registry("batch")
print(batch_state_registry.keys())

print(f"Categorical mapping: {batch_state_registry.categorical_mapping}")
print(f"Original key: {batch_state_registry.original_key}")

adata_manager.summary_stats

#%%
print(adata_manager.fields)

## Loading adata target
#%%
adata_target = sc.read_h5ad('/group/testa/Project/OvarianAtlas/Zheng2023/Metacells/cancer_seacells_hdg_patients.h5ad')
adata_target.obs['tissue-treatment'] = adata_target.obs['tissue'].astype('str') + '_' + adata_target.obs['treatment'].astype('str')
#%%
sc.pp.log1p(adata_target)
hdg = pd.read_csv(genes, index_col=0)
adata_target.var['highly_variable'] = hdg.highly_variable
adata_target.var.highly_variable = adata_target.var.highly_variable.fillna(False)
adata_target.raw = adata_target
adata_target = adata_target[:, adata_target.var.highly_variable]
sc.tl.pca(adata_target, use_highly_variable=True)
sc.pp.neighbors(adata_target, use_rep='X_pca')
sc.tl.umap(adata_target)
#%%
genes = pd.read_csv(scriptsPath + 'cancer.csv', index_col=0)
missing_gene = genes[~genes.index.isin(adata_target.var_names)].index
missing_gene
missing_gene = 'ZBTB20-AS2'
#%%
new_adata = adata_target.copy()
adataX=pd.DataFrame(adata_target.X.todense().T,index=adata_target.var_names,columns=adata_target.obs_names)
# adata_target
# list(adata_target.var_names) + ['ZBTB20-AS2']
adataX=adataX.T
adataX['ZBTB20-AS2']=adataX['MT-CYB']*0
adata_new=sc.AnnData(adataX)
adata_new.obs=adata_target.obs
adata_new

## Creating an adata manager also for the adata target resembling the one created for ad
#%%
adata_manager = adata_manager.transfer_fields(
    adata_target=adata_new, extend_categories=True
)
adata_manager.view_registry()

## Creating the model 
#%%
scgen.SCGEN.setup_anndata(ad, batch_key="paper_ID", labels_key="tissue-treatment", )

#%%
model = scgen.SCGEN(ad)
model.save(rawPath + "integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_HDG_oose.pt", overwrite=True)

#%%
model = scgen.SCGEN.load(rawPath + "integration/metacells/saved_models/cancer_batch_removal_tissuetreatment_HDG_oose.pt", adata=adata_new)

#%%
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

#%%
corrected_adata = model.batch_removal()
corrected_adata.write_h5ad(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_HDG.h5ad')

## Processing of integrated metacells in the same HDG space used to generate metacells
#%%
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')
sc.settings.figdir = figPath + "integration/metacells/cancer/"

#%%
adata = sc.read(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_HDG.h5ad')

#%%
cell_cycle_genes = [x.strip() for x in open(CCGenes)]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#%%
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#%%
sc.pl.umap(adata, color=["treatment"], frameon=False, save='cancer_seacells_HDG_treatm_HDG.png')
sc.pl.umap(adata, color=["tissue"], frameon=False, save='cancer_seacells_HDG_tissue_HDG.png')
sc.pl.umap(adata, color=["dataset"], frameon=False, save='cancer_seacells_HDG_dataset_HDG.png')
sc.pl.umap(adata, color=["paper_ID"], frameon=False, save='cancer_seacells_HDG_patients_HDG.png')
sc.pl.umap(adata, color=["phase"], frameon=False, save='cancer_seacells_HDG_cellcycle_HDG.png')
sc.pl.umap(adata, color=["tissue-treatment"], frameon=False, save='cancer_seacells_HDG_tissue-treat_HDG.png')

#%%
adata.write(outDir + 'seacells_hdg_patients_batch_corr_scgen_tissuetreat_embeddings_HDG.h5ad')