# Running Liana on cancer and immune cells divided by treatment type and tissue type

#%%
import scanpy as sc
# Only needed for visualization:
from scanpy.pl import umap
import liana as li
import pandas as pd
import plotnine as p9

## Inizialize directories
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + 'atlas_annotated/'
outDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/LR_interactions/cancer_immune/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/LR_interactions/cancer_immune/'


#%%
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled_with_ontologies.h5ad')
adata
adata.obs
adata = adata[adata.obs['ontologies_from_seacells'] != 'nan']
adata = adata[(adata.obs['max'] == 'CancerMSK') | (adata.obs['max'] == 'HematopoieticMSK')]


#%%
adata_cht = adata[(adata.obs['treatment'] == 'CHT')]
adata_naive = adata[(adata.obs['treatment'] == 'Naive')]
adata_nact = adata[(adata.obs['treatment'] == 'NACT')]

adata_primary = adata[(adata.obs['tissue'] == 'Primary')]
adata_metastasis = adata[(adata.obs['tissue'] == 'Metastasis')]
adata_ascites = adata[(adata.obs['tissue'] == 'Ascites')]

#%%
## Methods
li.mt.show_methods()

#%%
## import liana's rank_aggregate
from liana.mt import rank_aggregate
## resources
li.resource.show_resources()

#%%
## import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

#%%
# Run rank_aggregate
li.mt.rank_aggregate(adata_cht, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_naive, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_nact, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)

li.mt.rank_aggregate(adata_primary, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_metastasis, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)
li.mt.rank_aggregate(adata_ascites, groupby='ontologies_from_seacells', expr_prop=0.1, verbose=True, use_raw=False)

#%%
adata_cht.uns['liana_res'].head()
adata_naive.uns['liana_res'].head()
adata_nact.uns['liana_res'].head()

rank_aggregate.describe()

## Save the results
#%%
adata_cht.write_h5ad(outDir + 'liana_aggregaterank_cht.h5ad')
adata_naive.write_h5ad(outDir + 'liana_aggregaterank_naive.h5ad')
adata_nact.write_h5ad(outDir + 'liana_aggregaterank_nact.h5ad')

adata_primary.write_h5ad(outDir + 'liana_aggregaterank_primary.h5ad')
adata_metastasis.write_h5ad(outDir + 'liana_aggregaterank_metastasis.h5ad')
adata_ascites.write_h5ad(outDir + 'liana_aggregaterank_ascites.h5ad')

## Plotting the results

#%%
### Primary
adata_primary = sc.read(outDir + 'liana_aggregaterank_primary.h5ad')
set(adata_primary.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_primary,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'Extracellular_signaling',
                            'ILC',
                            'INF_mediated_signaling',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-cell_movement',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_primary'],
              target_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'Extracellular_signaling',
                            'ILC',
                            'INF_mediated_signaling',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-cell_movement',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_primary'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(80, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_primary.png', limitsize=False)
plot

#%%
### Ascites
adata_ascites = sc.read(outDir + 'liana_aggregaterank_ascites.h5ad')
set(adata_ascites.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_ascites,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'Response_to_extracellular_signals',
                            'Response_to_stress',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites'],
              target_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'Response_to_extracellular_signals',
                            'Response_to_stress',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(70, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_ascites.png', limitsize=False)
plot

#%%
### Metastasis
adata_metastasis = sc.read(outDir + 'liana_aggregaterank_metastasis.h5ad')
set(adata_metastasis.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_metastasis,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Cellular_metabolism',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ECM_shaping_cells',
                            'Extracellular_signaling-immune_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Plasma_cells',
                            'RNA_metabolism',
                            'Response_to_extracellular_signals',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_metastasis'],
              target_labels=['B_cells',
                            'Cellular_metabolism',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ECM_shaping_cells',
                            'Extracellular_signaling-immune_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Plasma_cells',
                            'RNA_metabolism',
                            'Response_to_extracellular_signals',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_metastasis'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(60, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_metastasis.png', limitsize=False)
plot

#%%
### Naive
adata_naive = sc.read(outDir + 'liana_aggregaterank_naive.h5ad')
set(adata_naive.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_naive,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ECM_shaping_cells',
                            'Extracellular_signaling',
                            'Extracellular_signaling-immune_cells',
                            'ILC',
                            'INF_mediated_signaling',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'RNA_metabolism',
                            'Response_to_extracellular_signals',
                            'Response_to_stress',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_metastasis',
                            'Unknown_cancer_primary'],
              target_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ECM_shaping_cells',
                            'Extracellular_signaling',
                            'Extracellular_signaling-immune_cells',
                            'ILC',
                            'INF_mediated_signaling',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'RNA_metabolism',
                            'Response_to_extracellular_signals',
                            'Response_to_stress',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_metastasis',
                            'Unknown_cancer_primary'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(80, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_naive.png', limitsize=False)
plot

#%%
### CHT
adata_cht = sc.read(outDir + 'liana_aggregaterank_cht.h5ad')
set(adata_cht.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_cht,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_movement',
                            'Plasma_cells',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_primary'],
              target_labels=['B_cells',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_movement',
                            'Plasma_cells',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_primary'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(60, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_cht.png', limitsize=False)
plot

#%%
### NACT
adata_nact = sc.read(outDir + 'liana_aggregaterank_nact.h5ad')
set(adata_nact.obs.ontologies_from_seacells)

my_plot = li.pl.dotplot(adata = adata_nact,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'Response_to_extracellular_signals',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_metastasis',
                            'Unknown_cancer_primary'],
              target_labels=['B_cells',
                            'Cellular_metabolism',
                            'Cellular_metabolism-extracellular_signaling',
                            'Ciliated_cancer_cells',
                            'Cycling_cells',
                            'Dendritic_cells',
                            'ILC',
                            'Immunoreactive_cells',
                            'M1_macrophages',
                            'Mast_cells',
                            'Myeloid_cells',
                            'NK_cells',
                            'Organelles_organization-cell_cycle',
                            'Organelles_organization-metabolism',
                            'Plasma_cells',
                            'Response_to_extracellular_signals',
                            'T_CD4_cells',
                            'T_CD8_ISG',
                            'T_CD8_cells',
                            'Unknown_cancer_ascites',
                            'Unknown_cancer_metastasis',
                            'Unknown_cancer_primary'],
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(60, 8)
             )

my_plot

plot = (my_plot +
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=8),
     # rotate labels on x axis
     axis_text_x=p9.element_text(rotation=60, hjust=1)
 )
)

plot.save(figDir + 'dotplot_rankaggr_nact.png', limitsize=False)
plot
