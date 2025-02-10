# Running Liana on fibroblasts and endothelial cells divided by treatment type and tissue type

#%%
import scanpy as sc
# Only needed for visualization:
from scanpy.pl import umap
import liana as li
import pandas as pd
import plotnine as p9
import os
## Inizialize directories
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'atlas_annotated/'
outDir = rawPath + 'downstream/LR_interactions/fibroblasts_endothelial/'
if not os.path.exists(outDir):
    os.makedirs(outDir)
sc.settings.figdir = figPath + 'LR_interactions/fibroblasts_endothelial/'
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
adata = sc.read(initDir + 'atlas_embeddings_cell_labelled_with_ontologies.h5ad')
adata
adata.obs
adata = adata[adata.obs['ontologies_from_seacells'] != 'nan']
adata = adata[(adata.obs['max'] == 'FibroblastsMSK') | (adata.obs['max'] == 'EndothelialMSK')]


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

#%%
# my_plot = li.pl.dotplot(adata = adata_cht,
#               colour='magnitude_rank',
#               size='specificity_rank',
#               inverse_size=True,
#               inverse_colour=True,
#               source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               top_n=10,
#               orderby='magnitude_rank',
#               orderby_ascending=True,
#               figure_size=(27, 10)
#              )

# my_plot

# plot = (my_plot +
#  p9.theme(
#      # adjust facet size
#      strip_text=p9.element_text(size=11)
#  )
# )

# plot
# plot.save('Figures/dotplot_rankaggr_cht.png', limitsize=False)

# #%%
# my_plot = li.pl.dotplot(adata = adata_naive,
#               colour='magnitude_rank',
#               size='specificity_rank',
#               inverse_size=True,
#               inverse_colour=True,
#               source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               top_n=10,
#               orderby='magnitude_rank',
#               orderby_ascending=True,
#               figure_size=(27, 10)
#              )

# my_plot

# plot = (my_plot +
#  p9.theme(
#      # adjust facet size
#      strip_text=p9.element_text(size=11)
#  )
# )

# plot
# plot.save('Figures/dotplot_rankaggr_naive.png', limitsize=False)

# #%%
# my_plot = li.pl.dotplot(adata = adata_nact,
#               colour='magnitude_rank',
#               size='specificity_rank',
#               inverse_size=True,
#               inverse_colour=True,
#               source_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               target_labels=['Cancer_primary', 'Cancer_ascites', 'Cancer_metastasis', 'Plasma_cells','Mast_cells', 'B_cells', 
#                              'M1_macrophages', 'Myeloid_cells', 'Dendritic_cells', 'T_CD4_cells', 'T_CD8_cells', 'NK_cells'],
#               top_n=10,
#               orderby='magnitude_rank',
#               orderby_ascending=True,
#               figure_size=(27, 10)
#              )

# my_plot

# plot = (my_plot +
#  p9.theme(
#      # adjust facet size
#      strip_text=p9.element_text(size=11)
#  )
# )

# plot
# plot.save('Figures/dotplot_rankaggr_nact.png', limitsize=False)

