# Circos plot to visualize cell-to-cell interactions (cancer vs endothelial)

#%%
## Imports
import scanpy as sc
# Only needed for visualization:
from scanpy.pl import umap
import liana as li
import pandas as pd
import plotnine as p9
import numpy as np

## Initialize directories
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/LR_interactions/cancer_endothelial/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/LR_interactions/cancer_endothelial/'

## Primary
### Load the data
#%%
adata_primary = sc.read(initDir + 'liana_aggregaterank_primary.h5ad')
adata_primary
adata_primary.obs
adata_primary.uns['liana_res']

### Check for specific ligands and receptors
np.any(adata_primary.uns['liana_res'].ligand_complex == 'VEGFA')

np.any(adata_primary.uns['liana_res'].receptor_complex == 'VEGFR1')

np.any(adata_primary.uns['liana_res'].receptor_complex == 'VEGFR2')

np.any(adata_primary.uns['liana_res'].receptor_complex == 'VEGFR3')

### Subset for pairs of interest
adata_primary.uns['liana_res'] = adata_primary.uns['liana_res'][adata_primary.uns['liana_res'].ligand_complex == 'VEGFA']

### Plot bubble
my_plot = li.pl.dotplot(adata = adata_primary,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              top_n=20,
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

plot
plot.save(figDir + 'dotplot_rankaggr_primary_specific.png', limitsize=False)

## Metastasis
### Load the data
#%%
adata_metastasis = sc.read(initDir + 'liana_aggregaterank_metastasis.h5ad')
adata_metastasis
adata_metastasis.obs
adata_metastasis.uns['liana_res']

### Check for specific ligands and receptors
np.any(adata_metastasis.uns['liana_res'].ligand_complex == 'VEGFA')

np.any(adata_metastasis.uns['liana_res'].receptor_complex == 'VEGFR1')

np.any(adata_metastasis.uns['liana_res'].receptor_complex == 'VEGFR2')

np.any(adata_metastasis.uns['liana_res'].receptor_complex == 'VEGFR3')

### Subset for pairs of interest
adata_metastasis.uns['liana_res'] = adata_metastasis.uns['liana_res'][adata_metastasis.uns['liana_res'].ligand_complex == 'VEGFA']

### Plot bubble
my_plot = li.pl.dotplot(adata = adata_metastasis,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              top_n=20,
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

plot
plot.save(figDir + 'dotplot_rankaggr_metastasis_specific.png', limitsize=False)

## Ascites





