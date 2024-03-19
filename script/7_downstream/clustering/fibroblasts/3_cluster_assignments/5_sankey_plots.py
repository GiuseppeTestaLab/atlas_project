

## Imports
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import kaleido
import matplotlib
import plotly.express as px
import plotly.graph_objects as go

## Initialize directiories
tissueDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/fibroblasts/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/'

#%%
## Loading data
primary = sc.read_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

#%%
## Plotting a sankey for primary metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = primary.obs.groupby(['cell_types', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(primary.obs['cell_types']).union(set(primary.obs['cell_states'])))
node_colors = {
    'Unknown_primary': 'rgba(229, 190, 221, 0.8)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.8)',
    'Protein_metabolism-cell_death': 'rgba(189, 206, 244, 0.8)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.8)',
    'Immunoreactive_cells': 'rgba(169, 230, 227, 0.8)',
    'Epithelium_development': 'rgba(159, 217, 151, 0.8)',
    'Response_to_stress-ROS': 'rgba(197, 207, 132, 0.8)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.8)',
    'Smooth_muscle_cells_development': 'rgba(232, 186, 134, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_types'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['cell_types']]})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors[node] for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo = 'none'
    )
)])

# Update layout
fig.update_layout(title_text="Sankey Diagram of CAFs - Primary", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_gc/caf_primary.png')


#%%
## Plotting a sankey for ascites metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = ascites.obs.groupby(['cell_types', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(ascites.obs['cell_types']).union(set(ascites.obs['cell_states'])))
node_colors = {
    'Unknown_ascites': 'rgba(229, 190, 221, 0.8)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.8)',
    'RNA_metabolism': 'rgba(189, 206, 244, 0.8)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.8)',
    'Cellular_metabolism-ECM': 'rgba(169, 230, 227, 0.8)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.8)',
    'Extracellular_tissue_development': 'rgba(232, 186, 134, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_types'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['cell_types']]})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors[node] for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo = 'none'
    )
)])

# Update layout
fig.update_layout(title_text="Sankey Diagram of CAFs - Ascites", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_gc/caf_ascites.png')

#%%
## Plotting a sankey for metastasis metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = metastasis.obs.groupby(['cell_types', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(metastasis.obs['cell_types']).union(set(metastasis.obs['cell_states'])))
node_colors = {
    'Unknown_metastasis': 'rgba(229, 190, 221, 0.8)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.8)',
    'Protein_catabolism': 'rgba(189, 206, 244, 0.8)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.8)',
    'Immunoreactive_cells-T_cells': 'rgba(169, 230, 227, 0.8)',
    'Epithelium_development-cell_division': 'rgba(159, 217, 151, 0.8)',
    'Collagen_degradation': 'rgba(197, 207, 132, 0.8)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.8)',
    'Smooth_muscle_cells_development': 'rgba(232, 186, 134, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)',
    'Vascular_processes_regulation': 'rgba(211, 189, 173, 0.8)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_types'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['cell_types']]})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors[node] for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo = 'none'
    )
)])

# Update layout
fig.update_layout(title_text="Sankey Diagram of CAFs - Metastasis", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/home/marta.sallese/ov_cancer_atlas/atlas_project/plots_gc/caf_metastasis.png')
# %%
