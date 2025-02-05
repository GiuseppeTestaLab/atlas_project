## Imports
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import kaleido
import matplotlib
import plotly.express as px
import plotly.graph_objects as go
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")

## Initialize directiories
tissueDir = rawPath + '
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/'

#%%
## Loading data
primary = sc.read_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

## Plotting a sankey for primary metacells
# %%
# Generate sankey_data with treatment, cell_states, and cell_types
sankey_data = primary.obs.groupby(['treatment', 'cell_states', 'cell_types']).size().reset_index(name='count')

# Define nodes
treatments = list(sankey_data['treatment'].unique())
cell_states = list(sankey_data['cell_states'].unique())
cell_types = list(sankey_data['cell_types'].unique())

# Combine all nodes
nodes = treatments + cell_states + cell_types

# Node colors (based on your provided color scheme, extended for treatments)
node_colors = {
    'Unknown_primary': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Protein_metabolism-cell_death': 'rgba(143, 140, 254, 0.4)',
    'Epithelium_development': 'rgba(161, 237, 206, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.6)',
    'Immunoreactive_cells': 'rgba(169, 230, 227, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'Response_to_stress-ROS': 'rgba(197, 207, 132, 0.6)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'Smooth_muscle_cells_development': 'rgba(219, 217, 247, 0.6)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)'
}

# Define links
links = []

# Add links from treatments to cell states
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['treatment'], 'rgba(128, 128, 128, 0.6)')})

# Add links from cell states to cell types
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_states'])
    target_index = nodes.index(row['cell_types'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['cell_states'], 'rgba(128, 128, 128, 0.6)')})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors.get(node, 'rgba(128, 128, 128, 0.6)') for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo='none'
    )
)])

fig.update_layout(title_text="Tripartite Sankey Diagram", font_size=10, width=800, height=600)
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/caf_primary_triple.png')
# %%
## Plotting a sankey for metastasis metacells
# %%
# Generate sankey_data with treatment, cell_states, and cell_types
sankey_data = metastasis.obs.groupby(['treatment', 'cell_states', 'cell_types']).size().reset_index(name='count')

# Define nodes
treatments = list(sankey_data['treatment'].unique())
cell_states = list(sankey_data['cell_states'].unique())
cell_types = list(sankey_data['cell_types'].unique())

# Combine all nodes
nodes = treatments + cell_states + cell_types

# Node colors (based on your provided color scheme, extended for treatments)
node_colors = {
    'Unknown_metastasis': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'Protein_catabolism': 'rgba(189, 206, 244, 0.6)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.6)',
    'Immunoreactive_cells-T_cells': 'rgba(169, 230, 227, 0.6)',
    'Epithelium_development-cell_division': 'rgba(159, 217, 151, 0.6)',
    'Collagen_degradation': 'rgba(197, 207, 132, 0.6)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.6)',
    'Smooth_muscle_cells_development': 'rgba(232, 186, 134, 0.6)',
    'Vascular_processes_regulation': 'rgba(211, 189, 173, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)',
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)'
}

# Define links
links = []

# Add links from treatments to cell states
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['treatment'], 'rgba(128, 128, 128, 0.6)')})

# Add links from cell states to cell types
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_states'])
    target_index = nodes.index(row['cell_types'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['cell_states'], 'rgba(128, 128, 128, 0.6)')})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors.get(node, 'rgba(128, 128, 128, 0.6)') for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo='none'
    )
)])

fig.update_layout(title_text="Tripartite Sankey Diagram", font_size=10, width=800, height=600)
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/caf_metastasis_triple.png')

## Plotting a sankey for ascites metacells
# %%
# Generate sankey_data with treatment, cell_states, and cell_types
sankey_data = ascites.obs.groupby(['treatment', 'cell_states', 'cell_types']).size().reset_index(name='count')

# Define nodes
treatments = list(sankey_data['treatment'].unique())
cell_states = list(sankey_data['cell_states'].unique())
cell_types = list(sankey_data['cell_types'].unique())

# Combine all nodes
nodes = treatments + cell_states + cell_types

# Node colors (based on your provided color scheme, extended for treatments)
node_colors = {
    'Unknown_ascites': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'RNA_metabolism': 'rgba(189, 206, 244, 0.6)',
    'Angiogenesis': 'rgba(181, 220, 249, 0.6)',
    'Cellular_metabolism-ECM': 'rgba(169, 230, 227, 0.6)',
    'ECM_shaping_cells': 'rgba(210, 200, 151, 0.6)',
    'Extracellular_tissue_development': 'rgba(232, 186, 134, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'vCAF': 'rgba(161, 237, 206, 0.8)',
    'mCAF': 'rgba(250, 224, 135, 0.8)',
    'iCAF': 'rgba(254, 168, 184, 0.8)',
    'starCAF': 'rgba(219, 217, 247, 0.8)'
}

# Define links
links = []

# Add links from treatments to cell states
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['treatment'], 'rgba(128, 128, 128, 0.6)')})

# Add links from cell states to cell types
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['cell_states'])
    target_index = nodes.index(row['cell_types'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors.get(row['cell_states'], 'rgba(128, 128, 128, 0.6)')})

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=nodes,
        color=[node_colors.get(node, 'rgba(128, 128, 128, 0.6)') for node in nodes]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        color=[link["color"] for link in links],
        hoverinfo='none'
    )
)])

fig.update_layout(title_text="Tripartite Sankey Diagram", font_size=10, width=800, height=600)
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/fibroblasts/caf_ascites_triple.png')
# %%
