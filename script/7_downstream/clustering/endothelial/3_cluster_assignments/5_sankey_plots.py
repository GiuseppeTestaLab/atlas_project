

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
import plotly.io as pio
pio.kaleido.scope.chromium_args = tuple([arg for arg in pio.kaleido.scope.chromium_args if arg != "--disable-dev-shm-usage"])

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")

## Initialize directiories
tissueDir = rawPath + 'downstream/clustering/endothelial/'
figDir = figPath + ''

#%%
## Loading data
primary = sc.read_h5ad(tissueDir + 'adata_primary_embeddings_cellstates.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings_cellstates.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings_cellstates.h5ad')

#%%
## Plotting a sankey for primary metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = primary.obs.groupby(['treatment', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(primary.obs['treatment']).union(set(primary.obs['cell_states'])))
node_colors = {
    'Unknown_primary': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Immunoreactive_cells-T_cells': 'rgba(250, 224, 135, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Immunoreactive_cells-neutrophils': 'rgba(181, 220, 249, 0.6)',
    'Immunoreactive_cells': 'rgba(169, 230, 227, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'RNA_metabolism': 'rgba(210, 200, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'Angiogenesis': 'rgba(254, 168, 184, 0.6)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['treatment']]})

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
fig.update_layout(title_text="Sankey Diagram of Cell States - Primary", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image(figPath + 'endothelial_primary_treatment.png')

#%%
## Plotting a sankey for ascites metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = ascites.obs.groupby(['treatment', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(ascites.obs['treatment']).union(set(ascites.obs['cell_states'])))
node_colors = {
    'Unknown_ascites': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Immunoreactive_cells-T_cells': 'rgba(250, 224, 135, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Immunoreactive_cells-neutrophils': 'rgba(181, 220, 249, 0.6)',
    'Phagocytic_cells': 'rgba(169, 230, 227, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'RNA_metabolism': 'rgba(210, 200, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'Angiogenesis': 'rgba(254, 168, 184, 0.6)',
    'Immunoreactive_cells-B_cells': 'rgba(219, 217, 247, 0.6)',
    '13': 'rgba(255, 255, 255, 0.6)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['treatment']]})

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
fig.update_layout(title_text="Sankey Diagram of Cell States - Ascites", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image(figPath + 'endothelial_ascites_treatment.png')

#%%
## Plotting a sankey for metastasis metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = metastasis.obs.groupby(['treatment', 'cell_states']).size().reset_index(name='count')

# Define nodes
nodes = list(set(metastasis.obs['treatment']).union(set(metastasis.obs['cell_states'])))
node_colors = {
    'Unknown_metastasis': 'rgba(229, 190, 221, 0.6)',
    'Cellular_metabolism': 'rgba(247, 216, 194, 0.6)',
    'Immunoreactive_cells-T_cells': 'rgba(250, 224, 135, 0.6)',
    'Cycling_cells': 'rgba(255, 200, 194, 0.6)',
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Immunoreactive_cells-neutrophils': 'rgba(181, 220, 249, 0.6)',
    'Immunoreactive_cells': 'rgba(169, 230, 227, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'RNA_metabolism': 'rgba(210, 200, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'Angiogenesis': 'rgba(254, 168, 184, 0.6)',
    'Immunoreactive_cells-B_cells': 'rgba(219, 217, 247, 0.6)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_states'])
    links.append({"source": source_index, "target": target_index, "value": row['count'], "color": node_colors[row['treatment']]})

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
fig.update_layout(title_text="Sankey Diagram of Cell States - Metastasis", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image(figPath + 'endothelial_metastasis_treatment.png')
# %%
