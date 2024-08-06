

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
tissueDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream/clustering/immune/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/immune/'

#%%
## Loading data
primary = sc.read_h5ad(tissueDir + 'adata_primary_embeddings.h5ad')
metastasis = sc.read(tissueDir + 'adata_metastasis_embeddings.h5ad')
ascites = sc.read(tissueDir + 'adata_ascites_embeddings.h5ad')

#%%
## Plotting a sankey for primary metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = primary.obs.groupby(['treatment', 'cell_subtypes']).size().reset_index(name='count')

# Define nodes
nodes = list(set(primary.obs['treatment']).union(set(primary.obs['cell_subtypes'])))
node_colors = {
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'B_cells': 'rgba(128, 177, 211, 0.6)',
    'Dendritic_cells': 'rgba(251, 128, 114, 0.6)',
    'ILC': 'rgba(190, 186, 218, 0.6)',
    'M1_macrophages': 'rgba(253, 180, 98, 0.6)',
    'Mast_cells': 'rgba(179, 222, 105, 0.6)',
    'Myeloid_cells': 'rgba(252, 205, 229, 0.6)',
    'NK_CD56': 'rgba(217, 217, 217, 0.6)',
    'NK_cytotoxic': 'rgba(188, 128, 189, 0.6)',
    'Plasma_cells': 'rgba(204, 235, 197, 0.6)',
    'T_CD4_CXCL13': 'rgba(255, 237, 111, 0.6)',
    'T_CD4_naive': 'rgba(255, 174, 201, 0.6)',
    'T_CD4_reg': 'rgba(170, 240, 209, 0.6)',
    'T_CD8_CXCL13': 'rgba(203, 153, 201, 0.6)',
    'T_CD8_ISG': 'rgba(253, 192, 134, 0.6)',
    'T_CD8_cytotoxic': 'rgba(158, 218, 229, 0.6)',
    # 'T_cells': 'rgba(255, 165, 0, 0.6)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_subtypes'])
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
fig.update_layout(title_text="Sankey Diagram of Cell Types - Primary", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/immune/immune_primary_treatment.png')

#%%
## Plotting a sankey for ascites metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = ascites.obs.groupby(['treatment', 'cell_subtypes']).size().reset_index(name='count')

# Define nodes
nodes = list(set(ascites.obs['treatment']).union(set(ascites.obs['cell_subtypes'])))
node_colors = {
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'B_cells': 'rgba(128, 177, 211, 0.6)',
    'Dendritic_cells': 'rgba(251, 128, 114, 0.6)',
    'ILC': 'rgba(190, 186, 218, 0.6)',
    'M1_macrophages': 'rgba(253, 180, 98, 0.6)',
    'Mast_cells': 'rgba(179, 222, 105, 0.6)',
    'Myeloid_cells': 'rgba(252, 205, 229, 0.6)',
    'NK_CD56': 'rgba(217, 217, 217, 0.6)',
    'NK_cytotoxic': 'rgba(188, 128, 189, 0.6)',
    'Plasma_cells': 'rgba(204, 235, 197, 0.6)',
    'T_CD4_CXCL13': 'rgba(255, 237, 111, 0.6)',
    'T_CD4_naive': 'rgba(255, 174, 201, 0.6)',
    'T_CD4_reg': 'rgba(170, 240, 209, 0.6)',
    'T_CD8_CXCL13': 'rgba(203, 153, 201, 0.6)',
    'T_CD8_ISG': 'rgba(253, 192, 134, 0.6)',
    'T_CD8_cytotoxic': 'rgba(158, 218, 229, 0.6)',
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_subtypes'])
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
fig.update_layout(title_text="Sankey Diagram of Cell Subtypes - Ascites", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/immune/immune_ascites_treatment.png')

#%%
## Plotting a sankey for metastasis metacells

# Group by cell_types and cell_states and count occurrences
sankey_data = metastasis.obs.groupby(['treatment', 'cell_subtypes']).size().reset_index(name='count')

# Define nodes
nodes = list(set(metastasis.obs['treatment']).union(set(metastasis.obs['cell_subtypes'])))
node_colors = {
    'CHT': 'rgba(189, 206, 244, 0.6)',
    'Naive': 'rgba(159, 217, 151, 0.6)',
    'NACT': 'rgba(232, 186, 134, 0.6)',
    'B_cells': 'rgba(128, 177, 211, 0.6)',
    'Dendritic_cells': 'rgba(251, 128, 114, 0.6)',
    'ILC': 'rgba(190, 186, 218, 0.6)',
    'M1_macrophages': 'rgba(253, 180, 98, 0.6)',
    'Mast_cells': 'rgba(179, 222, 105, 0.6)',
    'Myeloid_cells': 'rgba(252, 205, 229, 0.6)',
    'NK_CD56': 'rgba(217, 217, 217, 0.6)',
    'NK_cytotoxic': 'rgba(188, 128, 189, 0.6)',
    'Plasma_cells': 'rgba(204, 235, 197, 0.6)',
    'T_CD4_CXCL13': 'rgba(255, 237, 111, 0.6)',
    'T_CD4_naive': 'rgba(255, 174, 201, 0.6)',
    'T_CD4_reg': 'rgba(170, 240, 209, 0.6)',
    'T_CD8_CXCL13': 'rgba(203, 153, 201, 0.6)',
    'T_CD8_ISG': 'rgba(253, 192, 134, 0.6)',
    'T_CD8_cytotoxic': 'rgba(158, 218, 229, 0.6)'
}

# Define links
links = []
for index, row in sankey_data.iterrows():
    source_index = nodes.index(row['treatment'])
    target_index = nodes.index(row['cell_subtypes'])
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
fig.update_layout(title_text="Sankey Diagram of Cell Subtypes - Metastasis", width=800, height=600)

# Show plot
fig.show()

# Save plot
fig.write_image('/group/testa/Project/OvarianAtlas/atlas_project/plots_def/cluster_assignments/immune/immune_metastasis_treatment.png')
# %%
