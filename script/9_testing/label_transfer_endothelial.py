
#%%
import scanpy as sc
import sys
from sklearn.neighbors import KNeighborsClassifier

#%%
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'metacells_step0/endothelial/'
outDir = rawPath + 'integration/metacells/endothelial_testing/'
ooseDir = rawPath + 'out_of_sample_extension/endothelial_testing/'
genes = scriptsPath + '4_hdg/Tables/atlas_hdg_dispersion_patients_endothelial.csv'

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")


#%%
adata_path = ooseDir +  "integrated_query_seacells_scarches_tissuetreat.h5ad"
mcDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream_backup/downstream/clustering/endothelial/'

#%%
## Setting fig parameteres
sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

#%%
adata = sc.read_h5ad(adata_path)
adata

#%%
adata.obs
adata.obs["ref"] = ~adata.obs_names.str.startswith("new")
adata_query = adata[adata.obs_names.str.startswith("new")]
adata_ref = adata[~adata.obs_names.str.startswith("new")]
#%%
sc.pl.umap(adata, color='cell_type')

#%%
sc.pp.neighbors(adata, use_rep="latent_corrected")
sc.tl.umap(adata)
sc.pl.umap(adata, color='cell_type')

#%%
# Extract latent space embeddings
X_train = adata_ref.obsm["latent_corrected"]
y_train = adata_ref.obs["cell_type"].to_numpy()

X_test = adata_query.obsm["latent_corrected"]

#%%
# Train kNN classifier
knn = KNeighborsClassifier(n_neighbors=5, metric="cosine")  # Use cosine similarity for high-dimensional data
knn.fit(X_train, y_train)

#%%
# Predict labels for new dataset
adata_query.obs["predicted_cell_types"] = knn.predict(X_test)

#%%
sc.pl.umap(adata_query, color=["predicted_cell_types"], frameon=False, save="_oose_predicted_cell_states.png")

#%%
from sklearn.metrics import confusion_matrix
confusion_matrix(adata_query.obs["cell_type"], adata_query.obs["predicted_cell_types"])

# %%
adata.obs["predicted_cell_types"] = adata_query.obs["predicted_cell_types"]
# %%

#%%
# Save adata
adata.write_h5ad(ooseDir + 'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates.h5ad')

# %%
