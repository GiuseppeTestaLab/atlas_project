#%%
import scanpy as sc
import sys
from sklearn.neighbors import KNeighborsClassifier
import pandas as pd


import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")
figPath = config.get("DEFAULT", "figPath")

major_cell_types = ["fibroblasts", "endothelial", "cancer"]
cell_type_names = {"fibroblasts": "FibroblastsMSK", "immune": "ImmuneMSK", "endothelial": "EndothelialMSK", "cancer": "CancerMSK"}

for major_cell_type in major_cell_types:
    initDir = rawPath + f'metacells_step0/{major_cell_type}/'
    outDir = rawPath + f'integration/metacells/{major_cell_type}_testing_ovca/'
    ooseDir = rawPath + f'out_of_sample_extension/{major_cell_type}_testing_ovca/'
    genes = scriptsPath + f'4_hdg/Tables/atlas_hdg_dispersion_patients_{major_cell_type}.csv'

    utilsPath = config.get("DEFAULT", "utilsPath")
    rawPath = config.get("DEFAULT", "rawPath")
    scriptsPath = config.get("DEFAULT", "scriptsPath")

    tissues = ["metastasis", "ascites", "primary"]
    for tissue in tissues:
    
        adata_path = ooseDir +  f'integrated_query_seacells_scarches_tissuetreat_{major_cell_type}_{tissue}_raw.h5ad'
        mcDir = f'/group/testa/Project/OvarianAtlas/atlas_project/raw_data/downstream_backup/downstream/clustering/{major_cell_type}/'
        
        ## Setting fig parameteres
        sc.settings.set_figure_params(dpi_save=300, frameon=False, format='png')

        
        adata = sc.read_h5ad(adata_path)
        adata


        #primary = sc.read_h5ad(mcDir + 'adata_primary_embeddings_cellstates.h5ad')
        #ascites = sc.read_h5ad(mcDir + 'adata_ascites_embeddings_cellstates.h5ad')
        tiusse_adata = sc.read_h5ad(mcDir + f'adata_{tissue}_embeddings_cellstates.h5ad')

        
        # Extract cell states
        #primary_states = primary.obs[['cell_states']].rename(columns={'cell_states': 'primary_state'})
        #ascites_states = ascites.obs[['cell_states']].rename(columns={'cell_states': 'ascites_state'})
        # metastasis_states = metastasis.obs[['cell_states']].rename(columns={'cell_states': 'metastasis_state'})

        
        # Merge with adata.obs
        #adata.obs = adata.obs.join(primary_states, how='left').join(ascites_states, how='left').join(metastasis_states, how='left')

        
        # Combine all categories
        # all_categories = pd.api.types.union_categoricals([
        #     adata.obs['primary_state'],
        #     adata.obs['ascites_state'],
        #     adata.obs['metastasis_state']
        # ]).categories

        # # Align the categories of 'ascites_state' and 'metastasis_state' to match the combined categories
        # adata.obs['primary_state'] = adata.obs['primary_state'].cat.set_categories(all_categories)
        # adata.obs['ascites_state'] = adata.obs['ascites_state'].cat.set_categories(all_categories)
        # adata.obs['metastasis_state'] = adata.obs['metastasis_state'].cat.set_categories(all_categories)


        
        # Create a final 'cell_states' column combining all sources, prioritizing primary > ascites > metastasis
        adata.obs_names = ["-".join(n.split("-")[0:2]) if n.startswith(n) else n for n in adata.obs_names]
        adata.obs['cell_states'] = tiusse_adata.obs['cell_states']


        
        sc.pl.umap(adata, color='cell_states')

        
        sc.pp.neighbors(adata, use_rep="latent_corrected")
        sc.tl.umap(adata)
        sc.pl.umap(adata, color='cell_states')



        
        adata.obs
        adata.obs["ref"] = ~adata.obs_names.str.startswith("new")
        adata_query = adata[adata.obs_names.str.startswith("new")]
        adata_ref = adata[~adata.obs_names.str.startswith("new")]
        


        sc.pl.umap(adata, color='cell_states')

        
        sc.pp.neighbors(adata, use_rep="latent_corrected")
        sc.tl.umap(adata)
        sc.pl.umap(adata, color='cell_type')

        
        # Extract latent space embeddings
        X_train = adata_ref.obsm["latent_corrected"]
        y_train = adata_ref.obs["cell_states"].to_numpy()

        X_test = adata_query.obsm["latent_corrected"]

        
        # Train kNN classifier
        knn = KNeighborsClassifier(n_neighbors=5, metric="cosine")  # Use cosine similarity for high-dimensional data
        knn.fit(X_train, y_train)

        
        # Predict labels for new dataset
        adata_query.obs["predicted_cell_states"] = knn.predict(X_test)

        
        sc.pl.umap(adata_query, color=["predicted_cell_states"], frameon=False, save=f"_oose_predicted_cell_states_{major_cell_type}_{tissue}.png")

        # 
        # from sklearn.metrics import confusion_matrix
        # confusion_matrix(adata_query.obs["cell_states"], adata_query.obs["predicted_cell_states"])

        adata.obs["predicted_cell_states"] = adata_query.obs["predicted_cell_states"]

        
        # Save adata
        adata.write_h5ad(ooseDir + f'integrated_query_seacells_scarches_tissuetreat_predicted_cellstates_{major_cell_type}_{tissue}_raw.h5ad')

        # %%
