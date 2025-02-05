# atlas total cell labelling strategy

# %%
import pandas as pd
import scanpy as sc
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")

sys.path.insert(1, utilsPath)

from cell_labeller import ( # type: ignore
    assign_scores,
    actual_labeller,
    create_cancer_adata,
    immune_labeller,
    fibroblast_labeller,
)

## Cell labelling

adata = sc.read(rawPath + "original_anndata/atlas_embeddings.h5ad")
adata.obs["dataset"] = adata.obs.paper_ID.str.split("_").str[0]

### Creating annotation columns with Vasquez-Garcia et al. signatures genes

adata = assign_scores(adata)

### Cell labelling strategy

adata = actual_labeller(adata)

### Removing columns that are not needed

adata.obs.drop(
    columns=[
        "ID",
        "sample_name",
        "patient_id",
        "cell_type",
        "cell_subtype",
        "sample_ID",
        "cell_labels_ratio",
        "assignment",
        "leiden-1.8",
    ],
    inplace=True,
)

### Creating cell labels for cell subtypes

adata_immune = immune_labeller(adata)

immune_cells = pd.DataFrame(
    adata_immune.obs["cell_types"], index=adata_immune.obs.index
)
immune_cells

adata.obs["cell_types"] = adata.obs["max"].astype("str")

adata.obs["cell_types"][adata.obs.index.isin(immune_cells.index)] = immune_cells[
    "cell_types"
]

adata_fibroblast = fibroblast_labeller(adata)

fibro_cells = pd.DataFrame(
    adata_fibroblast.obs["cell_types"], index=adata_fibroblast.obs.index
)
fibro_cells

adata.obs["cell_types"][adata.obs.index.isin(fibro_cells.index)] = fibro_cells[
    "cell_types"
]

adata.obs.cell_types.astype("category")

adata.write(rawPath + "atlas_annotated/atlas_embeddings_cell_labelled.h5ad")

# %%
