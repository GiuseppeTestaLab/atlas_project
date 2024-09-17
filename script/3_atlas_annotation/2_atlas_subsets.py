import scanpy as sc
import pandas as pd

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")
rawPath = config.get("DEFAULT", "rawPath")

initDir = rawPath + "atlas_annotated/"

adata = sc.read(initDir + "atlas_embeddings_cell_labelled.h5ad")

# adata = adata.raw.to_adata()

## creating adata by cell type
#### cancer
adata_cancer = sc.AnnData(
    X=adata[adata.obs["max"] == "CancerMSK"].X,
    var=adata[adata.obs["max"] == "CancerMSK"].var,
    obs=adata[adata.obs["max"] == "CancerMSK"].obs,
)

adata_cancer.write_h5ad(initDir + "atlas_cancer_filt_norm_nolog.h5ad")

#### immune
adata_immune = sc.AnnData(
    X=adata[adata.obs["max"] == "HematopoieticMSK"].X,
    var=adata[adata.obs["max"] == "HematopoieticMSK"].var,
    obs=adata[adata.obs["max"] == "HematopoieticMSK"].obs,
)

adata_immune.write_h5ad(initDir + "atlas_immune_filt_norm_nolog.h5ad")

#### fibroblasts
adata_fibro = sc.AnnData(
    X=adata[adata.obs["max"] == "FibroblastsMSK"].X,
    var=adata[adata.obs["max"] == "FibroblastsMSK"].var,
    obs=adata[adata.obs["max"] == "FibroblastsMSK"].obs,
)

adata_fibro.write_h5ad(initDir + "atlas_fibroblast_filt_norm_nolog.h5ad")

#### endothelial
adata_endo = sc.AnnData(
    X=adata[adata.obs["max"] == "EndothelialMSK"].X,
    var=adata[adata.obs["max"] == "EndothelialMSK"].var,
    obs=adata[adata.obs["max"] == "EndothelialMSK"].obs,
)

adata_endo.write_h5ad(initDir + "atlas_endothelial_filt_norm_nolog.h5ad")
