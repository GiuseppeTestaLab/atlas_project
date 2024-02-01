import scanpy as sc
import numpy as np
import pandas as pd
import scgen

## Preprocessing
def preprocess_scgen(ad, batch, genes):
    #sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, batch_key=batch, n_top_genes=genes, inplace=True)
    sc.tl.pca(ad, use_highly_variable=True)
    sc.pp.neighbors(ad, use_rep='X_pca')
    sc.tl.umap(ad)
    return(ad)

def preprocess_scVI(adata, batch, genes):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=genes,
        batch_key=batch,
        subset=True
        )
    return(adata)