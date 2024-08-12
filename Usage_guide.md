# Ovarian Cancer Atlas (OvCA) and Single Cell transcriptomics Atlas Integration Pipeline (SCAIP)

This repository contains the code to reproduce SCAIP method used to generate OvCA described in the resource [paper]().  
This pipeline allows to integrate complex and very heterogenous cancer derived scRNA-seq data while preserving patients' heterogeneity.  

If just interested in accessing the atlas resource, you can reach it via a cellxgene web platform: [CellXGene](https://cellxgene.bioserver.ieo.it)  

![OvCa](https://github.com/GiuseppeTestaLab/atlas_project/blob/master/OvCA.png)

## Software requirements

Python 3.7 or higher.

## Documentation

You can build the documentation locally by following these steps:
```
git clone https://github.com/GiuseppeTestaLab/atlas_project/tree/master
```

## Files organization

We organized the scripts repository in folders following the order of the steps used to generate OvCA:
 - [1_original_counts](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/1_original_counts): contains all the scripts that allow to go from counts matrices to preprocessed data and harmonized metadata for every dataset of the atlas. You will find a `.csv` file used to set up the folder containing your data and the preprocessing parameteres. There are then two files that can be used to run all datasets together and then to concatenate the data.
 - [2_original_anndata](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/2_original_anndata): here on filtered and normalized data are computed the embeddings and `anndata var_names` are filtered to have common `var_names` to all datasets
 - [3_atlas_annotation](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/3_atlas_annotation)
- [4_hdg](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/4_hdg)





