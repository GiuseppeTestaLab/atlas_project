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

Every script comes with a `.sh` file that specifies the conda env to use to run the script.
We organized the scripts repository in folders following the order of the steps used to generate OvCA:
 - [1_original_counts](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/1_original_counts): contains all the scripts that allow to go from counts matrices to preprocessed data and harmonized metadata for every dataset of the atlas. You will find a `.csv` file used to set up the folder containing your data and the preprocessing parameteres. There are then two files that can be used to run all datasets together and then to concatenate the data.
 - [2_original_anndata](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/2_original_anndata): here on filtered and normalized data are computed the embeddings and `anndata var_names` are filtered to have common `var_names` to all datasets
 - [3_atlas_annotation](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/3_atlas_annotation): here we annotated the major cell types and divided the atlas into the four subsets
 - [4_hdg](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/4_hdg): definition of the highly defining genes. You can use the `.csv` file to set your working directories.
 - [5_metacells]I(https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/5_metacells): computation of metacells per cell type in the HDG space
 - [6_integration](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/6_integration): metacells and cells integration in the HDG space with scGen
 - [7_downstream](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/7_downstream): divided bt cell type, here we firstly divide the dataset by the three main tissue types (primary, ascites, metastasis), perform clusteirng and [differential expression analysis](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/7_downstream/clustering/cancer/1_differential_expression), [cluster stability analysis](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/7_downstream/clustering/cancer/2_sankey), and then [cluster assignment](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/7_downstream/clustering/cancer/3_cluster_assignments) defined by cell ontologies. Then we mapped the cell states defined in metacells to parental cells and performed [cell to cell interaction analyses](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/7_downstream/LR_interactions).
 - [8_out_of_sample_extension](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/8_out_of_sample_extension): here we used scArches to perform refernce mapping of a new dataset on the atlas. For the new dataset firstly we define [HDGs](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/8_out_of_sample_extension/Zheng_hdg), derive [metacells per cell type](https://github.com/GiuseppeTestaLab/atlas_project/tree/master/script/8_out_of_sample_extension/Zheng_metacells) and then used the model generated for atlas integration to integrate the new dataset.

