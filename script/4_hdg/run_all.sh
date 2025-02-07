#!/bin/bash
sbatch cancer/1_dispersion_table_all_datasets_cancer.sh
sbatch cancer/2_hvg_union_patients_dispersion_cancer.sh
sbatch immune/1_dispersion_table_all_datasets_immune.sh
sbatch immune/2_hvg_union_patients_dispersion_immune.sh
sbatch endothelial/1_dispersion_table_all_datasets_endothelial.sh
sbatch endothelial/2_hvg_union_patients_dispersion_endothelial.sh
sbatch fibroblasts/1_dispersion_table_all_datasets_fibroblast.sh 
sbatch fibroblasts/2_hvg_union_patients_dispersion_fibroblast.sh