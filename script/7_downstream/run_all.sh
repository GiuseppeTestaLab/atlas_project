#!/bin/bash
sbatch LR_interactions/liana_cancer_endothelial.sh
sbatch LR_interactions/liana_cancer_fibroblasts.sh
sbatch LR_interactions/liana_cancer_immune.sh
sbatch LR_interactions/liana_endothelial_immune.sh
sbatch LR_interactions/liana_fibroblasts_endothelial.sh
sbatch LR_interactions/liana_fibroblasts_immune.sh