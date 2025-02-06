#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=cpuq
#SBATCH --job-name=metastasis
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"/home/marta.sallese/miniconda3/envs/scvi_plots/bin/python /home/marta.sallese/ov_cancer_atlas/atlas_project/script/7_downstream/clustering/fibroblasts/2_sankey/metastasis_sankey.py"