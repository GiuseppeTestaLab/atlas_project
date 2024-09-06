#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=atlas
#SBATCH --mem=500GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${scriptsPath}2_original_anndata/2_run_atlas_preprocessing.py"

