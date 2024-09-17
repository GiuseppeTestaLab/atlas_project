#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=ren
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 

dataset=$1
datasetPy=$2

# Load configuration file
source ../../config.ini

# Set environment variables from the configuration file
scriptsPath=${DEFAULT_scriptsPath}
datasetPath=${scriptsPath}/2_original_anndata/${dataset}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

module load singularity

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${datasetPath}/${datasetPy}.py"