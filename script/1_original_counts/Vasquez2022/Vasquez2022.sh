#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=vasquez
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log

# Check if dataset name is provided
if [ -z "$1" ]; then
  echo "No dataset name provided"
  exit 1
fi

dataset=$1
datasetPy=$2

# Load configuration file
source config.ini

# Set environment variables from the configuration file
scriptsPath=${DEFAULT_scriptsPath}
datasetPath=${scriptsPath}/1_original_counts/${dataset}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

module load singularity

singularity run --nv -B $bindPaths -H $homePath $image \
                 "source ~/.bashrc && mamba activate ovarian && python3 ${datasetPath}/${datasetPy}.py