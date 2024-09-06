#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=varnames
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 

module load singularity

# Check if dataset name is provided
if [ -z "$1" ]; then
  echo "No dataset name provided"
  exit 1
fi

dataset=$1

source config.ini

# Set environment variables from the configuration file
initialPath=${DEFAULT_initialPath}
scriptsPath=${initialPath}/${scriptsPath}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${scriptsPath}2_original_anndata/1_common_var_name.py"