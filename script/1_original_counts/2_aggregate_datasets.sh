#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=aggregate
#SBATCH --mem=400GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

# Check if dataset name is provided
if [ -z "$1" ]; then
  echo "No dataset name provided"
  exit 1
fi

dataset=$1

# Load configuration file
source config.ini

# Set environment variables from the configuration file
initialPath=${DEFAULT_initialPath}
scriptsPath=${initialPath}/${scriptsPath}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

module load singularity

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${scriptsPath}1_original_counts/2_aggregate_datasets.py"