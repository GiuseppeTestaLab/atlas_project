#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=varnames
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 

# Load configuration file
source ../../config.ini

# Set environment variables from the configuration file
scriptsPath=${DEFAULT_scriptsPath}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

module load singularity

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${scriptsPath}4_hdg/common_var_names.py"