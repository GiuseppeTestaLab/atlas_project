#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=loret
#SBATCH --mem=200GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log

dataset=$1
datasetPy=$2


# Load configuration file
source ../../utils/bash_ini_parser/read_ini.sh
read_ini ../../utils/config.ini

# Set environment variables from the configuration file
scriptsPath=${INI__DEFAULT__scriptsPath}
datasetPath=${INI__DEFAULT__scriptsPath}2_original_anndata/${dataset}
bindPaths=${INI__SINGULARITY__bindPaths}
bindPaths=$(eval echo $bindPaths)
homePath=${INI__SINGULARITY__homePath}
image=${INI__SINGULARITY__image}

echo script=${scriptsPath} dataset=${datasetPath} bind=${bindPaths} home=${homePath} image=${image}

module load singularity

singularity exec -B $bindPaths -H $homePath $image \
                 /bin/bash -c "source ~/.bashrc && mamba activate ovarian && python3 ${datasetPath}/${datasetPy}.py"

