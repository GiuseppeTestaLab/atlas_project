#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=geistlinger
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log

# Check if dataset name is provided
if [ -z "$1" ]; then
  echo "No dataset name provided"
  exit 1
fi

dataset=$1
datasetPy=$2

# Load configuration file
source ../../utils/bash_ini_parser/read_ini.sh
read_ini ../../utils/config.ini

# Set environment variables from the configuration file
scriptsPath=${INI__DEFAULT__scriptsPath}
datasetPath=${INI__DEFAULT__scriptsPath}1_original_counts/${dataset}
bindPaths=${INI__SINGULARITY__bindPaths}
bindPaths=$(eval echo $bindPaths)
homePath=${INI__SINGULARITY__homePath}
image=${INI__SINGULARITY__image}

echo script=${scriptsPath} dataset=${datasetPath} bind=${bindPaths} home=${homePath} image=${image}

module load singularity

singularity exec -B $bindPaths -H $homePath $image \
                 /bin/bash -c "source ~/.bashrc && mamba activate ovarian && python3 ${datasetPath}/${datasetPy}.py"

