#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpuq
#SBATCH --job-name=obs_to_pandas
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log

module load singularity

# Load configuration file
source ./bash_ini_parser/read_ini.sh
read_ini ./config.ini

# Set environment variables from the configuration file
scriptsPath=${INI__DEFAULT__scriptsPath}
bindPaths=${INI__SINGULARITY__bindPaths}
bindPaths=$(eval echo $bindPaths)
homePath=${INI__SINGULARITY__homePath}
image=${INI__SINGULARITY__image}

singularity exec -B $bindPaths -H $homePath $image \
                 /bin/bash -c "eval \"\$(conda shell.bash hook)\" && conda activate downstream && python3 ${scriptsPath}../utils/to_pandas.py ${1}"