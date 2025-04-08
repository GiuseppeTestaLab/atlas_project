#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

#SBATCH --partition=cpuq
#SBATCH --job-name=fibroblasts_testing
#SBATCH --mem=128GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log
module load singularity

# Load configuration file
source ../../utils/bash_ini_parser/read_ini.sh
read_ini ../../utils/config.ini

# Set environment variables from the configuration file
scriptsPath=${INI__DEFAULT__scriptsPath}
bindPaths=${INI__SINGULARITY__bindPaths}
bindPaths=$(eval echo $bindPaths)
homePath=${INI__SINGULARITY__homePath}
image=${INI__SINGULARITY__image}

singularity exec --nv -B $bindPaths -H $homePath $image \
                 /bin/bash -c "eval \"\$(conda shell.bash hook)\" && conda activate scarches && \
                 python3 ${scriptsPath}9_testing/01_scarches_fibroblasts.py && \
                 python3 ${scriptsPath}9_testing/02_label_transfer_fibroblasts.py && \
                 conda activate downstream && \
                 python3 ${scriptsPath}9_testing/03_cell_states_fibroblasts.py"
