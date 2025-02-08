#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=cpuq
#SBATCH --job-name=scvi_fibro
#SBATCH --mem=300GB
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
                 /bin/bash -c "source ~/.bashrc && mamba activate integration && python3 ${scriptsPath}6_integration/cells/fibroblasts/fibroblasts_cells_scVI_scANVI_2500.py"
