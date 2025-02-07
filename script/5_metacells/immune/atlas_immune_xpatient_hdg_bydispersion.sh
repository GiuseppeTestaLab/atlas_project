#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=gpuq
#SBATCH --job-name=immune
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log
#SBATCH --gres=gpu:2
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
                 /bin/bash -c "source ~/.bashrc && mamba activate seacells && python3 ${scriptsPath}5_metacells/immune/atlas_immune_xpatient_hdg_bydispersion.py"