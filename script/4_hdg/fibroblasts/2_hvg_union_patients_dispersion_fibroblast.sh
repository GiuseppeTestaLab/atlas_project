#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=hvg_pat_fibro
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

# Load configuration file
source config.ini

# Set environment variables from the configuration file
initialPath=${DEFAULT_initialPath}
scriptsPath=${initialPath}/${scriptsPath}
bindPaths=$(echo ${SINGULARITY_bindPaths} | tr ',' ' ')
homePath=${SINGULARITY_homePath}
image=${SINGULARITY_image}

singularity run -B $bindPaths -H $homePath $image \
"/bin/python3 ${scriptsPath}4_hdg/cancer/2_hvg_union_patients_dispersion_fibroblast.py"

# singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
# -B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
# docker://testalab/downstream:covidiamo-3.1.0 \
# "/bin/python3 /home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/fibroblasts/2_hvg_union_patients_dispersion_fibroblast.py"