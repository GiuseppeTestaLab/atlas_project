#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=cpuq
#SBATCH --job-name=ascites
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j.log
module load singularity

singularity exec -B $bindPaths -H $homePath $image \
                 /bin/bash -c "eval \"\$(conda shell.bash hook)\" && conda activate downstream && python3 ${scriptsPath}7_downstream/clustering/endothelial/2_sankey/ascites_sankey.py"