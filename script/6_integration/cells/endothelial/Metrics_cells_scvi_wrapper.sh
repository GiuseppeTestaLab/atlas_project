#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=cpuq
#SBATCH --job-name=scib_endothelial
#SBATCH --mem=150GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"/home/marta.sallese/ov_cancer_atlas/miniconda3/envs/scvi/bin/python /home/marta.sallese/ov_cancer_atlas/atlas_project/script/6_integration/cells/endothelial/Metrics_cells_scvi_wrapper.py"