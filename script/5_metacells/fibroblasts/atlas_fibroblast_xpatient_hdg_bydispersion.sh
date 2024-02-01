#!/bin/bash
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=cpuq
#SBATCH --job-name=fibroblast
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"/home/marta.sallese/ov_cancer_atlas/miniconda3/envs/seacells/bin/python /home/marta.sallese/ov_cancer_atlas/atlas_project/script/5_metacells/fibroblasts/atlas_fibroblast_xpatient_hdg_bydispersion.py"