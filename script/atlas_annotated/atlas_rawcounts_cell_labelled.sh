#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=cpuq
#SBATCH --job-name=atlas_def
#SBATCH --mem=400GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"/bin/python3 /home/marta.sallese/ov_cancer_atlas/atlas_project/script/atlas_annotated/atlas_rawcounts_cell_labelled.py"