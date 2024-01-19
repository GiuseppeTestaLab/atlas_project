#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=cpuq
#SBATCH --job-name=seacells_scgen
#SBATCH --mem=200GB
#SBATCH --mail-type=ALL
#SBATCH --output=%x_%j.log 
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"/home/marta.sallese/ov_cancer_atlas/miniconda3/envs/scgen/bin/python /home/marta.sallese/ov_cancer_atlas/atlas_project/script/integration/metacells/Metacells_hdg_new/cancer_seacells_HDGbydispersion_scGen_tissuetreatment.py"