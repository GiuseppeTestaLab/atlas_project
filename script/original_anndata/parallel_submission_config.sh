#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition=cpuq
#SBATCH --mem=30GB 
module load singularity

singularity run -B /group/testa -B /run/user -B $TMPDIR:/tmp \
-B /home/marta.sallese -H /home/marta.sallese/ov_cancer_atlas \
docker://testalab/downstream:covidiamo-3.1.0 \
"$1 $2"


"conda activate scvi \n
/home/marta.sallese/ov_cancer_atlas/bin/python \
/home/marta.sallese/ov_cancer_atlas/atlas_project/script/original_anndata/Config_file.py $1"

'conda activate scvi \\\n
echo $CONDA_PREFIX'