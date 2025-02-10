#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=cpuq
#SBATCH --job-name=cluster_cancer
#SBATCH --mem=200GB
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

singularity exec -B $bindPaths -H $homePath $image \
                 /bin/bash -c "eval \"\$(conda shell.bash hook)\" && conda activate downstream && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/1_cancer_clusters_from_metacells.py && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/2_assigning_ontologies_to_clusters.py && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/3_plots_metacells.py && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/4_stats_plots.py && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/5_sankey_plots.py && \
                 python ${scriptsPath}7_downstream/clustering/cancer/3_cluster_assigments/6_global_embeddings.py"