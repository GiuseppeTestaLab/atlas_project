# %%
# imports
import subprocess
import logging
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

# Get datasets and initial path from the configuration file
datasets = config.get("DATASETS", "datasets").split(", ")
scriptsPath = config.get("DEFAULT", "scriptsPath")

# %%
# this will run the python files to generate raw/processed files for each dataset

for dataset in datasets:
    path = scriptsPath + "1_original_counts/" + dataset
    print(path)
    sbatch = f"sbatch {path}/{dataset}.sh {dataset} {dataset[0:-4]}"
    print(sbatch)
    subprocess.check_output(sbatch, shell=True)


# %%
