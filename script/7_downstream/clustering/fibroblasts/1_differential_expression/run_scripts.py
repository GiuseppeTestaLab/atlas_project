#%%
#imports
import subprocess
import logging
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

scriptsPath = config.get("DEFAULT", "scriptsPath")

# put here the scripts that you want to run
#%%
scripts = ["primary_dea", "ascites_dea", "metastasis_dea"]

# to be changed in relative path 
initialPath = scriptsPath + "7_downstream/clustering/fibroblasts/1_differential_expression/"

# this will run the python files to generate gene ontologies for each tissue

for script in scripts:
  path=initialPath+script
  print(path)
  sbatch="sbatch {}".format(path + '.sh')
  subprocess.check_output(sbatch, shell = True)
  

# %%
