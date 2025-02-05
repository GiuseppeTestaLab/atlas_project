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
scripts = ["primary_sankey", "ascites_sankey", "metastasis_sankey"]

# to be changed in relative path 
initialPath = scriptsPath + "7_downstream/clustering/endothelial/2_sankey/"

# this will run the python files to generate gene ontologies for each tissue

for script in scripts:
  path=initialPath+script
  print(path)
  sbatch="sbatch {}".format(path + '.sh')
  subprocess.check_output(sbatch, shell = True)
  

# %%
