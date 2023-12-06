#%%
#imports
import subprocess
import logging


# full list of implemented dataset are:
# ["Geistlinger2020",
# "Loret2022",
# "Olbrecht2021",
# "Qian2020",
# "Regner2021",
# "Ren2022",
# "Vasquez2022",
# "Xu2022",
# "Zhang2022"]

# put here the dataset that you want to rerun

datasets=["Geistlinger2020",
            "Loret2022",
            "Olbrecht2021",
            "Qian2020",
            "Regner2021",
            "Ren2022",
            "Vasquez2022",
            "Xu2022",
            "Zhang2022"]


# to be changed in relative path 
initialPath="/home/marta.sallese/ov_cancer_atlas/atlas_project/script/1_original_counts/"
#%%
# this will run the python files to generate raw/processed files for each dataset


for dataset in datasets:
  path=initialPath+dataset
  print(path)
  sbatch="sbatch {}".format(path + '/' + dataset + '.sh')
  subprocess.check_output(sbatch, shell = True)
  
