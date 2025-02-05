# %%

import os
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

scriptsPath = config.get("DEFAULT", "scriptsPath")

# to be changed in relative path
initialPath = (
    scriptsPath + "1_original_counts/"
)

datasets = [
    "Geistlinger2020",
    "Loret2022",
    "Olbrecht2021",
    "Qian2020",
    "Regner2021",
    "Ren2022",
    "Vasquez2022",
    "Xu2022",
    "Zhang2022",
]

for dataset in datasets:
    path = initialPath + dataset
    for name in os.listdir(path):
        initial_name = path + "/" + name
        final_name = path + "/" + dataset + "." + name.split(".")[-1]
        print("original file name {}".format(initial_name))
        print("renamed file name {}".format(final_name))
        if "sh" in final_name:
            os.rename(initial_name, final_name)
# %%
