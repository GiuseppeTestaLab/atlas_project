import pandas as pd
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

scriptsPath = config.get("DEFAULT", "scriptsPath")

dir = scriptsPath + '4_hdg/Tables/'

hdg_old = pd.read_csv(dir + 'atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
hdg_new = pd.read_csv(dir + 'atlas_hdg_dispersion_patients_cancer_zheng.csv', index_col=0)

hdg_new.index.isin(hdg_old.index).sum() # 2158

