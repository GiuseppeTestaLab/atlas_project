import pandas as pd

dir = '/home/marta.sallese/ov_cancer_atlas/atlas_project/script/4_hdg/Tables/'

hdg_old = pd.read_csv(dir + 'atlas_hdg_dispersion_patients_cancer.csv', index_col=0)
hdg_new = pd.read_csv(dir + 'atlas_hdg_dispersion_patients_cancer_zheng.csv', index_col=0)

hdg_new.index.isin(hdg_old.index).sum() # 2158

