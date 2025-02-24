#%%
import scanpy as sc
import sys
import os

dest_dir = r"/group/testa/Project/OvarianAtlasTestStep0/obs"
path = sys.argv[1]
obs_path = os.path.join(dest_dir ,"_".join(path.split(os.path.sep)[4:]) + ".csv")
adata = sc.read_h5ad(path)
adata.obs.to_csv(obs_path)

