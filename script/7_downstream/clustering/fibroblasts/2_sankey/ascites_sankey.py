# Visualizing leiden clusters with sankey plots

## Imports
#%%
import scanpy as sc
import pandas as pd
import numpy as np
import kaleido
import matplotlib
import plotly.express as px
import os
#kaleido cluster flag error fix see (https://github.com/plotly/Kaleido/issues/90)
import plotly.io as pio
pio.kaleido.scope.chromium_args = tuple([arg for arg in pio.kaleido.scope.chromium_args if arg != "--disable-dev-shm-usage"])


sc.logging.print_versions()
#%%
## Inizializing folders
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")

initDir = rawPath + 'downstream/clustering/fibroblasts/'
figDir = figPath + 'sankey/'
if not os.path.exists(figDir):
    os.makedirs(figDir)

## Loading data
#%%
adata = sc.read(initDir + 'adata_ascites_embeddings.h5ad')

#%%
leidenTotal=[]
for i in np.arange(0.01, 2.0, 0.1):
    sc.tl.leiden(adata,resolution = i,key_added="leiden-{}".format(round(i,2)))
    leidenTotal.append("leiden-{}".format(round(i,2)))

#%%
leidenTotal

#%%
LeidenDF = adata.obs[[i for i in adata.obs.columns if "leiden-" in i]]

#%%
LeidenDF["color"] = LeidenDF["leiden-0.91"].astype(int)

LeidenDF["treatment"] = adata.obs["treatment"]
LeidenDF["treatment"] = LeidenDF["treatment"].replace({"Naive":0, "NACT":1, "CHT":2})

LeidenDF["tumor_stage"] = adata.obs["tumor_stage"]
LeidenDF["tumor_stage"] = LeidenDF["tumor_stage"].replace({"IIIB":0, "IIIC":1, "IVA":2,"IVB":3})

LeidenDF["phase"] = adata.obs["phase"]
LeidenDF["phase"] = LeidenDF["phase"].replace({"G1":0, "G2M":1, "S":2})

LeidenDF["recurrence"] = adata.obs["recurrence"]
LeidenDF["recurrence"] = LeidenDF["recurrence"].replace({"Recurrence":0, "Unknown":1})

#%%
fig = px.parallel_categories(LeidenDF, 
                             dimensions=[i for i in adata.obs.columns if "leiden-" in i],
                             color_continuous_scale=px.colors.sequential.Agsunset,color="color",
                             width=2000,height=700)
fig.write_image(figDir + 'fibroblasts_ascites_leiden.png')
fig.write_html(figDir + 'fibroblasts_ascites_leiden.html')
fig.show()

#%%
fig = px.parallel_categories(LeidenDF, 
                             dimensions=[i for i in adata.obs.columns if "leiden-" in i],
                             color_continuous_scale=px.colors.sequential.Agsunset,color="treatment",
                             width=2000,height=700)
fig.write_image(figDir + 'fibroblasts_ascites_treatment.png')
fig.write_html(figDir + 'fibroblasts_ascites_treatment.html')
fig.show()

#%%
fig = px.parallel_categories(LeidenDF, 
                             dimensions=[i for i in adata.obs.columns if "leiden-" in i],
                             color_continuous_scale=px.colors.sequential.Agsunset,color="phase",
                             width=2000,height=700)
fig.write_image(figDir + 'fibroblasts_ascites_cellcycle.png')
fig.write_html(figDir + 'fibroblasts_ascites_cellcycle.html')
fig.show()

#%%
fig = px.parallel_categories(LeidenDF, 
                             dimensions=[i for i in adata.obs.columns if "leiden-" in i],
                             color_continuous_scale=px.colors.sequential.Agsunset,color="tumor_stage",
                             width=2000,height=700)
fig.write_image(figDir + 'fibroblasts_ascites_tumorstage.png')
fig.write_html(figDir + 'fibroblasts_ascites_tumorstage.html')
fig.show()

#%%
fig = px.parallel_categories(LeidenDF, 
                             dimensions=[i for i in adata.obs.columns if "leiden-" in i],
                             color_continuous_scale=px.colors.sequential.Agsunset,color="recurrence",
                             width=2000,height=700)
fig.write_image(figDir + 'fibroblasts_ascites_recurrence.png')
fig.write_html(figDir + 'fibroblasts_ascites_recurrence.html')
fig.show()
