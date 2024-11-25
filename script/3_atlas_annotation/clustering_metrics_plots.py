# %%
import numpy as np
from sklearn import metrics
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
figPath = config.get("DEFAULT", "figPath")
# %%
initDir = rawPath + "atlas_annotated/"
figDir = figPath + "atlas_annotated/"

# %%
palette = {
    "Adjusted Rand Index": "#cd2867",
    "Adjusted Mutual Information": "#8930f0",
    "Fowlkes-Mallows Index": "#1bab91",
}

# %%
data = [
    "clustering_metrics_cancer.csv",
    "clustering_metrics_endothelial.csv",
    "clustering_metrics_fibroblasts.csv",
    "clustering_metrics_immune.csv",
]

# %%
all_data = []
# Load all data and find the global y-axis limits
for i in data:
    df = pd.read_csv(initDir + i)
    df.columns = [
        "Leiden",
        "Adjusted Rand Index",
        "Adjusted Mutual Information",
        "Fowlkes-Mallows Index",
    ]
    all_data.append(df.set_index("Leiden"))
# %%
# Find the global y-axis limits
ymin = min([df.min().min() for df in all_data])
ymax = max([df.max().max() for df in all_data])
# %%
# Make plot with the same y-axis limits
for i, df in enumerate(all_data):
    ax = sns.lineplot(data=df, palette=palette)
    ax.set(xlabel="Leiden Resolution", ylabel="Score")
    ax.set_ylim(ymin, ymax)
    output_filename = os.path.join(figDir, data[i].replace(".csv", "_plot_consistent_axis.pdf"))
    ax.figure.savefig(output_filename)
    ax.figure.clear()
