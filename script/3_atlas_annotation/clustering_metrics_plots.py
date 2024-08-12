#%%
import numpy as np
from sklearn import metrics
import scanpy as sc
import pandas as pd
import seaborn as sns

#%%
initDir = '/group/testa/Project/OvarianAtlas/atlas_project/raw_data/atlas_annotated/'
figDir = '/group/testa/Project/OvarianAtlas/atlas_project/plots_def/atlas_annotated/'

#%%
palette = {
    'Adjusted Rand Index': '#cd2867',            
    'Adjusted Mutual Information': '#8930f0',    
    'Fowlkes-Mallows Index': '#1bab91'           
}

#%%
data = ['clustering_metrics_cancer.csv', 'clustering_metrics_endothelial.csv', 'clustering_metrics_fibroblasts.csv', 'clustering_metrics_immune.csv']

for i in data:
    df = pd.read_csv(initDir + i)
    df.columns = ['Leiden', 'Adjusted Rand Index', 'Adjusted Mutual Information', 'Fowlkes-Mallows Index']
    df = df.set_index('Leiden')
    ax = sns.lineplot(data=df, palette=palette)
    ax.set(xlabel='Leiden Resolution', ylabel='Score')
    ax.set_title(f'Clustering Metrics for {i.split("_")[2]} Cells')
    plt.savefig(figDir + 'clustering_metrics_' + i.split('_')[2] + '.png')