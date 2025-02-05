# Functions to create gene ontolgies annotations

## Imports
import os
import logging
import shutil
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams
from gprofiler import GProfiler
import sys
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")

utilsPath = config.get("DEFAULT", "utilsPath")
rawPath = config.get("DEFAULT", "rawPath")
scriptsPath = config.get("DEFAULT", "scriptsPath")

sys.path.insert(1, utilsPath)
from plotting_bubble import scale_data_5_75, plot_enrich 


## Functions

def annotate_ontolgies(adata, directory_root, leidenTotal, dedf, log_file):
    errorDic={}
    
    # Check if the log exists
    if os.path.isfile(log_file):
    # If it doesn't exist, create it
        os.remove(log_file)

    logging.basicConfig(filename=log_file, level=logging.INFO) 
                        # format='%(primary)s %(levelname)s %(name)s %(message)s')
    logger=logging.getLogger(__name__)

    for lei in leidenTotal:
        errorDic[lei]={}
        directory = directory_root + lei
        # Check if the directory exists
        if not os.path.exists(directory):
            # If it doesn't exist, create it
            os.makedirs(directory)
        else:
            shutil.rmtree(directory)
            os.makedirs(directory)
        for cl in adata.obs[lei].unique():
            if len(adata.obs[lei].unique()) == 1:
                continue
            dedf[lei][cl] = sc.get.rank_genes_groups_df(adata, group=cl, key ='wilcoxon_'+lei)
            dedf[lei][cl].to_csv(directory + '/rank_gene_groups_df_' + cl + '.csv')
            try:
                # print(dedf[lei][cl])
                threshold1 = dedf[lei][cl].dropna(axis='rows')
                threshold2 = threshold1.loc[(threshold1['logfoldchanges'] > 1) & (threshold1['logfoldchanges'] < 100), :]
                threshold3 = threshold2.loc[threshold2['scores'] > 5, :]
                threshold4 = threshold3.loc[threshold3['pvals_adj'] < 0.05, :]
                # break
                print('{}_{}_{}_{}'.format(lei, cl, threshold4.shape[0], threshold4.scores.min()))
                # ontologia di (dedf[lei][cl])
                gp = GProfiler(return_dataframe=True)
                query = threshold4.names.to_list()
                ontology = gp.profile(organism='hsapiens', 
                                        query=query, 
                                        no_evidences=False, 
                                        background=adata.var_names.to_list(),
                                        sources=['GO:CC', 'GO:BP', 'GO:MF','REAC','KEGG'])
                ontology.to_csv(directory + '/gprofiler_' + cl + '.csv')
                if ontology.shape[0] > 0:
                    plot_enrich(ontology, filename=directory + '/ontology_' + cl + '.png')
                else:
                    logger.info('leiden {} cluster {}'.format(lei, cl))
                    logger.info('ontology empty')
                
                
                errorDic[lei][cl]={}

            except Exception as e:
                errorDic[lei][cl]=e
                # with open(directory + '/ontology_' + cl + '.log', 'w') as f:
                #     f.write('error')
                logger.info('leiden {} cluster {}'.format(lei, cl))
                logger.error(e)
                continue
    
    return(adata)
