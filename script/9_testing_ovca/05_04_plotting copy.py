#%%
import seaborn as sns
import pandas as pd
import os
# %%
folder = "enrichment/"
files = os.listdir(folder)
files = [file for file in files if file.endswith(".csv")]
# %%
results = [pd.read_csv(f"{folder}{file}").assign(file=file) for file in files]
results = pd.concat(results)
# %%
results["major_cell_type"] = results["file"].apply(lambda x: x.split("_")[0])

#%%
# def getTissue(file):
#     if "type" in file:
#         return pd.NA
#     return file.split("_")[1]
results["tissue"] = results["cluster"].apply(lambda x: x.split("_")[0])
# %%
results["genes"] = results["file"].apply(lambda x: x.split("_")[-2])
results["is_tissue"] = results["file"].apply(lambda x: True if "tissue" in x else False)
#results["cluster"] = results["cluster"].apply(lambda x: "_".join(x.split("_")[1:]))
results = results.drop(columns=["Unnamed: 0"])
#%%
results["origin"] = results.apply(lambda x: x["genes"] + "_" + ("type" if x["is_tissue"] else "tissue"), axis=1)
results["is_raw"] = results["file"].apply(lambda x: True if "raw" in x else False)
results["is_sea"] = results["file"].apply(lambda x: True if "sea" in x else False)
#results["sea"] = results["file"].apply(lambda x: True if "sea" in x else False)
# %%
nonZero = results[results["percentage"] != 0]
# %%
sns.boxplot(data=nonZero, y="percentage", x="origin", hue="is_raw")
# %%
sns.scatterplot(data=nonZero, y="percentage", x="cluster", hue="is_raw")
# %%
nonRaw = results[results["is_raw"] == False]
sns.boxplot(data=nonRaw, y="percentage", x="origin")

# %%
sea = results[results["is_sea"] == True]
sns.boxenplot(data=sea, y="percentage", x="tissue")
sns.boxenplot(data=sea, y="percentage", x="origin", hue="major_cell_type")

# %%
ovca_type = results[results["origin"] == "ovca_type"]
sns.boxenplot(data=ovca_type, y="percentage", x="tissue")
#%%
sns.boxenplot(data=ovca_type, y="percentage", x="major_cell_type", hue="tissue")

# %%
sns.boxplot(data=sea, y="percentage", x="origin")
# %%
sns.scatterplot(data=sea, y="percentage", x="cluster", hue="origin")
# %%
sns.boxplot(data=sea, y="percentage", hue="origin")
# %%
nonZeroOvca = ovca_type[ovca_type["percentage"] != 0]
# %%
sns.boxenplot(data=nonZeroOvca, y="percentage", x="tissue")
#%%
sns.boxenplot(data=nonZeroOvca, y="percentage", x="major_cell_type", hue="tissue")
# %%
sns.boxplot(data=nonZeroOvca, y="percentage", x="origin")
# %%
ovcaZero = ovca_type[ovca_type["percentage"] == 0]
# %%
degs_query = pd.read_csv("enrichment/degs/endothelial/query_raw_tissue_ascites_Angiogenesis_de_ovca_tissue.csv")
from gprofiler import GProfiler
def run_gprof(query):
    gp = GProfiler(return_dataframe=True)
    enrichment_results = gp.profile(
        organism='hsapiens', 
        query=query,
        no_evidences=False, 
        sources=['GO:CC', 'GO:BP', 'GO:MF', 'REAC', 'KEGG'])
    return enrichment_results

enrichment = run_gprof(degs_query["names"].tolist())

# %%
import scanpy as sc
ad = sc.read_h5ad("/group/testa/Project/OvarianAtlasTestStep0/OvCA_umap_tiled.h5ad")
def run_gprof_bg(query, background):
    gp = GProfiler(return_dataframe=True)
    enrichment_results = gp.profile(
        organism='hsapiens', 
        query=query,
        background=background,
        no_evidences=False, 
        sources=['GO:CC', 'GO:BP', 'GO:MF', 'REAC', 'KEGG'])
    return enrichment_results

enrichment_bg = run_gprof_bg(degs_query["names"].tolist(), ad.var_names.tolist())

# %%
