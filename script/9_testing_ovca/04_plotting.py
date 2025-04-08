#%%
import seaborn as sns
import pandas as pd
# %%
paths = ["zheng_cancer_overlap_filter.csv", "zheng_endothelial_overlap_filter.csv", "zheng_fibroblasts_overlap_filter.csv"]
labels = ["cancer", "endothelial", "fibroblasts","cancer", "endothelial", "fibroblasts"]

# Read the CSV files and add a column indicating the origin
results = [pd.read_csv(path).assign(origin=label) for path, label in zip(paths, labels)]
# %%
can = results[0]
# %%
combined_results = pd.concat(results)
# %%
sns.boxplot(data=combined_results, y="percentage", x="origin")
# %%
dir = "enrichment"
files = ["cancer_overlap_filter_full_raw.csv",
    "endothelial_overlap_filter_full_raw.csv", 
    "fibroblasts_overlap_filter_full_raw.csv", 
    "cancer_overlap_filter_full_raw.csv", 
    "endothelial_overlap_filter_full.csv", 
    "fibroblasts_overlap_filter_full.csv"]
paths = [f"{dir}/{file}" for file in files]
paths.extend(["../9_testing/cancer_overlap_filter.csv",
             "../9_testing/endothelial_overlap_filter.csv", 
             "../9_testing/fibroblasts_overlap_filter.csv"])

labels = ["cancer", "endothelial", "fibroblasts","cancer", "endothelial", "fibroblasts", "cancer", "endothelial", "fibroblasts"]
raws = [True, True, True, False, False, False, False, False, False]
seaCell = [False, False, False, False, False, False, True, True, True]

# %%

results = [pd.read_csv(path).assign(origin=label, raw=raw, sea=sea) for path, label, raw, sea in zip(paths, labels, raws, seaCell)]

# %%
res = pd.concat(results)
# %%
sns.boxplot(data=res, y="percentage", x="origin", hue="raw")
# %%
sns.boxplot(data=res, y="percentage", x="cluster", hue="raw")

# %%
nonRaw = res[res["raw"] == False]
sns.boxplot(data=nonRaw, y="percentage", x="origin", hue="sea")

# %%
sns.scatterplot(data=nonRaw, y="percentage", x="origin", hue="sea")

# %%
sea = res[res["sea"] == True]
# %%
sns.boxplot(data=sea, y="percentage", x="origin")
# %%
sns.scatterplot(data=sea, y="percentage", x="cluster", hue="origin")
# %%
sns.boxplot(data=sea, y="percentage", hue="origin")
# %%
