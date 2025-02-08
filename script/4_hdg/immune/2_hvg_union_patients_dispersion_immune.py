# Computing HDG by patient based on dispersion values

# %%
import pandas as pd
import configparser

# Read configuration file
config = configparser.ConfigParser()
config.read("../../utils/config.ini")


scriptsPath = config.get("DEFAULT", "scriptsPath")


dir = scriptsPath + "4_hdg/Tables/"

# %%
dispersion_table = pd.read_csv(dir + "dispersion_table_immune.csv", index_col=0)

# %%
dispersion_table

# %%
max_values_per_row = dispersion_table.max(axis=1)

print(max_values_per_row)

# %%
# filtered_df = dispersion_table[dispersion_table > 2.05].dropna(how='all')  # 5088 genes
filtered_df = dispersion_table[dispersion_table > 1.88].dropna(how="all")  # 5074 genes

# %%
list_of_genes = filtered_df.index.tolist()
print(len(list_of_genes))

# %%
hdg_table = pd.DataFrame(index=dispersion_table.index)

# %%
hdg_table["highly_variable"] = hdg_table.index.isin(list_of_genes)

# %%
hdg_table.to_csv(dir + "atlas_hdg_common_dispersion_patients_immune.csv")

# %%
hdg_true = pd.DataFrame(index=list_of_genes)
hdg_true["highly_variable"] = "True"
hdg_true.to_csv(dir + "atlas_hdg_dispersion_patients_immune.csv")
# %%
