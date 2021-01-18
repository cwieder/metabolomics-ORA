import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils

data = pd.read_excel("../COPD/Metabolomics_raw.xlsx", index_col=0, header=3).T

mat_proc = utils.data_processing(data, 18, 1)
labels = ["COPD" if i.startswith("C") else "Healthy" if i.startswith("H") else "Null" for i in data.index]
data['Groups'] = labels
#utils.plot_PCA(mat_proc, data.loc["C01":,"Groups"], "Hansbro COPD vs Healthy", 57)

ttest_res = utils.t_tests(mat_proc, data.loc["C01":,"Groups"], "fdr_bh")
DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()

print(len(DEM), "Differential metabolites")
# with open("../COPD/DEM.txt", "a") as infile:
#     for i in DEM:
#         infile.write(i+"\n")

KEGG_pathways = pd.read_csv("KEGG_reference_pathways_compounds.csv", dtype=str, index_col=0)
background_list = mat_proc.columns.tolist()
print(len(background_list))

KEGG_dict = data.loc["KEGG", :].to_dict()
KEGG_background = [KEGG_dict[k] for k, v in KEGG_dict.items() if str(v) != 'nan']
DEM_KEGG = [KEGG_dict[k] for k in DEM if KEGG_dict[k] != 'nan']
print(len(DEM_KEGG), "Differential metabolites mapping to KEGG pathways")
print(len(KEGG_background), "Background mapping to KEGG")

ora_res = utils.over_representation_analysis(DEM_KEGG, KEGG_background, KEGG_pathways)
print(ora_res)
print(len(ora_res[ora_res["P-value"] < 0.1]["P-adjust"].tolist()))
# ora_res.to_csv("../COPD/ORA.csv")
