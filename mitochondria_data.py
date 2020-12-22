import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import utils
import re

mat = pd.read_excel("../mitochondria/abundance.xlsx", sheet_name="Metabolomics", index_col=6).T
mat = mat.iloc[6:, :]
mat = mat.loc[:, ~mat.columns.duplicated(keep='first')]
groups = [i.split(".", 1)[0] for i in mat.index.tolist()]

mat['Group'] = groups
mat_selected_groups = mat.loc[mat['Group'].isin(['MB', 'FCCP'])]
# mat_selected_groups = mat.loc[mat['Group'].isin(['FCCP', 'Acti'])] 1 sig pathway

matrix_proc = utils.data_processing(mat_selected_groups.iloc[:, :-1], 0, 0)
print(matrix_proc.head)

utils.plot_PCA(matrix_proc.iloc[:, :-1], mat_selected_groups["Group"], "mitochondria", 10)

ttest_res = utils.t_tests(matrix_proc, mat_selected_groups["Group"], "fdr_bh")
DA_metabolites = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
print(len(DA_metabolites))

KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
background_list = matrix_proc.columns.tolist()

ORA_KEGG = utils.over_representation_analysis(DA_metabolites, background_list, KEGG_pathways)
print(ORA_KEGG)
print(len(ORA_KEGG[ORA_KEGG["P-adjust"] < 0.1]["P-adjust"].tolist()))
print(len(ORA_KEGG[ORA_KEGG["P-value"] < 0.1]["P-adjust"].tolist()))