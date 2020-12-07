# Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19
# Goldman et al

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils
import pickle

goldman_mat = pd.read_excel("../Goldman_data/Goldman_metabolites.xlsx")
print(goldman_mat.shape)

goldman_mat_proc = utils.data_processing(goldman_mat, firstrow=0, firstcol=2)
print(goldman_mat_proc.shape)

goldman_clin_data = pd.read_excel("../Goldman_data/Goldman_clinical_data.xlsx")
# covid_sample_id =
# print(goldman_clin_data)

# utils.plot_PCA(goldman_mat_proc, goldman_mat["Healthy donor sample or COVID19 sample"], "COVID-19 vs healthy")

t_test_res = utils.t_tests(goldman_mat_proc, goldman_mat["Healthy donor sample or COVID19 sample"], "bonferroni")
t_test_res.to_csv("Goldman_ttest_res.csv")


KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)


Reactome_pathways = pd.read_csv("Reactome_pathway_set.csv", dtype=str, index_col=0)
Reactome_human_ids = [i for i in Reactome_pathways.index if i.startswith("R-HSA")]
Reactome_human = Reactome_pathways[Reactome_pathways.index.isin(Reactome_human_ids)]

name_map = pd.read_csv("../Goldman_data/name_map.csv", dtype=str)
background_list = goldman_mat_proc.columns.tolist()
background_list_KEGG = name_map[name_map["Query"].isin(background_list)]['KEGG'].dropna().tolist()
background_list_chEBI = name_map[name_map["Query"].isin(background_list)]['ChEBI'].dropna().tolist()

print(len(background_list_KEGG), "in background list KEGG")
DA_metabolites = t_test_res[t_test_res["P-adjust"] < 0.0000000000000001]["Metabolite"].tolist()

name_map[name_map["Query"].isin(background_list)]['KEGG'].dropna().to_csv("Goldman_background_KEGG.csv", index=None)
# t_test_res[t_test_res["P-adjust"] < 0.05]["Metabolite"].to_csv("Goldman_DA_names.csv", index=None)
# pd.DataFrame(goldman_mat_proc.columns.to_list()).to_csv("Goldman_background_names.csv", index=None)

print(len(DA_metabolites), "DA metabolites (all)")

DA_metabolites_KEGG = name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().tolist()
DA_metabolites_chEBI = name_map[name_map["Query"].isin(DA_metabolites)]['ChEBI'].dropna().tolist()
print(len(DA_metabolites_chEBI), "DA ChEBI")

name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().to_csv("Goldman_DA_KEGG.csv", index=False)
print(len(DA_metabolites_KEGG), "DA in KEGG")
# name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().to_csv("Goldman_DA_metabolites.csv", index=None)

# with open('KEGG_compounds.pickle', 'rb') as handle:
#     all_KEGG_compounds = pickle.load(handle)

ORA_res = utils.over_representation_analysis(DA_metabolites_KEGG, background_list_KEGG, KEGG_pathways)
ORA_res.to_csv("../Goldman_data/goldman_ORA.csv")

ORA_reactome = utils.over_representation_analysis(DA_metabolites_chEBI, background_list_chEBI, Reactome_human)
# print(ORA_reactome)
print("length")
print(len(ORA_reactome[ORA_reactome["P-adjust"] < 0.2]["P-adjust"].tolist()))
ORA_reactome.to_csv("Goldman_ORA_reactome.csv")
