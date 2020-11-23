import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils

mat = pd.read_excel("../Haiwei/Haiwei_no_annotation_CRC_healthy.xlsx", index_col=0)

mat_proc = utils.data_processing(mat, 0, 1)

# utils.plot_PCA(mat_proc, mat['Groups'], "Colorectal cancer vs. healthy")
classes = mat["Groups"].to_dict()
ttest_res = utils.t_tests(mat_proc, mat["Groups"], "fdr_bh")

KEGG_pathways = pd.read_csv("KEGG_reference_pathways_compounds.csv", dtype=str, index_col=0)

Reactome_pathways = pd.read_csv("Reactome_pathway_set.csv", dtype=str, index_col=0)
Reactome_human_ids = [i for i in Reactome_pathways.index if i.startswith("R-HSA")]
Reactome_human = Reactome_pathways[Reactome_pathways.index.isin(Reactome_human_ids)]

name_map = pd.read_csv("../Haiwei/name_map.csv", dtype=str)

background_list = mat_proc.columns.tolist()
background_list_chEBI = name_map[name_map["Query"].isin(background_list)]['ChEBI'].dropna().tolist()
background_list_KEGG = name_map[name_map["Query"].isin(background_list)]['KEGG'].dropna().tolist()

DA_metabolites = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
print(DA_metabolites)
DA_metabolites_chEBI = name_map[name_map["Query"].isin(DA_metabolites)]['ChEBI'].dropna().tolist()
DA_metabolites_kegg = name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().tolist()
print(len(DA_metabolites_chEBI), "DA ChEBI")
print(len(DA_metabolites_kegg), "DA KEGG")

ORA_reactome = utils.over_representation_analysis(DA_metabolites_chEBI, background_list_chEBI, Reactome_human)
ORA_KEGG = utils.over_representation_analysis(DA_metabolites_kegg, background_list_KEGG, KEGG_pathways)
print(ORA_reactome)
print(len(ORA_reactome[ORA_reactome["P-adjust"] < 0.2]["P-adjust"].tolist()))

print(ORA_KEGG)
print(len(ORA_KEGG[ORA_KEGG["P-adjust"] < 0.2]["P-adjust"].tolist()))
# ORA_reactome.to_csv("Haiwei_ORA_reactome.csv")
ORA_KEGG.to_csv("../Haiwei/Haiwei_ORA_KEGG.csv")
