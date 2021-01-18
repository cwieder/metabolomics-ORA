import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import utils

mat = pd.read_csv("miyamoto_data.csv", index_col=0).T

md_raw = pd.read_csv("Miyamoto_metadata.csv")

cancer_status = dict(zip(md_raw['Sample name'], md_raw['Cancer status']))

metadata_list = list(zip(md_raw['Organ'], md_raw['Cancer status'], md_raw['Smoker'], md_raw['Gender']))
metadata_dict = dict(zip(md_raw['Sample name'].values, metadata_list))

organ_sample = dict(zip(md_raw['Sample name'], md_raw['Organ']))
mat["Group"] = mat.index.map(organ_sample)
mat_serum = mat[mat["Group"] == "Serum"]
mat_copy = mat_serum.copy()

mat_copy["Group"] = mat_copy.index.map(cancer_status)

matrix_proc = utils.data_processing(mat_serum.iloc[:,:-1], 0, 0)

utils.plot_PCA(matrix_proc, mat_copy["Group"], "Miyamoto (cancer - serum)", 86)
ttest_res = utils.t_tests(matrix_proc, mat_copy["Group"], "fdr_bh")
DA_metabolites = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
print(len(DA_metabolites))

name_map = pd.read_csv("../Miyamoto/name_map.csv", dtype=str)
KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)

Reactome_pathways = pd.read_csv("Reactome_pathway_set.csv", dtype=str, index_col=0)
Reactome_human_ids = [i for i in Reactome_pathways.index if i.startswith("R-HSA")]
Reactome_human = Reactome_pathways[Reactome_pathways.index.isin(Reactome_human_ids)]

background_list = matrix_proc.columns.tolist()

background_list_chEBI = name_map[name_map["Query"].isin(background_list)]['ChEBI'].dropna().tolist()
background_list_KEGG = name_map[name_map["Query"].isin(background_list)]['KEGG'].dropna().tolist()

# DA_metabolites = ttest_res[ttest_res["P-adjust"] < 0.1]["Metabolite"].tolist()

name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().to_csv("../Miyamoto/Miyamoto_DA.csv")

print(DA_metabolites)
DA_metabolites_chEBI = name_map[name_map["Query"].isin(DA_metabolites)]['ChEBI'].dropna().tolist()
DA_metabolites_kegg = name_map[name_map["Query"].isin(DA_metabolites)]['KEGG'].dropna().tolist()
print(len(DA_metabolites_chEBI), "DA ChEBI")
print(len(DA_metabolites_kegg), "DA KEGG")

ORA_reactome = utils.over_representation_analysis(DA_metabolites_chEBI, background_list_chEBI, Reactome_human)
ORA_KEGG = utils.over_representation_analysis(DA_metabolites_kegg, background_list_KEGG, KEGG_pathways)
# print(ORA_reactome)
# print(len(ORA_reactome[ORA_reactome["P-value"] < 0.1]["P-adjust"].tolist()))

print(ORA_KEGG)
print(len(ORA_KEGG[ORA_KEGG["P-value"] < 0.1]["P-adjust"].tolist()))
# ORA_reactome.to_csv("../Miyamoto/Miyamoto_ORA_reactome.csv")
# ORA_KEGG.to_csv("../Miyamoto/Miyamoto_ORA_KEGG.csv")
