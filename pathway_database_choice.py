import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import gridspec
import matplotlib.patches as mpatches
from bioservices import *


#Import the relevant datasets

#Import Reactome datasets
DEM_auwerx_r, background_auwerx_r, mat_auwerx_r = process_datasets.auwerx_data(db="Reactome")
DEM_yamada_r, background_yamada_r, mat_yamada_r = process_datasets.yamada_data(db="Reactome")
DEM_stevens_r, background_stevens_r, mat_stevens_r = process_datasets.stevens_data(db="Reactome")
DEM_brown_r, background_brown_r, mat_brown_r = process_datasets.brown_data(db="Reactome")
DEM_yfgM_r, background_yfgM_r, mat_yfgM_r = process_datasets.zamboni_data("yfgM", db="Reactome")
DEM_dcuS_r, background_dcuS_r, mat_dcuS_r = process_datasets.zamboni_data("dcuS", db="Reactome")

# Import KEGG datasets
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

# Import BioCyc datasets
DEM_auwerx_b, background_auwerx_b, mat_auwerx_b = process_datasets.auwerx_data(db="Cyc")
DEM_yamada_b, background_yamada_b, mat_yamada_b = process_datasets.yamada_data(db="Cyc")
DEM_stevens_b, background_stevens_b, mat_stevens_b = process_datasets.stevens_data(db="Cyc")
DEM_brown_b, background_brown_b, mat_brown_b = process_datasets.brown_data(db="Cyc")
DEM_yfgM_b, background_yfgM_b, mat_yfgM_b = process_datasets.zamboni_data("yfgM", db="Cyc")
DEM_dcuS_b, background_dcuS_b, mat_dcuS_b = process_datasets.zamboni_data("dcuS", db="Cyc")

# Import KEGG pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import Reactome pathway sets
Reactome_pathways = pd.read_csv("Reactome_pathway_set_all_levels.csv", dtype=str, index_col=0)
Reactome_human_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("HSA")]
Reactome_eco_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("ECO")]
Reactome_mouse_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("MMU")]
all_reactome_human_bg = list(set([x for x in Reactome_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_reactome_mouse_bg = list(set([x for x in Reactome_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import BioCyc pathway sets
BioCyc_human_pathways = pd.read_csv("Metacyc_human_pathways.csv")
BioCyc_eco_pathways = pd.read_csv("Metacyc_EColi_pathways.csv")
all_biocyc_human_bg = list(set([x for x in BioCyc_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_biocyc_eco_bg = list(set([x for x in BioCyc_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Labbé": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg],
            "Yachida": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg],
            "Quirós": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg],
            "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg],
            "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg]}

datasets_reactome = {
    "Labbé": [DEM_brown_r, background_brown_r, Reactome_mouse_pathways, all_reactome_mouse_bg],
    "Yachida": [DEM_yamada_r, background_yamada_r, Reactome_human_pathways, all_reactome_human_bg],
    "Stevens": [DEM_stevens_r, background_stevens_r, Reactome_human_pathways, all_reactome_human_bg],
    "Quirós": [DEM_auwerx_r, background_auwerx_r, Reactome_human_pathways, all_reactome_human_bg],
    "Fuhrer (yfgM)": [DEM_yfgM_r, background_yfgM_r, Reactome_human_pathways, all_reactome_human_bg],
    "Fuhrer (dcuS)": [DEM_dcuS_r, background_dcuS_r, Reactome_human_pathways, all_reactome_human_bg]}

datasets_biocyc = {"Labbé": [DEM_brown_b, background_brown_b, BioCyc_human_pathways, all_biocyc_human_bg],
                   "Yachida": [DEM_yamada_b, background_yamada_b, BioCyc_human_pathways, all_biocyc_human_bg],
                   "Stevens": [DEM_stevens_b, background_stevens_b, BioCyc_human_pathways, all_biocyc_human_bg],
                   "Quirós": [DEM_auwerx_b, background_auwerx_b, BioCyc_human_pathways, all_biocyc_human_bg],
                   "Fuhrer (yfgM)": [DEM_yfgM_b, background_yfgM_b, BioCyc_eco_pathways, all_biocyc_eco_bg],
                   "Fuhrer (dcuS)": [DEM_dcuS_b, background_dcuS_b, BioCyc_eco_pathways, all_biocyc_eco_bg]}

print("Data processing complete.")

db_dict = {"KEGG": datasets, "Reactome": datasets_reactome, "BioCyc": datasets_biocyc}
results_dicts = dict.fromkeys(db_dict.keys(), {})
for d in db_dict.keys():
    d_set = db_dict[d]
    res_datasets = {}
    for i in d_set.keys():
        ora_res = utils.over_representation_analysis(d_set[i][0], d_set[i][1], d_set[i][2])
        significant_ids = ora_res[ora_res["P-value"] <= 0.1]["Pathway_ID"]
        res_datasets[i] = significant_ids.tolist()
    results_dicts[d] = res_datasets

kegg_db = KEGG(verbose=False)

# get reactome ChEBI compounds for each pathway
# get reactome ChEBI compounds for each pathway
def process_paths(pathway_df):
    KEGG_pathways = pathway_df.dropna(axis=0, how='all', subset=pathway_df.columns.tolist()[1:])
    pathways = KEGG_pathways.index.tolist()
    pathway_names = KEGG_pathways["Pathway_name"].tolist()
    pathway_dict = {}

    for pathway in pathways:
        pathway_compounds = KEGG_pathways.loc[pathway, :].tolist()
        pathway_compounds = [str(i) for i in pathway_compounds if str(i) != "nan"]
        cpds = pathway_compounds[1:]
        if len(cpds) > 2:
            pathway_dict[pathway] = cpds
    return pathway_dict

Reactome_pathways = pd.read_csv("Reactome_pathway_set_all_levels.csv", dtype=str, index_col=0)
Reactome_pathways = process_paths(Reactome_pathways)

ch = ChEBI()
reactome_sig_path_cpds = {}

for name, v in results_dicts["Reactome"].items():
    pathways = v
    all_compounds_in_pathways = []
    for pathway in pathways:
        pathway_compounds = Reactome_pathways[pathway]
        pathway_compounds_kegg = []
        for cpd in pathway_compounds:
            try:
                kegg_id = ch.conv("CHEBI:" + cpd, "KEGG COMPOUND accession")
                pathway_compounds_kegg.append(kegg_id[0])
            except (ValueError, AttributeError):
                continue
        print(pathway, len(pathway_compounds), len(pathway_compounds_kegg))
        all_compounds_in_pathways = all_compounds_in_pathways + pathway_compounds_kegg
    reactome_sig_path_cpds[name] = list(set(all_compounds_in_pathways))

KEGG_human_pathways = process_paths(pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0))
KEGG_eco_pathways = process_paths(pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0))
KEGG_mouse_pathways = process_paths(pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0))
kegg_pathways_merged = {**KEGG_human_pathways, **KEGG_eco_pathways}
all_kegg_pathways_merged = {**kegg_pathways_merged, **KEGG_mouse_pathways}

kegg_sig_path_cpds = {}

for name, v in results_dicts["KEGG"].items():
    pathways = v
    all_compounds_in_pathways = []
    for pathway in pathways:
        pathway_compounds = all_kegg_pathways_merged[pathway]
        print(pathway, len(pathway_compounds))
        all_compounds_in_pathways = all_compounds_in_pathways + pathway_compounds
    kegg_sig_path_cpds[name] = list(set(all_compounds_in_pathways))

def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

res = {}
for dset in kegg_sig_path_cpds.keys():
    JI = jaccard_similarity(kegg_sig_path_cpds[dset], reactome_sig_path_cpds[dset])
    res[dset] = JI



# res_lists = []
# for d in [datasets, datasets_reactome, datasets_biocyc]:
#     num = 1
#     for i in d.keys():
#         ora_res = utils.over_representation_analysis(d[i][0], d[i][1], d[i][2])
#         n_significant_p_01 = len(ora_res[ora_res["P-value"] <= 0.1]["P-value"])
#         print(n_significant_p_01)
#         print([num, i, n_significant_p_01])
#         res_lists.append([num, i, n_significant_p_01])
#     num += 1
# df = pd.DataFrame(res_lists)
# df.to_csv("all_databases_significant_paths.csv")
#
# #
# # # def identify_significant_paths(db="KEGG"):
# # #     d_sets = datasets
# # #     if db == "Reactome":
# # #         d_sets = datasets_reactome
# # #     if db == "Cyc":
# # #         d_sets = datasets_biocyc
# # #     res_list = []
# # #     for i in d_sets.keys():
# # #         ora_res = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
# # #         n_significant_p_01 = len(ora_res[ora_res["P-value"] <= 0.1])
# # #         res = [i, n_significant_p_01]
# # #         res_list.append(res)
# # #     df = pd.DataFrame(res_list)
# # #     df.columns = ["Dataset", db]
# # #     print(df)
# # #     df.to_csv("significant_pathways_" + db + ".csv")
# # #
# # #
# # # identify_significant_paths(db="Cyc")
# #
# # bio_paths = pd.read_csv("significant_pathways_Cyc.csv", index_col=1)
# # bio_paths.rename(columns={'Cyc':'BioCyc'}, inplace=True)
# # kegg_paths = pd.read_csv("significant_pathways_KEGG.csv", index_col=1)
# # reactome_paths = pd.read_csv("significant_pathways_Reactome.csv", index_col=1)
# # dfs = pd.DataFrame([kegg_paths.iloc[:, 1], reactome_paths.iloc[:, 1], bio_paths.iloc[:,1]])
# # print(dfs)
# quit()
# yamada_ipa = pd.read_csv("yamada_IPA_results.csv")
# yamada_ipa["pvalue"] = np.exp(-yamada_ipa[" -log(p-value)"])
# print(yamada_ipa['pvalue'].tolist())
# print(yamada_ipa[yamada_ipa["pvalue"] <= 0.1])
#
# brown_ipa = pd.read_csv("brown_IPA_results.csv")
# brown_ipa["pvalue"] = np.exp(-brown_ipa[" -log(p-value)"])
# print(brown_ipa['pvalue'].tolist())
# print(brown_ipa[brown_ipa["pvalue"] <= 0.1])
#
# stevens_ipa = pd.read_csv("stevens_IPA_results.csv")
# stevens_ipa["pvalue"] = np.exp(-stevens_ipa[" -log(p-value)"])
# print(stevens_ipa['pvalue'].tolist())
# print(stevens_ipa[stevens_ipa["pvalue"] <= 0.1])
#
# quiros_ipa = pd.read_csv("quiros_IPA_results.csv")
# quiros_ipa["pvalue"] = np.exp(-quiros_ipa[" -log(p-value)"])
# print(quiros_ipa['pvalue'].tolist())
# print(quiros_ipa[quiros_ipa["pvalue"] <= 0.1])
#
# dcus_ipa = pd.read_csv("zamboni_dcuS_IPA_results.csv")
# dcus_ipa["pvalue"] = np.exp(-dcus_ipa[" -log(p-value)"])
# print(dcus_ipa['pvalue'].tolist())
# print(dcus_ipa[dcus_ipa["pvalue"] <= 0.1])
#
# yfgm_ipa = pd.read_csv("zamboni_yfgM_IPA_results.csv")
# yfgm_ipa["pvalue"] = np.exp(-yfgm_ipa[" -log(p-value)"])
# print(yfgm_ipa['pvalue'].tolist())
# print(yfgm_ipa[yfgm_ipa["pvalue"] <= 0.1])
#
# # plt.style.use("seaborn")
# #
# # ax = dfs.plot.barh(rot=0, figsize=(8,4), width=0.8)
# # plt.xlabel("Number of significant pathways at P ≤ 0.1", fontsize=13)
# # plt.ylabel("Pathway database (organism-specific)", fontsize=13)
# # plt.tight_layout()
# # # plt.savefig("../Figures/pathway_significant_comparison.png", dpi=600)
# #
# # plt.show()
#
