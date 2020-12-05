# Simulation of metabolite misidentification
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Import the relevant datasets
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data()
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data()
DEM_brown, background_brown, mat_brown = process_datasets.brown_data()
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS")

# Import pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown]}
            # "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg],
            # "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg]}

# Random misidentification
percentage_misidentifications = [0, 1, 5, 10, 20, 50, 70]

results_lists = []
for d in datasets.keys():
    print(d)
    for i in percentage_misidentifications:
        res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2])
        results_lists.append([d, i] + res)

res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage reduction", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
res_df.to_csv("Metabolite_misidentification_simulation.csv")

simulation_res = res_df
print(simulation_res.head)
plt.figure()
plt.style.use("seaborn")
for i in datasets.keys():
    plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                 simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                 yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                 label=i, fmt='o', linestyle="solid", capsize=20, markersize=5)
plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
plt.legend()
plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=14)
plt.xlabel("Percentage of metabolites misidentified", fontsize=14)
plt.savefig("metabolite_random_misidentification.png", dpi=300)
plt.show()

# def misidentify_metabolites(percentage, processed_matrix, organism_compounds, background_list, pathway_df):
#     '''
#     Randomly swaps a percentage of KEGG compounds and then performs ORA
#     :param percentage: percentage of compounds to be misidentified
#     :param processed_matrix: processed abundance matrix with KEGG compounds as columns
#     :param organism_compounds: all KEGG compounds for the organism
#     :param background_list: background list for dataset
#     :param pathway_df: list of KEGG pathways
#     :return: mean number of p-values significant at P <0.1, Q-values, and standard deviation
#     '''
#
#     mat_unannotated = processed_matrix.iloc[:, :-1]
#     metabolites = mat_unannotated.columns.tolist()
#     n_misidentified = int(len(metabolites)*(percentage/100))
#     p_vals = []
#     q_vals = []
#     for i in range(0, 10):
#         # Randomly replace n compounds
#         metabolites_to_replace = np.random.choice(metabolites, n_misidentified, replace=False)
#         replacement_compounds = np.random.choice(np.setdiff1d(organism_compounds, background_list), n_misidentified, replace=False)
#         replacement_dict = dict(zip(metabolites_to_replace, replacement_compounds))
#         misidentified_matrix = mat_unannotated.rename(columns=replacement_dict)
#
#         # Perform t-tests and ORA
#         ttest_res = utils.t_tests(misidentified_matrix, processed_matrix["Group"], "fdr_bh")
#         DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
#         ora_res = utils.over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
#         p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
#         q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
#     mean_p_signficant_paths = np.mean(p_vals)
#     mean_q_signficant_paths = np.mean(q_vals)
#     sd_p_signficant_paths = np.std(p_vals)
#     sd_q_signficant_paths = np.std(q_vals)
#     return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths]


