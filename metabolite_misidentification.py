# Simulation of metabolite misidentification
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import time
import requests
#
# # Import the relevant datasets
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data()
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data()
DEM_brown, background_brown, mat_brown = process_datasets.brown_data()
# # DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM")
# # DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS")
#
# # Import pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
#
datasets = {"Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown]}
            # "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg],
            # "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg]}

# Random misidentification
percentage_misidentifications = [0, 1, 5, 10, 20, 50, 70]
def random_misidentification():
    results_lists = []
    for d in datasets.keys():
        print(d)
        for i in percentage_misidentifications:
            res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][1])
            results_lists.append([d, i] + res)
            time.sleep(1)

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
    res_df.to_csv("Metabolite_misidentification_simulation.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure()
    plt.style.use("seaborn")
    for i in datasets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=14)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=14)
    plt.savefig("metabolite_random_misidentification.png", dpi=300)
    plt.show()

# Misidentification by mass
# Obtained exact masses for KEGG compounds
# Replace with a compound in a similar mass window

KEGG_compounds_masses = pd.read_csv("KEGG_compounds_exact_mass.csv", index_col=0)
# Filter for human compounds only
# print(len(np.setdiff1d(all_KEGG_human_bg, KEGG_compounds_masses.index.tolist())))

def misidentification_mass_plot():
    results_lists = []
    for d in datasets.keys():
        print(d)
        for i in percentage_misidentifications:
            print(i)
            res = utils.misidentify_metabolites_by_mass(i, datasets[d][4], datasets[d][2], KEGG_compounds_masses, datasets[d][1])
            results_lists.append([d, i] + res)

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
    res_df.to_csv("Metabolite_misidentification_by_mass_simulation.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure()
    plt.style.use("seaborn")
    for i in datasets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=14)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=14)
    plt.savefig("metabolite_misidentification_by_mass.png", dpi=300)
    plt.show()

def misidentify_metabolites_by_formula(percentage, processed_matrix, pathway_df):
    #TODO Ensure the replacement compounds are organism-specific
    '''
    Randomly swaps a percentage of KEGG compounds and then performs ORA
    :param percentage: percentage of compounds to be misidentified
    :param processed_matrix: processed abundance matrix with KEGG compounds as columns
    :param pathway_df: list of KEGG pathways
    :return: mean number of p-values significant at P <0.1, Q-values, and standard deviation
    '''

    mat_unannotated = processed_matrix.iloc[:, :-1]
    metabolites = mat_unannotated.columns.tolist()

    # %of compounds without a mass in KEGG
    # in_KEGG = all_compound_masses.index.tolist()
    # print((len(np.setdiff1d(organism_bg, in_KEGG))/len(organism_bg))*100)

    n_misidentified = int(len(metabolites)*(percentage/100))
    p_vals = []
    q_vals = []
    for i in range(0, 1):
        metabolites_to_replace = dict()
        while len(metabolites_to_replace) < n_misidentified:
            replacement_cpd = np.random.choice(metabolites, 1, replace=False)[0]
            cpd_url = 'http://rest.kegg.jp/get/' + replacement_cpd
            data = requests.get(cpd_url)
            cpd_info = data.text
            cpd_info = cpd_info.split("\n")
            for field in cpd_info:
                if field.startswith("FORMULA"):
                    formula = field.split()
                    replacement_compounds_url = 'http://rest.kegg.jp/find/compound/' + formula[1] + '/formula'
                    data = requests.get(replacement_compounds_url)
                    cpd_info = data.text.split("\n")
                    cpd_info = list(filter(None, cpd_info))
                    if len(cpd_info) > 1:
                        potential_replacements = {"compound_to_replace": "", "potential_replacements_ids": []}
                        potential_replacements["compound_to_replace"] = replacement_cpd
                        for x in cpd_info:
                            x = x.split("\t")
                            potential_replacements["potential_replacements_ids"].append(x[0].replace("cpd:", ""))
                        if len(np.setdiff1d(potential_replacements["potential_replacements_ids"], metabolites)) > 0:
                            metabolites_to_replace[replacement_cpd] = potential_replacements
        replacement_dict = {}

        for compound, entries in metabolites_to_replace.items():
            replacement_dict[compound] = np.random.choice(np.setdiff1d(entries["potential_replacements_ids"], metabolites), 1, replace=False)[0]
        print(replacement_dict)

        misidentified_matrix = mat_unannotated.rename(columns=replacement_dict)
        # Perform t-tests and ORA
        ttest_res = utils.t_tests(misidentified_matrix, processed_matrix["Group"], "fdr_bh")
        DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
        ora_res = utils.over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
        p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
        q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    sd_p_signficant_paths = np.std(p_vals)
    sd_q_signficant_paths = np.std(q_vals)
    return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths]

print(misidentify_metabolites_by_formula(50, mat_yamada, KEGG_human_pathways))