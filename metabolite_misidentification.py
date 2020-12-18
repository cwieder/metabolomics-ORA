# Simulation of metabolite misidentification
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import time
import requests

# Import the relevant datasets
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data()
# DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data()
DEM_brown, background_brown, mat_brown = process_datasets.brown_data()
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS")

# Import pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada],
            # "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS]}

# Random misidentification
percentage_misidentifications = [i for i in range(0, 80, 10)]
def random_misidentification():
    results_lists = []
    for d in datasets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in [i for i in range(0, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                                    zamboni=True)
                results_lists.append([d, i] + res[:-1])
                time.sleep(1)
        else:
            for i in percentage_misidentifications:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2])
                results_lists.append([d, i] + res[:-1])
                time.sleep(1)
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

def TPR_heatmap(TPR=False, FPR=False):
    results_lists = []
    for d in ["Yamada", "Brown", "Zamboni (dcuS)", "Zamboni (yfgM)"]:
        print(d)
        if d.startswith("Zamboni"):
            original_pathways = utils.misidentify_metabolites(0, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                                    zamboni=True)[4][0]

            for i in [i for i in range(10, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                                    zamboni=True)[4]
                misidentified_pathways = res
                pathway_fractions = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways)
                    number_common_paths = len([i for i in x if i in original_pathways])
                    if TPR:
                        fraction_pathways = number_common_paths/total_significant_paths
                        pathway_fractions.append(fraction_pathways)
                    elif FPR:
                        fraction_pathways = (total_significant_paths-number_common_paths)/total_significant_paths
                        pathway_fractions.append(fraction_pathways)
                avg_fraction = np.mean(pathway_fractions)
                results_lists.append([d, i, avg_fraction])

        else:
            original_pathways = \
            utils.misidentify_metabolites(0, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                          zamboni=False)[4][0]
            misidentified_pathways = []
            for i in [i for i in range(10, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2])[4]
                misidentified_pathways = res
                pathway_fractions = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways)
                    number_common_paths = len([i for i in x if i in original_pathways])
                    if TPR:
                        fraction_pathways = number_common_paths / total_significant_paths
                        pathway_fractions.append(fraction_pathways)
                    elif FPR:
                        fraction_pathways = (total_significant_paths - number_common_paths) / total_significant_paths
                        pathway_fractions.append(fraction_pathways)
                    time.sleep(1)
                avg_fraction = np.mean(pathway_fractions)
                results_lists.append([d, i, avg_fraction])
        time.sleep(1)
    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage misidentification", "Average fraction"])
    res_df = res_df.pivot(index='Percentage misidentification', columns='Dataset', values='Average fraction')
    res_df.to_csv("metabolite_misidentification_heatmap.csv")
    res_df = pd.read_csv("metabolite_misidentification_heatmap.csv", index_col=0)
    ax = sns.heatmap(res_df, annot=True, cmap="rocket")
    plt.subplots_adjust(bottom=0.28)
    plt.ylabel("Percentage metabolite misidentification (%)")
    plt.savefig("random_misidentification_heatmap_FPR.png", dpi=300)
    plt.show()

TPR_heatmap(FPR=True)


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
        if d.startswith("Zamboni"):
            for i in [i for i in range(0, 8, 2)]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, datasets[d][4], datasets[d][2], KEGG_compounds_masses,
                                                datasets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in [i for i in range(0, 12, 2)]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, datasets[d][4], datasets[d][2], KEGG_compounds_masses,
                                                            datasets[d][3], zamboni=False)
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
    plt.ylabel("Mean number of pathways significant at P < 0.1", fontsize=14)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=14)
    plt.savefig("metabolite_misidentification_by_mass.png", dpi=300)
    plt.show()

# misidentification_mass_plot()

KEGG_compounds_formula = pd.read_csv("KEGG_compound_formulae.csv", index_col=0)

def misidentification_formula_plot():
    results_lists = []
    for d in datasets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in [i for i in range(0, 6, 2)]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, datasets[d][4], datasets[d][2], KEGG_compounds_formula,
                                                datasets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in [i for i in range(0, 12, 2)]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, datasets[d][4], datasets[d][2], KEGG_compounds_formula,
                                                            datasets[d][3], zamboni=False)
                results_lists.append([d, i] + res)

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
    res_df.to_csv("Metabolite_misidentification_by_formula_simulation.csv")

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
    plt.ylabel("Mean number of pathways significant at P < 0.1", fontsize=14)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=14)
    plt.savefig("metabolite_misidentification_by_formula.png", dpi=300)
    plt.show()

# misidentification_formula_plot()