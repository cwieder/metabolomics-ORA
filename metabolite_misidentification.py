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
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data()
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

datasets = {"Auwerx": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx, [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)]],
            "Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            # "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS]}

# Random misidentification
percentage_misidentifications = [i for i in range(0, 100, 10)]
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
        else:
            for i in percentage_misidentifications:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2])
                results_lists.append([d, i] + res[:-1])

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
    # plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations")
    plt.xlabel("Percentage of metabolites misidentified")
    plt.savefig("metabolite_random_misidentification.png", dpi=300)
    plt.show()


def TPR_heatmap(TPR=False, FPR=False):
    """
    Plots precision/TRP/FPR heatmap
    :param TPR: plot TRP
    :param FPR: plot FPR
    :return: heatmap
    """
    results_TPR = []
    results_FPR = []
    for d in ["Auwerx", "Yamada", "Brown", "Zamboni (dcuS)", "Zamboni (yfgM)"]:
        print(d)
        if d.startswith("Zamboni"):
            original_pathways = utils.misidentify_metabolites(0, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                                    zamboni=True)[4][0]

            for i in [i for i in range(10, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                                    zamboni=True)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways) # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths/total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = (total_significant_paths-number_common_paths)/total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])

        else:
            original_pathways = \
            utils.misidentify_metabolites(0, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2],
                                          zamboni=False)[4][0]
            for i in [i for i in range(10, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, datasets[d][4], datasets[d][3], datasets[d][1], datasets[d][2])[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways) # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths/total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = (total_significant_paths-number_common_paths)/total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])
        time.sleep(1)
    res_df_TPR = pd.DataFrame(results_TPR,
                          columns=["Dataset", "Percentage misidentification", "Average fraction"])
    res_df_FPR = pd.DataFrame(results_FPR,
                              columns=["Dataset", "Percentage misidentification", "Average fraction"])

    res_df_TPR = res_df_TPR.pivot(index='Percentage misidentification', columns='Dataset', values='Average fraction')
    res_df_FPR = res_df_FPR.pivot(index='Percentage misidentification', columns='Dataset', values='Average fraction')
    # res_df.to_csv("metabolite_misidentification_heatmap.csv")
    # res_df = pd.read_csv("metabolite_misidentification_heatmap.csv", index_col=0)
    plt.style.use("seaborn")
    plt.figure(figsize=(9, 5))
    plt.subplot(121)
    plt.title('TPR')
    sns.heatmap(res_df_TPR, annot=True, cmap="mako", square=True)
    plt.subplot(122)
    plt.title('FPR')
    sns.heatmap(res_df_FPR, annot=True, cmap="rocket", square=True)
    plt.subplots_adjust(bottom=0.28)
        # plt.ylabel("Percentage metabolite misidentification (%)")
    plt.savefig("random_misidentification_heatmap.png", dpi=300)
    plt.show()

TPR_heatmap(TPR=True)

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
            for i in [i for i in range(0, 7, 1)]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, datasets[d][4], datasets[d][2], KEGG_compounds_masses,
                                                datasets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in datasets[d][5]:
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
    # plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1")
    plt.xlabel("Percentage of metabolites misidentified")
    plt.savefig("metabolite_misidentification_by_mass.png", dpi=300)
    plt.show()

# misidentification_mass_plot()

KEGG_compounds_formula = pd.read_csv("KEGG_compound_formulae.csv", index_col=0)

def misidentification_formula_plot():
    results_lists = []
    for d in datasets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in [i for i in range(0, 6, 1)]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, datasets[d][4], datasets[d][2], KEGG_compounds_formula,
                                                datasets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in datasets[d][6]:
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
    plt.ylabel("Mean number of pathways significant at P < 0.1")
    plt.xlabel("Percentage of metabolites misidentified")
    plt.savefig("metabolite_misidentification_by_formula.png", dpi=300)
    plt.show()

# misidentification_formula_plot()