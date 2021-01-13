# Simulation of metabolite misidentification
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import time

# Import the relevant datasets
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="Reactome")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="Reactome")
# DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="Reactome")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="Reactome")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="Reactome")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="Reactome")

# DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
# DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
# # DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
# DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

# Import pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import Reactome pathway sets
Reactome_pathways = pd.read_csv("Reactome_pathway_set.csv", dtype=str, index_col=0)
Reactome_human_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("HSA")]
Reactome_eco_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("ECO")]
Reactome_mouse_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("MMU")]
all_reactome_human_bg = list(set([x for x in Reactome_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_reactome_mouse_bg = list(set([x for x in Reactome_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import compounds and masses
KEGG_compounds_masses = pd.read_csv("KEGG_compounds_masses_estimated.csv", names=["compound", "formula", "mass"])
Reactome_compounds_masses = pd.read_csv("CHEBI_compounds_masses.csv", names=["compound", "formula", "mass"])

datasets = {"Auwerx": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx, [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)], KEGG_compounds_masses],
            "Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)], KEGG_compounds_masses],
            # "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)], KEGG_compounds_masses],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], KEGG_compounds_masses],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], KEGG_compounds_masses]}

datasets_reactome = {"Auwerx": [DEM_auwerx, background_auwerx, Reactome_human_pathways, all_reactome_human_bg, mat_auwerx, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], Reactome_compounds_masses],
                     "Yamada": [DEM_yamada, background_yamada, Reactome_human_pathways, all_reactome_human_bg, mat_yamada, [i for i in range(0, 12, 1)], [i for i in range(0, 9, 1)], Reactome_compounds_masses],
                     # "Stevens": [DEM_stevens, background_stevens, Reactome_human_pathways, all_reactome_human_bg],
                     "Brown": [DEM_brown, background_brown, Reactome_mouse_pathways, all_reactome_mouse_bg, mat_brown, [i for i in range(0, 12, 1)], [i for i in range(0, 10, 1)], Reactome_compounds_masses],
                     "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, Reactome_human_pathways, all_reactome_human_bg, mat_yfgM, [i for i in range(0, 6, 1)], [i for i in range(0, 5, 1)], Reactome_compounds_masses],
                     "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, Reactome_human_pathways, all_reactome_human_bg, mat_dcuS, [i for i in range(0, 6, 1)], [i for i in range(0, 5, 1)], Reactome_compounds_masses]}

print("Processing complete.")
# Random misidentification
percentage_misidentifications = [i for i in range(0, 100, 10)]
def random_misidentification(db="KEGG"):
    results_lists = []
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    for d in d_sets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in [i for i in range(0, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, d_sets[d][4], d_sets[d][3], d_sets[d][1], d_sets[d][2],
                                                    zamboni=True)
                results_lists.append([d, i] + res[:-1])
        else:
            for i in percentage_misidentifications:
                print(i)
                res = utils.misidentify_metabolites(i, d_sets[d][4], d_sets[d][3], d_sets[d][1], d_sets[d][2])
                results_lists.append([d, i] + res[:-1])

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
    # res_df.to_csv("Metabolite_misidentification_simulation_random_Reactome.csv")

    simulation_res = res_df
    plt.figure()
    plt.style.use("seaborn")
    for i in d_sets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    plt.title("Reactome", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations")
    plt.xlabel("Percentage of metabolites misidentified")
    plt.savefig("metabolite_random_misidentification_Reactome2.png", dpi=300)
    plt.show()

# random_misidentification(db="Reactome")

# Parameter grid for TPR/FPR heatmaps. Numbers correspond to indexes in datasets param grids.
param_grid_heatmaps = {"random": [utils.misidentify_metabolites, 4, 3, 1, 2, [i for i in range(10, 70, 10)]],
                       "mass": [utils.misidentify_metabolites_by_mass, 4, 2, 7, 3, [i for i in range(1, 6, 1)]],
                       "formula": [utils.misidentify_metabolites_by_formula, 4, 2, 7, 3, [i for i in range(1, 5, 1)]]}

def TPR_heatmap(pg, fname, db="KEGG"):
    """
    Plots TRP/FPR heatmap
    :param misidenetification_funct: function for misidentification
    :param db: Database, default is KEGG
    :return: heatmap
    """
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    results_TPR = []
    results_FPR = []
    for d in ["Auwerx", "Brown", "Zamboni (dcuS)", "Zamboni (yfgM)"]:
        print(d)
        if d.startswith("Zamboni"):
            original_pathways = pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                                                    zamboni=True)[4][0]

            for i in pg[5]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                                                    zamboni=True)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways) # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths/total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways])/total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])

        else:
            original_pathways = pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4][0]
            for i in pg[5]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways) # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths/total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways])/total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])
    res_df_TPR = pd.DataFrame(results_TPR,
                          columns=["Dataset", "Percentage misidentification", "Average fraction"])
    res_df_FPR = pd.DataFrame(results_FPR,
                              columns=["Dataset", "Percentage misidentification", "Average fraction"])

    res_df_TPR = res_df_TPR.pivot(index='Percentage misidentification', columns='Dataset', values='Average fraction')
    res_df_FPR = res_df_FPR.pivot(index='Percentage misidentification', columns='Dataset', values='Average fraction')
    # res_df = pd.read_csv("metabolite_misidentification_heatmap.csv", index_col=0)
    plt.style.use("seaborn")
    plt.figure(figsize=(9, 5))
    plt.subplot(121)
    plt.title('TPR')
    sns.heatmap(res_df_TPR, annot=True, cmap="mako", square=True)
    plt.subplot(122)
    plt.title('FPR')
    sns.heatmap(res_df_FPR, annot=True, cmap="rocket_r", square=True)
    plt.subplots_adjust(bottom=0.28)
        # plt.ylabel("Percentage metabolite misidentification (%)")
    plt.savefig(fname, dpi=300)
    plt.show()

# TPR_heatmap(param_grid_heatmaps["random"], "formula_misidentification_heatmap_Reactome_new.png", db="Reactome")
TPR_heatmap(param_grid_heatmaps["formula"], "formula_misidentification_heatmap_Reactome_new.png", db="Reactome")
TPR_heatmap(param_grid_heatmaps["mass"], "mass_misidentification_heatmap_Reactome_new.png", db="Reactome")



def misidentification_mass_plot(db="KEGG"):
    d_sets = datasets
    masses = KEGG_compounds_masses
    if db == "Reactome":
        d_sets = datasets_reactome
        masses = Reactome_compounds_masses
    results_lists = []
    for d in d_sets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in d_sets[d][5]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, d_sets[d][4], d_sets[d][2], masses,
                                                d_sets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in d_sets[d][5]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, d_sets[d][4], d_sets[d][2], masses,
                                                            d_sets[d][3], zamboni=False)
                results_lists.append([d, i] + res)

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std", "sig_paths"])
    # res_df.to_csv("Metabolite_misidentification_by_mass_simulation_Reactome.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure()
    plt.style.use("seaborn")
    for i in d_sets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    # plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1")
    plt.xlabel("Percentage of metabolites misidentified by mass")
    plt.title("KEGG", fontsize=14)
    # plt.savefig("metabolite_misidentification_by_mass_Reactome.png", dpi=300)
    plt.show()

# misidentification_mass_plot(db="KEGG")

def misidentification_formula_plot(db="KEGG"):
    d_sets = datasets
    compounds = KEGG_compounds_masses
    if db == "Reactome":
        d_sets = datasets_reactome
        compounds = Reactome_compounds_masses
    results_lists = []
    for d in d_sets.keys():
        print(d)
        if d.startswith("Zamboni"):
            for i in d_sets[d][6]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, d_sets[d][4], d_sets[d][2], compounds,
                                                d_sets[d][3], zamboni=True)
                results_lists.append([d, i] + res)
        else:
            for i in d_sets[d][6]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, d_sets[d][4], d_sets[d][2], compounds,
                                                            d_sets[d][3], zamboni=False)
                results_lists.append([d, i] + res)

    res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std", "sig_paths"])
    res_df.to_csv("Metabolite_misidentification_by_formula_simulation_Reactome.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure()
    plt.style.use("seaborn")
    for i in d_sets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    plt.title("Reactome", fontsize=14)
    plt.legend()
    plt.ylabel("Mean number of pathways significant at P < 0.1")
    plt.xlabel("Percentage of metabolites misidentified by formula")
    plt.savefig("metabolite_misidentification_by_formula_Reactome.png", dpi=300)
    plt.show()

# misidentification_formula_plot(db="Reactome")