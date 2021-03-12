# Simulation of metabolite misidentification
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

# Import the relevant datasets
# DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="Reactome")
# DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="Reactome")
# # DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="Reactome")
# DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="Reactome")
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="Reactome")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="Reactome")

DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

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

datasets = {"Labbé": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown,
                      [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)], KEGG_compounds_masses],
            "Yachida": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada,
                        [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)], KEGG_compounds_masses],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens,
                        [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)], KEGG_compounds_masses],
            "Quirós": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx,
                       [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)], KEGG_compounds_masses],
            "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM,
                              [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], KEGG_compounds_masses],
            "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS,
                              [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], KEGG_compounds_masses]}

datasets_reactome = {
    "Quirós": [DEM_auwerx, background_auwerx, Reactome_human_pathways, all_reactome_human_bg, mat_auwerx,
               [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)], Reactome_compounds_masses],
    "Yachida": [DEM_yamada, background_yamada, Reactome_human_pathways, all_reactome_human_bg, mat_yamada,
                [i for i in range(0, 12, 1)], [i for i in range(0, 9, 1)], Reactome_compounds_masses],
    # "Stevens": [DEM_stevens, background_stevens, Reactome_human_pathways, all_reactome_human_bg],
    "Labbé": [DEM_brown, background_brown, Reactome_mouse_pathways, all_reactome_mouse_bg, mat_brown,
              [i for i in range(0, 12, 1)], [i for i in range(0, 10, 1)], Reactome_compounds_masses],
    "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, Reactome_human_pathways, all_reactome_human_bg, mat_yfgM,
                      [i for i in range(0, 6, 1)], [i for i in range(0, 5, 1)], Reactome_compounds_masses],
    "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, Reactome_human_pathways, all_reactome_human_bg, mat_dcuS,
                      [i for i in range(0, 6, 1)], [i for i in range(0, 5, 1)], Reactome_compounds_masses]}

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
        if d.startswith("Fuhrer"):
            for i in [i for i in range(0, 70, 10)]:
                print(i)
                res = utils.misidentify_metabolites(i, d_sets[d][4], d_sets[d][3], d_sets[d][1], d_sets[d][2],
                                                    zamboni=True)
                results_lists.append([d, i] + res[:-2])
        else:
            for i in percentage_misidentifications:
                print(i)
                res = utils.misidentify_metabolites(i, d_sets[d][4], d_sets[d][3], d_sets[d][1], d_sets[d][2])
                results_lists.append([d, i] + res[:-2])

    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std",
                                   "q_std"])
    # res_df.to_csv("Metabolite_misidentification_simulation_random_Reactome.csv")

    simulation_res = res_df
    plt.figure(figsize=(5, 4), dpi=600)
    plt.style.use("seaborn")
    for i in d_sets.keys():
        if i in ["Quirós", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                         simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                         yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                         label=i, fmt='o', linestyle="--", capsize=5, markeredgewidth=2, markersize=4)

        else:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                         simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                         yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                         label=i, fmt='o', linestyle="solid", capsize=5, markeredgewidth=2, markersize=4)

    plt.title("Random misidentification", fontsize=14)
    plt.legend(fontsize=11)
    plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=13)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=13)
    plt.savefig("metabolite_random_misidentification_KEGG.png", dpi=600)
    plt.show()


# random_misidentification(db="KEGG")

# Parameter grid for TPR/FPR heatmaps. Numbers correspond to indexes in datasets param grids.
param_grid_heatmaps = {"random": [utils.misidentify_metabolites, 4, 3, 1, 2, [i for i in range(10, 70, 10)]],
                       "mass": [utils.misidentify_metabolites_by_mass, 4, 2, 7, 3, [i for i in range(1, 7, 1)]],
                       "formula": [utils.misidentify_metabolites_by_formula, 4, 2, 7, 3, [i for i in range(1, 6, 1)]]}


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
    for d in d_sets.keys():
        print(d)
        if d.startswith("Fuhrer"):
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
                    total_significant_paths = len(original_pathways)  # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths / total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways]) / total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])

        else:
            original_pathways = \
                pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4][0]
            for i in pg[5]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways)  # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths / total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways]) / total_significant_paths
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
    print(res_df_FPR)
    res_df_TPR = res_df_TPR.reindex(columns=list(d_sets.keys()))
    res_df_FPR = res_df_FPR.reindex(columns=list(d_sets.keys()))
    # res_df = pd.read_csv("metabolite_misidentification_heatmap.csv", index_col=0)
    plt.style.use("seaborn")
    plt.figure(figsize=(9, 5), dpi=600)
    plt.subplot(121)
    plt.title('Pathway re-discovery rate')
    sns.heatmap(res_df_TPR, annot=True, cmap="mako", square=True)
    plt.subplot(122)
    plt.title('Pathway mis-discovery rate')
    sns.heatmap(res_df_FPR, annot=True, cmap="rocket_r", square=True)
    plt.subplots_adjust(bottom=0.28)
    # plt.ylabel("Percentage metabolite misidentification (%)")
    plt.tight_layout()
    # plt.savefig(fname, dpi=600)
    plt.show()

# TPR_heatmap(param_grid_heatmaps["random"], "random_misidentification_heatmap_KEGG_new.png", db="KEGG")
# TPR_heatmap(param_grid_heatmaps["formula"], "formula_misidentification_heatmap_KEGG_new.png", db="KEGG")
# TPR_heatmap(param_grid_heatmaps["mass"], "mass_misidentification_heatmap_KEGG_new.png", db="KEGG")

def plot_ROC(pg, fname, db="KEGG"):
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
    # for d in ["Quirós", "Labbé", "Yachida", "Fuhrer (dcuS)", "Fuhrer (yfgM)"]:
    for d in d_sets.keys():
        print(d)
        if d.startswith("Fuhrer"):
            original_pathways = pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                                      zamboni=True)[4][0]

            for i in pg[5]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                            zamboni=True)[4]
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in res:
                    total_significant_paths = len(original_pathways)  # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths / total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways]) / total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR])
                results_FPR.append([d, i, avg_fraction_FPR])

        else:
            original_pathways = \
                pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4][0]
            non_signif_original_pathways = \
                pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[5][0]
            total_significant_paths = len(original_pathways)  # True positive
            total_non_significant_paths = len(non_signif_original_pathways)  # True negative
            for i in pg[5]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = number_common_paths / total_significant_paths
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len([i for i in x if i not in original_pathways]) / total_significant_paths
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
    print(res_df_FPR)
    print(res_df_TPR)
    plt.style.use("seaborn")
    plt.figure(figsize=(6, 6))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([1, 0], [1, 0], color='black', linestyle=':', label="y=x")
    for i in d_sets.keys():
        xs = res_df_FPR[i]
        ys = res_df_TPR[i]
        plt.plot(xs, ys, label=i, marker="o")

    plt.title('ROC curve (mass misidentification)')
    plt.legend()
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    # plt.savefig(fname, dpi=300)
    plt.show()

    def auc_from_fpr_tpr(fpr, tpr, trapezoid=False):
        inds = [i for (i, (s, e)) in enumerate(zip(fpr[: -1], fpr[1:])) if s != e] + [len(fpr) - 1]
        fpr, tpr = fpr[inds], tpr[inds]
        area = 0
        ft = zip(fpr, tpr)
        for p0, p1 in zip(ft[: -1], ft[1:]):
            area += (p1[0] - p0[0]) * ((p1[1] + p0[1]) / 2 if trapezoid else p0[1])
        return area



def misidentification_mass_plot(db="KEGG"):
    d_sets = datasets
    masses = KEGG_compounds_masses
    if db == "Reactome":
        d_sets = datasets_reactome
        masses = Reactome_compounds_masses
    results_lists = []
    for d in d_sets.keys():
        print(d)
        if d.startswith("Fuhrer"):
            for i in d_sets[d][5]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, d_sets[d][4], d_sets[d][2], masses,
                                                            d_sets[d][3], zamboni=True)
                results_lists.append([d, i] + res[:-1])
        else:
            for i in d_sets[d][5]:
                print(i)
                res = utils.misidentify_metabolites_by_mass(i, d_sets[d][4], d_sets[d][2], masses,
                                                            d_sets[d][3], zamboni=False)
                results_lists.append([d, i] + res[:-1])

    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std",
                                   "q_std", "sig_paths"])
    # res_df.to_csv("Metabolite_misidentification_by_mass_simulation_Reactome.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure(figsize=(5, 4), dpi=600)
    plt.style.use("seaborn")
    for i in d_sets.keys():
        if i in ["Quirós", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                         simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                         yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                         label=i, fmt='o', linestyle="--", capsize=5, markeredgewidth=2, markersize=4)
        else:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                         simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                         yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                         label=i, fmt='o', linestyle="solid", capsize=5, markeredgewidth=2, markersize=4)
    # plt.title("Number of pathways with P-values < 0.1 in response to \n varying levels of metabolite misidentification", fontsize=14)
    plt.legend(fontsize=11)
    plt.ylabel("Mean number of pathways significant at P < 0.1", fontsize=13)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=13)
    plt.title("Misidentification by mass", fontsize=14)
    # plt.savefig("metabolite_misidentification_by_mass_KEGG.png", dpi=600)
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
        if d.startswith("Fuhrer"):
            for i in d_sets[d][6]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, d_sets[d][4], d_sets[d][2], compounds,
                                                               d_sets[d][3], zamboni=True)
                results_lists.append([d, i] + res[:-1])
        else:
            for i in d_sets[d][6]:
                print(i)
                res = utils.misidentify_metabolites_by_formula(i, d_sets[d][4], d_sets[d][2], compounds,
                                                               d_sets[d][3], zamboni=False)
                results_lists.append([d, i] + res[:-1])

    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage misidentification", "n_p_less_0.1", "n_q_less_0.1", "p_std",
                                   "q_std", "sig_paths"])
    res_df.to_csv("Metabolite_misidentification_by_formula_simulation_Reactome.csv")

    simulation_res = res_df
    print(simulation_res.head)
    plt.figure(figsize=(5, 4), dpi=600)
    plt.style.use("seaborn")
    for i in d_sets.keys():
        if i in ["Quirós", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                     simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                     label=i, fmt='o', linestyle="--", capsize=5, markeredgewidth=2, markersize=4)
        else:
            plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage misidentification'],
                         simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                         yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                         label=i, fmt='o', linestyle="solid", capsize=5, markeredgewidth=2, markersize=4)

    plt.title("Misidentification by formula", fontsize=14)
    plt.legend(fontsize=11)
    plt.ylabel("Mean number of pathways significant at P < 0.1", fontsize=13)
    plt.xlabel("Percentage of metabolites misidentified", fontsize=13)
    plt.savefig("metabolite_misidentification_by_formula_KEGG.png", dpi=600)
    plt.show()


# misidentification_formula_plot(db="KEGG")

def misidentification_barplot(pg, db="KEGG"):

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
    for d in d_sets.keys():
        print(d)
        if d.startswith("Fuhrer"):
            original_pathways = pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                                      zamboni=True)[4][0]

            for i in [4]:
                print(i)
                res = pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]],
                            zamboni=True)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways)  # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = 1-(number_common_paths / total_significant_paths)
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len(
                        [i for i in x if i not in original_pathways]) / total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                sem_TPR = np.std(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                sem_FPR = np.std(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR, sem_TPR])
                results_FPR.append([d, i, avg_fraction_FPR, sem_FPR])

        else:
            original_pathways = \
                pg[0](0, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4][
                    0]
            for i in [4]:
                print(i)
                res = \
                pg[0](i, d_sets[d][pg[1]], d_sets[d][pg[2]], d_sets[d][pg[3]], d_sets[d][pg[4]], zamboni=False)[4]
                misidentified_pathways = res
                pathway_fractions_TPR = []
                pathway_fractions_FPR = []
                for x in misidentified_pathways:
                    total_significant_paths = len(original_pathways)  # True positive + false positive
                    number_common_paths = len([i for i in x if i in original_pathways])
                    fraction_pathways_TPR = 1-(number_common_paths / total_significant_paths)
                    pathway_fractions_TPR.append(fraction_pathways_TPR)
                    fraction_pathways_FPR = len(
                        [i for i in x if i not in original_pathways]) / total_significant_paths
                    pathway_fractions_FPR.append(fraction_pathways_FPR)
                avg_fraction_TPR = np.mean(pathway_fractions_TPR)
                sem_TPR = stats.sem(pathway_fractions_TPR)
                avg_fraction_FPR = np.mean(pathway_fractions_FPR)
                sem_FPR = stats.sem(pathway_fractions_FPR)
                results_TPR.append([d, i, avg_fraction_TPR, sem_TPR])
                results_FPR.append([d, i, avg_fraction_FPR, sem_FPR])
    res_df_TPR = pd.DataFrame(results_TPR,
                              columns=["Dataset", "Percentage misidentification", "Average fraction", "SEM"])
    res_df_FPR = pd.DataFrame(results_FPR,
                              columns=["Dataset", "Percentage misidentification", "Average fraction", "SEM"])

    return res_df_TPR, res_df_FPR


mass_TPR, mass_FPR = misidentification_barplot(param_grid_heatmaps["mass"], db="KEGG")
formula_TPR, formula_FPR = misidentification_barplot(param_grid_heatmaps["formula"], db="KEGG")

plt.style.use("seaborn-darkgrid")
plt.rc('xtick', labelsize=11)

fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 9), sharex=True, sharey=True)

ax1.bar(mass_FPR["Dataset"].tolist(), mass_FPR["Average fraction"],
        yerr=mass_FPR["SEM"], capsize=5, label="Pathway gain rate", color="cornflowerblue")
ax1.bar(mass_TPR["Dataset"].tolist(), -mass_TPR["Average fraction"],
        yerr=mass_TPR["SEM"], capsize=5, label="Pathway loss rate", color="tomato")
ax1.set_title("Misidentification by mass", fontsize=13)

ax2.bar(formula_FPR["Dataset"].tolist(), formula_FPR["Average fraction"],
        yerr=formula_FPR["SEM"], capsize=5, label="Pathway gain rate", color="cornflowerblue")
ax2.bar(formula_TPR["Dataset"].tolist(), -formula_TPR["Average fraction"],
        yerr=formula_TPR["SEM"], capsize=5, label="Pathway loss rate", color="tomato")
ax2.set_title("Misidentification by chemical formula", fontsize=13)

ax2.set_xlabel("Dataset", fontsize=13)
fig.text(0.06, 0.25, "Pathway gain (upper bars) and pathway loss (lower bars) rate", fontsize=13, ha='center',
         rotation="vertical")
ax1.legend(fontsize=11)
for tick in ax2.get_xticklabels():
    tick.set_rotation(45)

ticks = ax1.get_yticks()
ticks2 = ax2.get_yticks()

# set labels to absolute values and with integer representation
ax1.set_yticklabels([round((abs(x)), 1) for x in ax1.get_yticks()])

plt.savefig("pathway_gain_loss_mass_formula_4_pct_2.png", dpi=600)
plt.show()



