# Simulation of background list changes and reduction
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import gridspec
import matplotlib.patches as mpatches

# Import the relevant datasets
# Import Reactome datasets
# DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="Reactome")
# DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="Reactome")
# DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="Reactome")
# DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="Reactome")
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="Reactome")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="Reactome")

# Import KEGG datasets
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

# Import BioCyc datasets
# DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="Cyc")
# DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="Cyc")
# DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="Cyc")
# DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="Cyc")
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="Cyc")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="Cyc")

# Import KEGG pathway sets
KEGG_human_pathways = pd.read_csv("data/KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("data/KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("data/KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import Reactome pathway sets
Reactome_pathways = pd.read_csv("data/Reactome_pathway_set.csv", dtype=str, index_col=0)
Reactome_human_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("HSA")]
Reactome_eco_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("ECO")]
Reactome_mouse_pathways = Reactome_pathways[Reactome_pathways.index.str.contains("MMU")]
all_reactome_human_bg = list(set([x for x in Reactome_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_reactome_mouse_bg = list(set([x for x in Reactome_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# Import BioCyc pathway sets
BioCyc_human_pathways = pd.read_csv("data/Metacyc_human_pathways.csv")
BioCyc_eco_pathways = pd.read_csv("data/Metacyc_EColi_pathways.csv")
all_biocyc_human_bg = list(set([x for x in BioCyc_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_biocyc_eco_bg = list(set([x for x in BioCyc_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Labbé": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown],
            "Yachida": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Quirós": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx],
            "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM],
            "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS]}

# datasets_reactome = {"Quirós": [DEM_auwerx, background_auwerx, Reactome_human_pathways, all_reactome_human_bg],
#                      "Yachida": [DEM_yamada, background_yamada, Reactome_human_pathways, all_reactome_human_bg],
#                      "Stevens": [DEM_stevens, background_stevens, Reactome_human_pathways, all_reactome_human_bg],
#                      "Labbé": [DEM_brown, background_brown, Reactome_mouse_pathways, all_reactome_mouse_bg],
#                      "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, Reactome_human_pathways, all_reactome_human_bg],
#                      "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, Reactome_human_pathways, all_reactome_human_bg]}
#
# datasets_biocyc = {"Quirós": [DEM_auwerx, background_auwerx, BioCyc_human_pathways, all_biocyc_human_bg],
#                    "Yachida": [DEM_yamada, background_yamada, BioCyc_human_pathways, all_biocyc_human_bg],
#                    "Stevens": [DEM_stevens, background_stevens, BioCyc_human_pathways, all_biocyc_human_bg],
#                    "Labbé": [DEM_brown, background_brown, BioCyc_human_pathways, all_biocyc_human_bg],
#                    "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, BioCyc_eco_pathways, all_biocyc_eco_bg],
#                    "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, BioCyc_eco_pathways, all_biocyc_eco_bg]}

print("Data processing complete.")


def plot_log_pvalues(db="KEGG"):
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    if db == "Cyc":
        d_sets = datasets_biocyc
    plt_dict = {}
    for i in d_sets.keys():
        ora_res = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
        ora_res_all = utils.over_representation_analysis(d_sets[i][0], d_sets[i][3], d_sets[i][2])
        intersect = (set(ora_res["Pathway_ID"].tolist()) & set(ora_res_all["Pathway_ID"].tolist()))
        # Ensures pathways are the same in both results (whole background can have additional pathways)
        ora_res_all = ora_res_all[ora_res_all["Pathway_ID"].isin(intersect)]
        ora_res_pvals = np.negative(np.log10(ora_res["P-value"].tolist()))
        ora_res_all_pvals = np.negative(np.log10(ora_res_all["P-value"].tolist()))
        plt_dict[i] = [ora_res_pvals, ora_res_all_pvals]

    plt.figure(figsize=(6, 6), dpi=600)
    sns.set_style("darkgrid")
    sns.set_palette("muted")
    for i in plt_dict.keys():
        x = plt_dict[i][0]
        y = plt_dict[i][1]
        # jittered_y = y + 0.1 * np.random.rand(len(y)) - 0.05
        # jittered_x = x + 0.1 * np.random.rand(len(x)) - 0.05
        ax = sns.regplot(x=x, y=y,
                         ci=95,
                         scatter_kws={'s': 5})
    # set dotted lines for TOF-MS data
    ax.lines[3].set_linestyle("--")
    ax.lines[4].set_linestyle("--")
    ax.lines[5].set_linestyle("--")
    ax.set_xlabel("Specified background list (-log10 P-value)",
                  fontsize=13)
    ax.set_ylabel("All " + "KEGG" + " compounds (organism-specific) (-log10 P-value)",
                  fontsize=13)
    ax.set(ylim=(0, 10), xlim=(0, 10))
    ax.legend(plt_dict.keys(), fontsize=11)
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', linestyle=':')
    ax.axhline(y=1, linewidth=1, color='black', linestyle='--')
    ax.axvline(x=1, linewidth=1, color='black', linestyle='--')
    # plt.title("KEGG")
    # plt.savefig("../Figures/logp_plot_KEGG.png", dpi=600)
    plt.show()

    # fig, ax = plt.subplots(3,2)
    # l = list(plt_dict.keys())
    # ax = ax.flatten()
    # colours = sns.color_palette("muted")
    # for i in range(0, 6):
    #     plt.sca(ax[i])
    #     sns.regplot(x=plt_dict[l[i]][0], y=plt_dict[l[i]][1],
    #                  ci=None,
    #                  scatter_kws={'s': 5}, color=colours[i])
    #     lims = [
    #         np.min([plt.xlim(), plt.ylim()]),  # min of both axes
    #         np.max([plt.xlim(), plt.ylim()]),  # max of both axes
    #     ]
    #
    #     # now plot both limits against eachother
    #     plt.plot(lims, lims, 'k-', alpha=0.75, zorder=0, linestyle=':', color='black')
    #     plt.title(l[i])
    # plt.tight_layout()
    # # plt.savefig("datasets_log_p_subplots.png", dpi = 300)
    # plt.show()


plot_log_pvalues(db="KEGG")


def plot_grouped_stacked_bar(db="KEGG"):
    dataframes = []
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    if db == "Cyc":
        d_sets = datasets_biocyc
    for i in d_sets.keys():
        ora_res = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
        ora_res_all = utils.over_representation_analysis(d_sets[i][0], d_sets[i][3], d_sets[i][2])
        intersect = (set(ora_res["Pathway_ID"].tolist()) & set(ora_res_all["Pathway_ID"].tolist()))
        # Ensures pathways are the same in both results (whole background can have additional pathways)
        ora_res_all = ora_res_all[ora_res_all["Pathway_ID"].isin(intersect)]
        n_p_less_01 = len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01 = len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist())
        n_p_less_01_all = len(ora_res_all[ora_res_all["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01_all = len(ora_res_all[ora_res_all["P-adjust"] < 0.1]["P-adjust"].tolist())
        df = pd.DataFrame([[n_q_less_01, n_p_less_01], [n_q_less_01_all, n_p_less_01_all]],
                          index=["Specified background list", "All " + db + " compounds"], columns=["P", "Q"])
        df["Name"] = i
        dataframes.append(df)

    dfall = pd.concat([pd.melt(i.reset_index(),
                               id_vars=["Name", "index"])  # transform in tidy format each df
                       for i in dataframes],
                      ignore_index=True)

    dfall.set_index(["Name", "index", "variable"], inplace=True)
    dfall["vcs"] = dfall.groupby(level=["Name", "index"]).cumsum()
    dfall.reset_index(inplace=True)
    print(dfall)
    sns.set_style("dark")
    sns.set_palette("muted")
    plt.figure(figsize=(6, 6), dpi=600)
    for i, g in enumerate(dfall.groupby("variable")):
        ax = sns.barplot(data=g[1],
                         x="Name",
                         y="vcs",
                         hue="index",
                         zorder=-i,  # so first bars stay on top
                         edgecolor="k",
                         ci=None)
    ax.set_xlabel('Background list used in ORA', fontsize=13)
    # ax.set_xticks(x, labels, rotation='vertical')
    ax.set_ylabel('Number of significant pathways at \n P < 0.1 (solid bars) and Q < 0.1 (hatched bars)', fontsize=13)
    # labels = ["Quirós", "Yachida", "Stevens", "Labbé", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]

    # Set hatches for q-values bars
    # plt.subplots_adjust(right=0.8)
    plt.xticks(rotation=45)
    bars = ax.patches
    print(len(bars))
    for i in range(6, 12):
        bars[i].set_color(sns.color_palette("deep", 6)[i - 6])
        bars[i].set_edgecolor("k")
    for i in range(0, 6):
        bars[i].set_color(sns.color_palette("pastel", 6)[i])
        bars[i].set_edgecolor("k")
    for i in range(12, 18):
        bars[i].set_color(sns.color_palette("pastel", 6)[i - 12])
        bars[i].set_edgecolor("k")
    for i in range(18, 24):
        bars[i].set_color(sns.color_palette("deep", 6)[i - 18])
        bars[i].set_edgecolor("k")
    for i in range(0, 12, 1):
        bars[i].set_hatch('//')

    specified_patch = mpatches.Patch(color='darkgray', label='Light colours: specified background set')
    unspecified_patch = mpatches.Patch(color='dimgray', label='Dark colours: all KEGG compounds')
    hatched = mpatches.Patch(facecolor='white', hatch="//", label='Hatched bars: pathways significant at Q ≤ 0.1',
                             edgecolor='k', alpha=0.5)
    solid = mpatches.Patch(facecolor='white', label='Solid bars: pathways significant at P ≤ 0.1', edgecolor='k',
                           alpha=0.5)
    plt.legend(handles=[specified_patch, unspecified_patch, hatched, solid], fontsize=11)

    # plt.title(db, fontsize=14)
    plt.tight_layout()
    # plt.savefig("../Figures/all_vs_experimental_barchart_KEGG.png", dpi=600)
    plt.show()


# plot_grouped_stacked_bar(db="KEGG")

# Reducing background set
def reduce_background_set(db="KEGG"):
    percentage_reductions_keep_DEM = [i for i in range(100, 5, -5)]
    percentage_reductions = [i for i in range(100, 5, -5)]
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    results_lists_keep_DEM = []
    results_lists = []
    for d in d_sets.keys():
        print(d)
        if d.startswith("Fuhrer"):
            for i in percentage_reductions:
                res = utils.reduce_background_list_ora(d_sets[d][1], d_sets[d][4], i, d_sets[d][0], d_sets[d][2], keep_DEM=False, Zamboni=True)
                results_lists.append([d, i] + res)
            for i in percentage_reductions_keep_DEM:
                print(i)
                res_keep_DEM = utils.reduce_background_list_ora(d_sets[d][1], d_sets[d][4], i, d_sets[d][0], d_sets[d][2], keep_DEM=True, Zamboni=True)
                results_lists_keep_DEM.append([d, i] + res_keep_DEM)
        else:
            for i in percentage_reductions:
                res = utils.reduce_background_list_ora(d_sets[d][1], d_sets[d][4], i, d_sets[d][0], d_sets[d][2], keep_DEM=False, Zamboni=False)
                results_lists.append([d, i] + res)
            for i in percentage_reductions_keep_DEM:
                print(i)
                res_keep_DEM = utils.reduce_background_list_ora(d_sets[d][1], d_sets[d][4], i, d_sets[d][0], d_sets[d][2], keep_DEM=True, Zamboni=False)
                results_lists_keep_DEM.append([d, i] + res_keep_DEM)



    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage reduction", "n_p_less_0.1",
                                   "n_q_less_0.1", "mean_proportion_p_vals", "p_std",
                                   "q_std", "sd_proportion_p_vals"])
    res_df_keep_DEM = pd.DataFrame(results_lists_keep_DEM,
                                   columns=["Dataset", "Percentage reduction", "n_p_less_0.1",
                                            "n_q_less_0.1", "mean_proportion_p_vals", "p_std",
                                            "q_std", "sd_proportion_p_vals"])
    res_df_keep_DEM.to_csv("Background_reduction_simulation_keep_DEM.csv")
    res_df.to_csv("Background_reduction_simulation.csv")
    # simulation_res = pd.read_csv("Background_reduction_simulation.csv")
    # simulation_res_keep_DEM = pd.read_csv("Background_reduction_simulation_keep_DEM.csv")
    simulation_res = res_df
    simulation_res_keep_DEM = res_df_keep_DEM
    with plt.style.context('seaborn-darkgrid'):
        sns.set_palette("deep")
        fig = plt.figure(figsize=(10, 6), dpi=600)
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.set_title("Random background list reduction", fontsize=14)
        for i in d_sets.keys():
            if i in ["Quirós", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]:
                ax1.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage reduction'],
                             simulation_res[simulation_res["Dataset"] == i]['mean_proportion_p_vals'],
                             yerr=simulation_res[simulation_res["Dataset"] == i]['sd_proportion_p_vals'],
                             label=i, fmt='o', linestyle="--", capsize=5, markeredgewidth=2, markersize=4)
            else:
                ax1.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage reduction'],
                             simulation_res[simulation_res["Dataset"] == i]['mean_proportion_p_vals'],
                             yerr=simulation_res[simulation_res["Dataset"] == i]['sd_proportion_p_vals'],
                             label=i, fmt='o', linestyle="solid", capsize=5, markeredgewidth=2, markersize=4)
        ax1.set_xlim(100, 10)
        ax2.set_title("No DA metabolite removal", fontsize=14)
        for i in d_sets.keys():
            if i in ["Quirós", "Fuhrer (yfgM)", "Fuhrer (dcuS)"]:
                ax2.errorbar(simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i]['Percentage reduction'],
                             simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i]['mean_proportion_p_vals'],
                             yerr=simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i][
                                 'sd_proportion_p_vals'],
                             label=i, fmt='o', linestyle="--", capsize=5, markeredgewidth=2, markersize=4)
            else:
                ax2.errorbar(simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i]['Percentage reduction'],
                             simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i]['mean_proportion_p_vals'],
                             yerr=simulation_res_keep_DEM[simulation_res_keep_DEM["Dataset"] == i][
                                 'sd_proportion_p_vals'],
                             label=i, fmt='o', linestyle="solid", capsize=5, markeredgewidth=2, markersize=4)
        ax2.set_xlim(100, 10)
        # fig.suptitle("Reactome", fontsize=14)
        handles, labels = ax1.get_legend_handles_labels()
        # plt.subplots_adjust(right=0.7)
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)

        plt.ylabel("Proportion of pathways significant at P ≤ 0.1 \n compared to at baseline (original background list)",
                   fontsize=13)
        plt.xlabel("Percentage of original compounds present in background list", fontsize=13)
        plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=11)
        plt.tight_layout()
        plt.savefig("background_list_reduction_proportion_KEGG_new2.png", dpi=600)
        plt.show()


# reduce_background_set(db="KEGG")

# Mind the gap set
