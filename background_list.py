# Simulation of background list changes and reduction
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Import the relevant datasets
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data()
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="Reactome")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data()
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="Reactome")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS")

# Import Reactome datasets

# Import KEGG pathway sets
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
all_reactome_eco_bg = list(set([x for x in Reactome_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_reactome_mouse_bg = list(set([x for x in Reactome_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Auwerx": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg],
            "Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg],
            # "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg]}

datasets_reactome = {"Yamada": [DEM_yamada, background_yamada, Reactome_human_pathways, all_reactome_human_bg],
                     "Brown": [DEM_brown, background_brown, Reactome_mouse_pathways, all_reactome_mouse_bg]}

def plot_log_pvalues(db="KEGG"):
    d_sets = datasets
    if db == "Reactome":
        d_sets = datasets_reactome
    plt_dict = {}
    for i in d_sets.keys():
        ora_res = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
        ora_res_all = utils.over_representation_analysis(d_sets[i][0], d_sets[i][3], d_sets[i][2])

        intersect = (set(ora_res["Pathway_ID"].tolist()) & set(ora_res_all["Pathway_ID"].tolist()))
        ora_res_all = ora_res_all[ora_res_all["Pathway_ID"].isin(intersect)]
        ora_res_pvals = np.negative(np.log10(ora_res["P-value"].tolist()))
        ora_res_all_pvals = np.negative(np.log10(ora_res_all["P-value"].tolist()))
        plt_dict[i] = [ora_res_pvals, ora_res_all_pvals]
    plt.figure(figsize=(6, 6))
    sns.set_style("darkgrid")
    sns.set_palette("muted")
    for i in plt_dict.keys():
        ax = sns.regplot(x=plt_dict[i][0], y=plt_dict[i][1],
                         ci=95,
                         scatter_kws={'s':5})
    ax.set_xlabel("Experimental background list (-log10 P-value)",
                  fontsize=12)
    ax.set_ylabel("All KEGG compounds (organism-specific) (-log10 P-value)",
                  fontsize=12)
    ax.set(ylim=(0, 10), xlim=(0, 10))
    ax.legend(plt_dict.keys())
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', linestyle=':')
    ax.axhline(y=1, linewidth=1, color='black',linestyle='--')
    ax.axvline(x=1, linewidth=1, color='black',linestyle='--')
    plt.show()
    #
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
    # plt.savefig("datasets_log_p_subplots.png", dpi = 300)
    # plt.show()

plot_log_pvalues(db="Reactome")

def plot_grouped_stacked_bar():
    dataframes = []
    for i in datasets.keys():
        ora_res = utils.over_representation_analysis(datasets[i][0], datasets[i][1], datasets[i][2])
        ora_res_all = utils.over_representation_analysis(datasets[i][0], datasets[i][3], datasets[i][2])
        n_p_less_01 = len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01 = len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist())
        n_p_less_01_all = len(ora_res_all[ora_res_all["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01_all = len(ora_res_all[ora_res_all["P-adjust"] < 0.1]["P-adjust"].tolist())
        df = pd.DataFrame([[n_p_less_01, n_q_less_01], [n_p_less_01_all, n_q_less_01_all]],
                          index=["Specified background list", "All KEGG compounds"], columns=["P", "Q"])
        df["Name"] = "df"+i
        dataframes.append(df)

    dfall = pd.concat([pd.melt(i.reset_index(),
                               id_vars=["Name", "index"]) # transform in tidy format each df
                       for i in dataframes],
                       ignore_index=True)

    dfall.set_index(["Name", "index", "variable"], inplace=True)
    dfall["vcs"] = dfall.groupby(level=["Name", "index"]).cumsum()
    dfall.reset_index(inplace=True)
    sns.set_style("dark")
    sns.set_palette("muted")
    for i, g in enumerate(dfall.groupby("variable")):
        ax = sns.barplot(data=g[1],
                         x="index",
                         y="vcs",
                         hue="Name",
                         zorder=-i, # so first bars stay on top
                         edgecolor="k")
    ax.set_xlabel('Background list used in ORA', fontsize=14)
    ax.set_ylabel('Number of significant pathways at \n P < 0.1 (solid bars) and Q < 0.1 (hatched bars)',
           fontsize=14)
    labels = ["Auwerx", "Yamada", "Stevens", "Brown", "Zamboni (yfgM)", "Zamboni (dcuS)"]
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[0:6], labels, title="Dataset")
    # Set hatches for q-values bars
    bars = ax.patches
    for i in range(12, 24, 1):
        bars[i].set_hatch('//')
    # plt.savefig("all_vs_experimental_barchart.png", dpi=300)
    plt.show()

# plot_grouped_stacked_bar()

# Reducing background set
def reduce_background_set():
    percentage_reductions = [i for i in range(100, 45, -5)]

    results_lists = []
    for d in datasets.keys():
        print(d)
        for i in percentage_reductions:
            res = utils.reduce_background_list_ora(datasets[d][1], i, datasets[d][0], datasets[d][2], keep_DEM=True)
            results_lists.append([d, i] + res)

    res_df = pd.DataFrame(results_lists,
                          columns=["Dataset", "Percentage reduction", "n_p_less_0.1",
                                   "n_q_less_0.1", "mean_proportion_p_vals", "p_std",
                                   "q_std", "sd_proportion_p_vals"])
    # res_df.to_csv("Background_reduction_simulation_keep_DEM.csv")

    # simulation_res = pd.read_csv("Background_reduction_simulation.csv")
    simulation_res = res_df
    print(simulation_res.head)
    plt.figure()
    plt.style.use("seaborn")
    for i in datasets.keys():
        plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage reduction'],
                     simulation_res[simulation_res["Dataset"] == i]['mean_proportion_p_vals'],
                     yerr=simulation_res[simulation_res["Dataset"] == i]['sd_proportion_p_vals'],
                     label=i, fmt='o', linestyle="solid", capsize=5,  markeredgewidth=2, markersize=4)
    # plt.title("Number of pathways with P-values < 0.1 in \n response to varying background list size", fontsize=14)
    plt.xlim(100, 50)
    # plt.legend()
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.subplots_adjust(right=0.7)
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylabel("Proportion of pathways significant at P < 0.1 \n compared to at baseline (original background set)")
    # plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=14)
    plt.xlabel("Percentage of original background list")
    plt.savefig("background_list_reduction_proportion_noDEMremoval.png", dpi=300)
    plt.show()

# reduce_background_set()
# Mind the gap set