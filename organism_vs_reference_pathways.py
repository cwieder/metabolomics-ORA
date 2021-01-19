import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Import the relevant datasets
DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

# Import pathway sets
KEGG_reference_pathways = pd.read_csv("KEGG_reference_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

# param grid
datasets = {"Auwerx": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx, [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)]],
            "Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]]}

print("Data import complete")

def organism_vs_reference(db="KEGG"):
    d_sets = datasets
    # if db == "Reactome":
    #     d_sets = datasets_reactome
    # if db == "Cyc":
    #     d_sets = datasets_biocyc
    plt_dict = {}
    for i in d_sets.keys():
        ora_res_org = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
        ora_res_ref = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], KEGG_reference_pathways)
        intersect = (set(ora_res_org["Pathway_ID"].str.slice(start=-5).tolist()) & set(ora_res_ref["Pathway_ID"].str.slice(start=-5).tolist()))
        # Ensures pathways are the same in both results (whole background can have additional pathways)
        ora_res_ref = ora_res_ref[ora_res_ref["Pathway_ID"].str.slice(start=-5).isin(intersect)]
        ora_res_org = ora_res_org[ora_res_org["Pathway_ID"].str.slice(start=-5).isin(intersect)]
        ora_res_org_pvals = np.negative(np.log10(ora_res_org["P-value"].tolist()))
        ora_res_ref_pvals = np.negative(np.log10(ora_res_ref["P-value"].tolist()))
        plt_dict[i] = [ora_res_org_pvals, ora_res_ref_pvals]
    for k, v in plt_dict.items():
        print(len(v[0]), len(v[1]))

    plt.figure(figsize=(6, 6))
    sns.set_style("darkgrid")
    sns.set_palette("muted")
    for i in plt_dict.keys():
        x = plt_dict[i][0]
        y = plt_dict[i][1]
        # jittered_y = y + 0.1 * np.random.rand(len(y)) - 0.05
        # jittered_x = x + 0.1 * np.random.rand(len(x)) - 0.05
        ax = sns.regplot(x=x, y=y,
                         ci=95,
                         scatter_kws={'s': 3})
    ax.set_xlabel("Organism-specific pathway set (-log10 P-value)",
                  fontsize=12)
    ax.set_ylabel("Reference pathway set (-log10 P-value)",
                  fontsize=12)
    ax.set(ylim=(0, 8), xlim=(0, 8))
    ax.legend(plt_dict.keys())
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', linestyle=':')
    ax.axhline(y=1, linewidth=1, color='black', linestyle='--')
    ax.axvline(x=1, linewidth=1, color='black', linestyle='--')
    plt.title("KEGG: Organism-specific vs. reference pathways")
    plt.savefig("../Figures/organism_vs_reference_pathways_pvals.png", dpi=300)
    plt.show()

# organism_vs_reference()

def organism_vs_reference_pathways_bar(db="KEGG"):
    dataframes = []
    d_sets = datasets
    # if db == "Reactome":
    #     d_sets = datasets_reactome
    # if db == "Cyc":
    #     d_sets = datasets_biocyc
    for i in d_sets.keys():
        ora_res_org = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], d_sets[i][2])
        ora_res_ref = utils.over_representation_analysis(d_sets[i][0], d_sets[i][1], KEGG_reference_pathways)
        intersect = (set(ora_res_org["Pathway_ID"].str.slice(start=-5).tolist()) & set(ora_res_ref["Pathway_ID"].str.slice(start=-5).tolist()))
        # Ensures pathways are the same in both results (whole background can have additional pathways)
        ora_res_ref = ora_res_ref[ora_res_ref["Pathway_ID"].str.slice(start=-5).isin(intersect)]
        ora_res_org = ora_res_org[ora_res_org["Pathway_ID"].str.slice(start=-5).isin(intersect)]

        n_p_less_01_ref = len(ora_res_ref[ora_res_ref["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01_ref = len(ora_res_ref[ora_res_ref["P-adjust"] < 0.1]["P-adjust"].tolist())
        n_p_less_01_org = len(ora_res_org[ora_res_org["P-value"] < 0.1]["P-value"].tolist())
        n_q_less_01_org = len(ora_res_org[ora_res_org["P-adjust"] < 0.1]["P-adjust"].tolist())
        df = pd.DataFrame([[n_p_less_01_ref, n_q_less_01_ref], [n_p_less_01_org, n_q_less_01_org]],
                          index=["General pathway set", "Organism-specific pathway set"], columns=["P", "Q"])
        df["Name"] = "df" + i
        dataframes.append(df)

    dfall = pd.concat([pd.melt(i.reset_index(),
                               id_vars=["Name", "index"])  # transform in tidy format each df
                       for i in dataframes],
                      ignore_index=True)

    dfall.set_index(["Name", "index", "variable"], inplace=True)
    dfall["vcs"] = dfall.groupby(level=["Name", "index"]).cumsum()
    dfall.reset_index(inplace=True)
    sns.set_style("dark")
    sns.set_palette("muted")
    plt.figure(figsize=(7, 5))
    for i, g in enumerate(dfall.groupby("variable")):
        ax = sns.barplot(data=g[1],
                         x="index",
                         y="vcs",
                         hue="Name",
                         zorder=-i,  # so first bars stay on top
                         edgecolor="k")
    ax.set_xlabel('Pathway set used')
    ax.set_ylabel('Number of significant pathways at \n P < 0.1 (solid bars) and Q < 0.1 (hatched bars)')
    labels = ["Auwerx", "Yamada", "Stevens", "Brown", "Zamboni (yfgM)", "Zamboni (dcuS)"]
    h, l = ax.get_legend_handles_labels()
    plt.legend(h[0:6], labels, title="Dataset", bbox_to_anchor=(1.4, 1), loc="upper right")
    # Set hatches for q-values bars
    plt.subplots_adjust(right=0.75)
    bars = ax.patches
    for i in range(12, 24, 1):
        bars[i].set_hatch('//')
    plt.title("KEGG: Organism-specific vs. general pathway set", fontsize=14)
    plt.tight_layout()
    plt.savefig("../Figures/pathway_set_barchart.png", dpi=300)
    plt.show()

# organism_vs_reference_pathways_bar(db="KEGG")

# Calculate pathway stats
