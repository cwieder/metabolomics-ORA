# Simulation of background list changes and reduction
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Import the relevant datasets
DEM_yamada, background_yamada = process_datasets.yamada_data()
DEM_stevens, background_stevens = process_datasets.stevens_data()
DEM_brown, background_brown = process_datasets.brown_data()
DEM_yfgM, background_yfgM = process_datasets.zamboni_data("yfgM")
DEM_dcuS, background_dcuS = process_datasets.zamboni_data("dcuS")

# Import pathway sets
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
all_KEGG_human_bg = list(set([x for x in KEGG_human_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_eco_bg = list(set([x for x in KEGG_eco_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
all_KEGG_mouse_bg = list(set([x for x in KEGG_mouse_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))

datasets = {"Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg]}

dataframes = []
for i in datasets.keys():
    ora_res = utils.over_representation_analysis(datasets[i][0], datasets[i][1], datasets[i][2])
    ora_res_all = utils.over_representation_analysis(datasets[i][0], datasets[i][3], datasets[i][2])
    n_p_less_01 = len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())
    n_q_less_01 = len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist())
    n_p_less_01_all = len(ora_res_all[ora_res_all["P-value"] < 0.1]["P-value"].tolist())
    n_q_less_01_all = len(ora_res_all[ora_res_all["P-adjust"] < 0.1]["P-adjust"].tolist())
    df = pd.DataFrame([[n_p_less_01, n_q_less_01], [n_p_less_01_all, n_q_less_01_all]],
                      index=["Experimental", "All"], columns=["P", "Q"])
    df["Name"] = "df"+i
    dataframes.append(df)

dfall = pd.concat([pd.melt(i.reset_index(),
                           id_vars=["Name", "index"]) # transform in tidy format each df
                   for i in dataframes],
                   ignore_index=True)

dfall.set_index(["Name", "index", "variable"], inplace=True)
dfall["vcs"] = dfall.groupby(level=["Name", "index"]).cumsum()
dfall.reset_index(inplace=True)
c = ["blue", "purple", "red", "green", "pink"]
sns.set_theme()
for i, g in enumerate(dfall.groupby("variable")):
    ax = sns.barplot(data=g[1],
                     x="index",
                     y="vcs",
                     hue="Name",
                     color=c[i],
                     zorder=-i, # so first bars stay on top
                     edgecolor="k")
ax.legend_.remove() # remove the redundant legends
ax.set(xlabel='Background list used in ORA',
       ylabel='Number of significant pathways at \n P < 0.1 (lower bars) and Q < 0.1 (upper bars)')
plt.show()
quit()
# Reducing background set
percentage_reductions = [0, 1, 5, 10, 20, 50, 70]

results_lists = []
for d in datasets.keys():
    print(d)
    for i in percentage_reductions:
        res = utils.reduce_background_list_ora(datasets[d][1], i, datasets[d][0], datasets[d][2])
        results_lists.append([d, i] + res)

res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage reduction", "n_p_less_0.1", "n_q_less_0.1", "p_std", "q_std"])
res_df.to_csv("Background_reduction_simulation.csv")

# simulation_res = pd.read_csv("Background_reduction_simulation.csv")
simulation_res = res_df
print(simulation_res.head)
plt.figure()
plt.style.use("seaborn")
for i in datasets.keys():
    plt.errorbar(simulation_res[simulation_res["Dataset"] == i]['Percentage reduction'],
                 simulation_res[simulation_res["Dataset"] == i]['n_p_less_0.1'],
                 yerr=simulation_res[simulation_res["Dataset"] == i]['p_std'],
                 label=i, fmt='o', linestyle="solid", capsize=20, markersize=5)
plt.title("Number of pathways with P-values < 0.1 in \n response to varying background list size", fontsize=14)
plt.legend()
plt.ylabel("Mean number of pathways significant at P < 0.1 \n based on 100 random permutations", fontsize=14)
plt.xlabel("Percentage reduction of background list", fontsize=14)
plt.savefig("background_list_reduction.png", dpi=300)
plt.show()

# Mind the gap set