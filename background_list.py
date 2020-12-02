# Simulation of background list changes and reduction
import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Import the relevant datasets
DEM_yamada, background_yamada = process_datasets.yamada_data()
DEM_stevens, background_stevens = process_datasets.stevens_data()
DEM_brown, background_brown = process_datasets.brown_data()
DEM_yfgM, background_yfgM = process_datasets.zamboni_data("yfgM")
DEM_dcuS, background_dcuS = process_datasets.zamboni_data("dcuS")
# Import pathway set
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)

# Using default KEGG background set (whole KEGG)

# Reducing background set
datasets = {"Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways],
            "Zamboni (dcuS)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways]}
#
percentage_reductions = [0, 1, 5, 10, 20, 50, 70]

results_lists = []
for d in datasets.keys():
    print(d)
    for i in percentage_reductions:
        res = utils.reduce_background_list_ora(datasets[d][1], i, datasets[d][0], datasets[d][2])
        results_lists.append([d, i] + res)

res_df = pd.DataFrame(results_lists, columns=["Dataset", "Percentage reduction", "n_p_less_0.1", "p_std", "n_q_less_0.1", "q_std"])
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

plt.show()
# Mind the gap set