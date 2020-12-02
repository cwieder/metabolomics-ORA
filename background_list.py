# Simulation of background list changes and reduction
import utils
import process_datasets
import pandas as pd

# Import the relevant datasets
DEM_yamada, background_yamada = process_datasets.yamada_data()
DEM_stevens, background_stevens = process_datasets.stevens_data()
DEM_brown, background_brown = process_datasets.brown_data()
print(len(background_brown))
print(len(background_stevens))
print(len(background_yamada))
# Import pathway set
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
# Usin

# Using default KEGG background set (whole KEGG)

# Reducing background set
datasets = {"yamada": [DEM_yamada, background_yamada],
            "stevens": [DEM_stevens, background_stevens],
            "brown": [DEM_brown, background_brown]}

percentage_reductions = [1]

results_lists = []
for d in ["yamada"]:
    print(d)
    for i in percentage_reductions:
        res = utils.reduce_background_list_ora(datasets[d][1], i, datasets[d][0], KEGG_human_pathways)
        results_lists.append([d, i] + res)

print(results_lists)

# Mind the gap set