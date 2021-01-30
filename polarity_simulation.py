import pickle
import pandas as pd
import utils
import process_datasets
import matplotlib.pyplot as plt

logp_all = pd.read_csv("hmdb_logp_all.csv", index_col=0)
print(logp_all.shape)

# DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
# DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
# DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
# DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
# DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")


KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)


cols = mat_brown.columns.tolist()
print(logp_all.columns)


matching_id = logp_all[logp_all["kegg_id"].isin(cols)]

# plt.style.use("seaborn")
# plt.hist(matching_id["logp"], bins=20)
# plt.xlabel("LogP partition coefficient of KEGG compounds")
# plt.ylabel("Frequency")
# plt.show()

# Brown - split at -1 or -2
cutoff = -1
polar = matching_id[matching_id["logp"] < cutoff]
nonpolar = matching_id[matching_id["logp"] > cutoff]



polar_mat = mat_brown.filter(polar["kegg_id"])
polar_mat = polar_mat.iloc[:, ~polar_mat.columns.duplicated()]

ttest_polar = utils.t_tests(polar_mat, mat_brown["Group"], "fdr_bh")
dem_polar = ttest_polar[ttest_polar["P-adjust"] <= 0.05]["Metabolite"].tolist()

ora_polar = utils.over_representation_analysis(dem_polar, polar_mat.columns.tolist(), KEGG_mouse_pathways)
nonpolar_mat = mat_brown.filter(nonpolar["kegg_id"])

print(ora_polar[ora_polar["P-value"] <= 0.1])



