import pickle
import pandas as pd
import utils
import process_datasets
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import numpy as np
logp_all = pd.read_csv("hmdb_logp_all.csv", index_col=0)

DEM_auwerx, background_auwerx, mat_auwerx = process_datasets.auwerx_data(db="KEGG")
DEM_yamada, background_yamada, mat_yamada = process_datasets.yamada_data(db="KEGG")
DEM_stevens, background_stevens, mat_stevens = process_datasets.stevens_data(db="KEGG")
DEM_brown, background_brown, mat_brown = process_datasets.brown_data(db="KEGG")
DEM_yfgM, background_yfgM, mat_yfgM = process_datasets.zamboni_data("yfgM", db="KEGG")
DEM_dcuS, background_dcuS, mat_dcuS = process_datasets.zamboni_data("dcuS", db="KEGG")

KEGG_mouse_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_human_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
KEGG_eco_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)


cols = mat_brown.columns.tolist()
print(logp_all.columns)


matching_id = logp_all[logp_all["kegg_id"].isin(cols)]

print(np.median(matching_id['logp']))
print(np.mean(matching_id['logp']))
print(np.std(matching_id['logp']))
print(np.min(matching_id['logp']))
print(np.max(matching_id['logp']))

plt.style.use("seaborn")
plt.hist(matching_id["logp"], bins=20)
plt.xlabel("LogP partition coefficient of KEGG compounds")
plt.ylabel("Frequency")
plt.show()

# Brown - split at -1 or -2
def log_p_venn(cutoff, mat, pathways):
    polar = matching_id[matching_id["logp"] < cutoff]
    nonpolar = matching_id[matching_id["logp"] > cutoff]

    polar_mat = mat.filter(polar["kegg_id"])
    polar_mat = polar_mat.iloc[:, ~polar_mat.columns.duplicated()]
    ttest_polar = utils.t_tests(polar_mat, mat["Group"], "fdr_bh")
    dem_polar = ttest_polar[ttest_polar["P-adjust"] <= 0.05]["Metabolite"].tolist()
    ora_polar = utils.over_representation_analysis(dem_polar, polar_mat.columns.tolist(), pathways)

    nonpolar_mat = mat.filter(nonpolar["kegg_id"])
    nonpolar_mat = nonpolar_mat.iloc[:, ~nonpolar_mat.columns.duplicated()]
    ttest_nonpolar = utils.t_tests(nonpolar_mat, mat["Group"], "fdr_bh")
    dem_nonpolar = ttest_nonpolar[ttest_nonpolar["P-adjust"] <= 0.05]["Metabolite"].tolist()
    ora_nonpolar = utils.over_representation_analysis(dem_nonpolar, nonpolar_mat.columns.tolist(), pathways)

    print(ora_polar[ora_polar["P-value"] <= 0.1])
    print(ora_nonpolar[ora_nonpolar["P-value"] <= 0.1])

    polar_pathways = ora_polar[ora_polar["P-value"] <= 0.1]["Pathway_ID"].tolist()
    nonpolar_pathways = ora_nonpolar[ora_nonpolar["P-value"] <= 0.1]["Pathway_ID"].tolist()
    common_paths = list(set(polar_pathways) & set(nonpolar_pathways))

    pp_names = [' '.join(i.split(' ')[:-4]) for i in ora_polar[ora_polar["P-value"] <= 0.1]["Pathway_name"].tolist()]
    np_names = [' '.join(i.split(' ')[:-4]) for i in ora_nonpolar[ora_nonpolar["P-value"] <= 0.1]["Pathway_name"].tolist()]
    common_names = list(set(pp_names) & set(np_names))

    plt.style.use("seaborn")
    venn = venn2(subsets = (set(polar_pathways),
                            set(nonpolar_pathways)), set_labels = ('Significant pathways:\n polar compounds', 'Significant pathways:\n non-polar compounds'),
                 set_colors=('tab:blue', 'tab:orange'), alpha = 0.3)
    venn.get_label_by_id('100').set_text("\n\n".join([i for i in pp_names if i not in common_names]))
    # venn.get_label_by_id('110').set_text("\n".join(common_names))
    venn.get_label_by_id('010').set_text("\n\n".join([i for i in np_names if i not in common_names]))
    # venn2_circles(subsets = (set(polar_pathways),
    #                         set(nonpolar_pathways)),
    #               linestyle='solid', linewidth=0.5, color='k')
    plt.annotate("\n".join(common_names), xy=venn.get_label_by_id('110').get_position() - np.array([-0.01, -0.25]), xytext=(100,70),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray', linewidth=1))

    plt.tight_layout()
    # plt.savefig("../Figures/logp_brown_venn.png", dpi=300)
    plt.show()

log_p_venn(-2, mat_yamada, KEGG_human_pathways)

def log_p_venn_zamboni(cutoff, mat, pathways):
    polar = matching_id[matching_id["logp"] < cutoff]
    nonpolar = matching_id[matching_id["logp"] > cutoff]

    polar_mat = mat[mat.columns.intersection(polar["kegg_id"])]
    DEM_polar = []
    for x in polar_mat.T.itertuples():
        if x[1] > 6 or x[1] < -6:
            DEM_polar.append(x[0])
    ora_polar = utils.over_representation_analysis(DEM_polar, polar_mat.columns.tolist(), pathways)
    print(ora_polar[ora_polar["P-value"] <= 0.1])

    nonpolar_mat = mat[mat.columns.intersection(nonpolar["kegg_id"])]
    DEM_nonpolar = []
    for x in nonpolar_mat.T.itertuples():
        if x[1] > 6 or x[1] < -6:
            DEM_nonpolar.append(x[0])
    ora_nonpolar = utils.over_representation_analysis(DEM_nonpolar, nonpolar_mat.columns.tolist(), pathways)
    print(ora_nonpolar[ora_nonpolar["P-value"] <= 0.1])

    polar_pathways = ora_polar[ora_polar["P-value"] <= 0.1]["Pathway_ID"].tolist()
    nonpolar_pathways = ora_nonpolar[ora_nonpolar["P-value"] <= 0.1]["Pathway_ID"].tolist()
    common_paths = list(set(polar_pathways) & set(nonpolar_pathways))

    pp_names = [' '.join(i.split(' ')[:-5]) for i in ora_polar[ora_polar["P-value"] <= 0.1]["Pathway_name"].tolist()]
    np_names = [' '.join(i.split(' ')[:-5]) for i in ora_nonpolar[ora_nonpolar["P-value"] <= 0.1]["Pathway_name"].tolist()]
    common_names = list(set(pp_names) & set(np_names))

    plt.style.use("seaborn")
    venn = venn2(subsets = (set(polar_pathways),
                            set(nonpolar_pathways)), set_labels = ('Significant pathways:\n polar compounds', 'Significant pathways:\n non-polar compounds'),
                 set_colors=('tab:blue', 'tab:orange'), alpha = 0.3)
    venn.get_label_by_id('100').set_text("\n".join([i for i in pp_names if i not in common_names]))
    # venn.get_label_by_id('110').set_text("\n".join(common_names))
    venn.get_label_by_id('010').set_text("\n".join([i for i in np_names if i not in common_names]))
    # venn2_circles(subsets = (set(polar_pathways),
    #                         set(nonpolar_pathways)),
    #               linestyle='solid', linewidth=0.5, color='k')
    # plt.annotate("\n".join(common_names), xy=venn.get_label_by_id('110').get_position() - np.array([-0.01, -0.25]), xytext=(100,70),
    # ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    # arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray', linewidth=1))

    plt.tight_layout()
    # plt.savefig("../Figures/logp_brown_venn.png", dpi=600)
    plt.show()

# log_p_venn_zamboni(0, mat_brown, KEGG_mouse_pathways)