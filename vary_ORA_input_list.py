import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle



# Import the relevant datasets
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

# param grid
datasets = {"Auwerx": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx, [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)]],
            "Yamada": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Brown": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Zamboni (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]],
            "Zamboni (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]]}

print("Data import complete")

def vary_p_val_cutoff():
    p_value_cutoffs = [0.001, 0.01, 0.05, 0.1, 0.2]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in datasets.keys():
        for p in p_value_cutoffs:
            for m in multiple_test_options:
                t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                DA_metabolites = t_test_res[t_test_res["P-adjust"] < p]["Metabolite"].tolist()
                ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                res_list.append([d, p, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    with plt.style.context('seaborn'):
        fig, ax1 = plt.subplots(1, 1)
        plt.style.use("seaborn")
        cols = sns.color_palette("deep", 8)

        for num, i in enumerate(datasets.keys()):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_numpy().astype('str'),
                         res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" Bonferroni", linestyle='-', color=cols[num])
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_numpy().astype('str'),
                         res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" FDR BH", linestyle='--', color=cols[num])
        # plt.title("Reactome", fontsize=14)

        ax1.legend()
        ax1.set_ylabel("Number of pathways significant at P < 0.1")
        ax1.set_xlabel("P-value cutoff for differentially abundant metabolite list")
        plt.savefig("../Figures/KEGG_p_value_cutoffs.png", dpi=300)
        plt.show()

# vary_p_val_cutoff()

def vary_dam_size():
    proportion_of_bg = [i for i in range(0, 85, 5)]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []

    # Import files for zamboni data
    n_zscore = pd.read_csv("../Zamboni/mod_zscore_neg_CW.csv", index_col=0)
    p_zscore = pd.read_csv("../Zamboni/mod_zscore_pos_CW.csv", index_col=0)
    # remove unannotated
    n_zscore = n_zscore[n_zscore.index.notnull()]
    p_zscore = p_zscore[p_zscore.index.notnull()]
    with open('zamboni_pos_annotation_dict.pickle', 'rb') as handle:
        annotations_pos = pickle.load(handle)
    with open('zamboni_neg_annotation_dict.pickle', 'rb') as handle:
        annotations_neg = pickle.load(handle)

    for d in datasets.keys():
        print(len(datasets[d][1]))
        pass
        print(d)
        for p in proportion_of_bg:
            print(p)
            if d.startswith("Zamboni"):
                strain = d[d.find("(")+1:d.find(")")]
                neg_col = n_zscore.loc[:, strain].to_frame(name="z-score")
                neg_col["annotation"] = neg_col.index.map(annotations_neg)
                pos_col = p_zscore.loc[:, strain].to_frame(name="z-score")
                pos_col["annotation"] = pos_col.index.map(annotations_pos)
                all_modes = pd.concat([neg_col, pos_col])
                all_modes["z-score"] = all_modes["z-score"].abs()
                all_modes_sorted = all_modes.sort_values(by=['z-score'], ascending=False)
                DA_ions = all_modes_sorted.head(int(len(all_modes_sorted) * (p / 100)))
                DA_annotations = DA_ions["annotation"].tolist()
                DA_annotations = [j for i in DA_annotations for j in i]
                DA_metabolites = list(filter(None, list(set(DA_annotations))))
                ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                res_list.append([d, p, "none", len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])

            else:
                for m in multiple_test_options:
                    t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                    t_test_res = t_test_res.sort_values(by=["P-adjust"])
                    DA_metabolites = t_test_res.head(int(len(t_test_res) * (p/100)))["Metabolite"].tolist()
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, p, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    res_df_Zamboni = res_df[res_df["Multiple_correction_method"] == "none"]
    with plt.style.context('seaborn'):
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 8))
        plt.style.use("seaborn")
        cols = sns.color_palette("deep", 8)

        for num, i in enumerate(["Auwerx", "Yamada", "Stevens", "Brown"]):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_list(),
                         res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" Bonferroni", linestyle='-', color=cols[num])
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_list(),
                         res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" FDR BH", linestyle='--', color=cols[num])
        for num, i in enumerate(["Zamboni (dcuS)", "Zamboni (yfgM)"]):
            ax1.plot(res_df_Zamboni[res_df_Zamboni["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_Zamboni[res_df_Zamboni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i,
                     linestyle='-', color=cols[4+num])

        ax1.legend()
        ax1.set_ylabel("Number of pathways significant at P < 0.1")
        ax1.set_xlabel("Top N% of differentially abundant metabolites (ranked by Q-value)")
        plt.savefig("../Figures/vary_DAM_size_percentage.png", dpi=300)
        plt.show()

vary_dam_size()