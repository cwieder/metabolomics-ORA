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
datasets = {"Quirós": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx, [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)]],
            "Yachida": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg, mat_stevens],
            "Labbé": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
            "Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]],
            "Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS, [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]]}

print("Data import complete")

def vary_dam_size():
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in ["Quirós", "Yachida", "Stevens", "Labbé"]:
        for m in multiple_test_options:
            t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
            t_test_res = t_test_res.sort_values(by=["P-adjust"])
            DA_max = len(t_test_res[t_test_res["P-adjust"] < 0.10]["Metabolite"].tolist())
            size_range = [i for i in range(0, 105, 5)]
            for s in size_range:
                # DA_metabolites = t_test_res[t_test_res["P-adjust"] < s]["Metabolite"].tolist()
                DA_metabolites = t_test_res.head(int(DA_max * (s / 100)))["Metabolite"].tolist()
                ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                res_list.append([d, s, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    print(res_df_bonferroni)
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    with plt.style.context('seaborn'):
        fig, ax1 = plt.subplots(1, 1)
        plt.style.use("seaborn")
        cols = sns.color_palette("deep", 8)

        for num, i in enumerate(["Quirós", "Yachida", "Stevens", "Labbé"]):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_numpy().astype('int'),
                         res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" Bonferroni", linestyle='-', color=cols[num])
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_numpy().astype('int'),
                         res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" FDR BH", linestyle='--', color=cols[num])
        # plt.title("Reactome", fontsize=14)

        ax1.legend()
        ax1.axvline(100, color='black')
        ax1.set_ylabel("Number of pathways significant at P < 0.1")
        ax1.set_xlabel("Top N% of Q-values (100% is equivalent to Q-value < 0.10)")
        plt.savefig("../Figures/KEGG_pq_value_0_10_cutoff.png", dpi=300)
        plt.show()


def vary_p_val_cutoff():
    proportion_of_bg = [i for i in range(0, 85, 5)]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []

    # Import files for Fuhrer data
    n_zscore = pd.read_csv("../Fuhrer/mod_zscore_neg_CW.csv", index_col=0)
    p_zscore = pd.read_csv("../Fuhrer/mod_zscore_pos_CW.csv", index_col=0)
    # remove unannotated
    n_zscore = n_zscore[n_zscore.index.notnull()]
    p_zscore = p_zscore[p_zscore.index.notnull()]
    with open('zamboni_pos_annotation_dict.pickle', 'rb') as handle:
        annotations_pos = pickle.load(handle)
    with open('zamboni_neg_annotation_dict.pickle', 'rb') as handle:
        annotations_neg = pickle.load(handle)

    for d in datasets.keys():
        print(d)
        for p in proportion_of_bg:
            print(p)
            if d.startswith("Fuhrer"):
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

        for num, i in enumerate(["Quirós", "Yachida", "Stevens", "Labbé"]):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_list(),
                         res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" Bonferroni", linestyle='-', color=cols[num])
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_list(),
                         res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i+" FDR BH", linestyle='--', color=cols[num])
        for num, i in enumerate(["Fuhrer (dcuS)", "Fuhrer (yfgM)"]):
            ax1.plot(res_df_Zamboni[res_df_Zamboni["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_Zamboni[res_df_Zamboni["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i,
                     linestyle='-', color=cols[4+num])

        ax1.legend()
        ax1.set_ylabel("Number of pathways significant at P < 0.1")
        ax1.set_xlabel("Top N% of differentially abundant metabolites (ranked by Q-value)")
        # plt.savefig("../Figures/vary_DAM_size_percentage.png", dpi=300)
        # plt.show()

# vary_dam_size()

def vary_pval():
    cutoffs = [0.001, 0.005, 0.01, 0.05, 0.1]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in ["Quirós", "Labbé", "Yachida", "Stevens"]:
        for c in cutoffs:
            for m in multiple_test_options:
                t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                DA_metabolites = t_test_res[t_test_res["P-adjust"] < c]["Metabolite"].tolist()
                ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                res_list.append([d, c, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])

    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])
    res_df.to_csv("p-value_cutoffs.csv")
    res_df = pd.read_csv("p-value_cutoffs.csv")
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    res_df_FDR_BH = res_df_FDR_BH.drop("Multiple_correction_method", axis=1)
    res_df_FDR_BH = res_df_FDR_BH.pivot(index='Dataset', columns=['Cutoff_P'], values='n_p_less_01')
    print(res_df_FDR_BH)
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_bonferroni = res_df_bonferroni.drop("Multiple_correction_method", axis=1)
    res_df_bonferroni = res_df_bonferroni.pivot(index='Dataset', columns=['Cutoff_P'],
                                        values='n_p_less_01')

    #
    with plt.style.context('seaborn'):
        sns.set_palette("Blues_r")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

        res_df_FDR_BH.plot.bar(ax=ax1)
        ax1.set_title('Benjamini-Hochberg FDR')
        ax1.set_ylabel("Number of pathways significant at P < 0.1")
        # plt.ylabel("Number of significant pathway at P < 0.1")
        # plt.xlabel("Dataset")
        ax1.legend(title="Q-value threshold")

        res_df_bonferroni.plot.bar(ax=ax2)
        ax2.set_title('Bonferroni')
        ax2.legend(title="Q-value threshold")
        # plt.title('Bonferroni')
        plt.tight_layout()
        plt.savefig("../Figures/vary_pvalue_cutoff.png", dpi=300)
        plt.show()


vary_pval()

