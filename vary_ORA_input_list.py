import utils
import process_datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm

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
datasets = {"Labbé": [DEM_brown, background_brown, KEGG_mouse_pathways, all_KEGG_mouse_bg, mat_brown,
                      [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
"Yachida": [DEM_yamada, background_yamada, KEGG_human_pathways, all_KEGG_human_bg, mat_yamada,
            [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
"Stevens": [DEM_stevens, background_stevens, KEGG_human_pathways, all_KEGG_human_bg,
            mat_stevens, [i for i in range(0, 40, 5)], [i for i in range(0, 35, 5)]],
"Quirós": [DEM_auwerx, background_auwerx, KEGG_human_pathways, all_KEGG_human_bg, mat_auwerx,
           [i for i in range(0, 14, 1)], [i for i in range(0, 12, 1)]],
"Fuhrer (yfgM)": [DEM_yfgM, background_yfgM, KEGG_eco_pathways, all_KEGG_eco_bg, mat_yfgM,
                  [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]],
"Fuhrer (dcuS)": [DEM_dcuS, background_dcuS, KEGG_eco_pathways, all_KEGG_eco_bg, mat_dcuS,
                  [i for i in range(0, 7, 1)], [i for i in range(0, 6, 1)]]}

print("Data import complete")


def vary_pval():
    cutoffs = [0.001, 0.005, 0.01, 0.05, 0.1]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in datasets:
        for c in cutoffs:
            for m in multiple_test_options:
                if d.startswith("Fuhrer"):
                    p_vals_dict = {}
                    for x in datasets[d][4].T.itertuples():
                        zscore = x[1]
                        pval = stats.norm.cdf(zscore)
                        p_vals_dict[x[0]] = pval
                    pvals_df = pd.DataFrame(p_vals_dict, columns=["Metabolite", "P-value"])
                    padj = sm.stats.multipletests(pvals_df["P-value"], 0.05, method=m)
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, c, m, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())])
                else:
                    t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                    DA_metabolites = t_test_res[t_test_res["P-adjust"] <= c]["Metabolite"].tolist()
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, c, m, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())])

    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])
    res_df.to_csv("p-value_cutoffs.csv")
    res_df = pd.read_csv("p-value_cutoffs.csv")
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    res_df_FDR_BH = res_df_FDR_BH.drop("Multiple_correction_method", axis=1)
    res_df_FDR_BH = res_df_FDR_BH.pivot(index='Dataset', columns=['Cutoff_P'], values='n_p_less_01')
    res_df_FDR_BH = res_df_FDR_BH.reindex(list(datasets.keys()))
    print(res_df_FDR_BH)
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_bonferroni = res_df_bonferroni.drop("Multiple_correction_method", axis=1)
    res_df_bonferroni = res_df_bonferroni.pivot(index='Dataset', columns=['Cutoff_P'],
                                                values='n_p_less_01')
    res_df_bonferroni = res_df_bonferroni.reindex(list(datasets.keys()))

    with plt.style.context('seaborn'):
        sns.set_palette("Blues_r")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

        res_df_FDR_BH.plot.bar(ax=ax1)
        ax1.set_title('Benjamini-Hochberg FDR', fontsize=14)
        ax1.set_ylabel("Number of pathways significant at P ≤ 0.1", fontsize=13)
        # plt.ylabel("Number of significant pathway at P < 0.1")
        # plt.xlabel("Dataset")
        ax1.legend(title="Q-value threshold (≤)", fontsize=11)

        res_df_bonferroni.plot.bar(ax=ax2)
        ax2.set_title('Bonferroni', fontsize=14)
        ax2.legend(title="Q-value threshold (≤)", fontsize=11)
        # plt.title('Bonferroni')
        plt.tight_layout()
        plt.savefig("../Figures/vary_pvalue_cutoff.png", dpi=600)
        plt.show()


# vary_pval()

def vary_dam_size():
    proportion_of_bg = [i for i in range(0, 101, 1)]
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in datasets.keys():
        print(d)
        for p in proportion_of_bg:
            for m in multiple_test_options:
                if d.startswith("Fuhrer"):
                    p_vals_dict = {}
                    for x in datasets[d][4].T.itertuples():
                        zscore = x[1]
                        pval = stats.norm.cdf(zscore)
                        p_vals_dict[x[0]] = pval
                    pvals_df = pd.DataFrame.from_dict(p_vals_dict, orient='index', columns=["P-value"])
                    padj = sm.stats.multipletests(pvals_df["P-value"], 0.05, method=m)[1]
                    pvals_df["P-adjust"] = padj
                    pvals_df["ranked_padj"] = stats.rankdata(pvals_df["P-adjust"], method='min')
                    pvals_df = pvals_df.sort_values(by="ranked_padj")
                    select_up_to = int((len(pvals_df["P-adjust"]) - 1) * (p / 100))
                    selected_rank = pvals_df.iloc[select_up_to, 2]
                    DA_metabolites = pvals_df[pvals_df["ranked_padj"] <= selected_rank].index.tolist()
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, p, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
                else:
                    t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                    t_test_res["ranked_padj"] = stats.rankdata(t_test_res["P-adjust"], method='min')
                    t_test_res = t_test_res.sort_values(by="ranked_padj")
                    select_up_to = int((len(t_test_res["Metabolite"])-1) * (p/100))
                    selected_rank = t_test_res.iloc[select_up_to, 3]
                    DA_metabolites = t_test_res[t_test_res["ranked_padj"] <= selected_rank]["Metabolite"].tolist()
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, p, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])

    # add annotations
    cutoffs = [0.05]
    annotations = []
    for d in datasets.keys():
        for c in cutoffs:
            for m in multiple_test_options:
                if d.startswith("Fuhrer"):
                    p_vals_dict = {}
                    for x in datasets[d][4].T.itertuples():
                        zscore = x[1]
                        pval = stats.norm.cdf(zscore)
                        p_vals_dict[x[0]] = pval
                    pvals_df = pd.DataFrame.from_dict(p_vals_dict, orient='index', columns=["P-value"])
                    padj = sm.stats.multipletests(pvals_df["P-value"], 0.05, method=m)[1]
                    pvals_df["P-adjust"] = padj
                    DA_metabolites = pvals_df[pvals_df["P-adjust"] <= c].index.tolist()
                    percentage = (len(DA_metabolites) / len(datasets[d][1])) * 100
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_df.loc[len(res_df), :] = [d, percentage, m,
                                                  len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())]

                    annotations.append(
                        [c, percentage, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist()), d, m])
                else:
                    t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                    DA_metabolites = t_test_res[t_test_res["P-adjust"] <= c]["Metabolite"].tolist()
                    percentage = (len(DA_metabolites) / len(t_test_res["Metabolite"])) * 100
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_df.loc[len(res_df), :] = [d, percentage, m,
                                                  len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())]

                    annotations.append(
                        [c, percentage, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist()), d, m])

    cutoffs_df = pd.DataFrame(annotations)
    cutoffs_bonferroni = cutoffs_df[cutoffs_df[4] == "bonferroni"]
    cutoffs_FDR_BH = cutoffs_df[cutoffs_df[4] == "fdr_bh"]

    res_df = res_df.sort_values(by='Cutoff_P')
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    print(res_df)
    with plt.style.context('seaborn-darkgrid'):
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8), dpi=600)
        cols = sns.color_palette("muted", 8)
        for num, i in enumerate(datasets.keys()):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o',
                     label=i + " Bonferroni", linestyle='-', color=cols[num], markersize=3)
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i + " FDR BH",
                     linestyle='dotted', color=cols[num], markersize=3)
        ax1.legend(fontsize=11)
        ax1.set_ylabel("Number of pathways significant at P ≤ 0.1", fontsize=13)
        ax1.set_xlabel("Top N% of metabolites (ranked by adjusted P or Q-value)", fontsize=13)
        for num, row in enumerate(cutoffs_bonferroni.itertuples()):
            plt.text(row[2], row[3], row[1], bbox=dict(boxstyle="round",
                                                       ec=cols[num],
                                                       fc='white'))
        for num, row in enumerate(cutoffs_FDR_BH.itertuples()):
            plt.text(row[2], row[3], row[1], bbox=dict(boxstyle="round",
                                                       ec=cols[num],
                                                       fc='white',
                                                       linestyle="dotted"))
        plt.tight_layout()
        plt.savefig("vary_input_metabolites_all_new.png", dpi=600)
        plt.show()


vary_dam_size()

def vary_dam_size_single(dataset):
    multiple_test_options = ["bonferroni", "fdr_bh"]
    res_list = []
    for d in [dataset]:
        print(d)
        for m in multiple_test_options:
            if d.startswith("Fuhrer"):
                p_vals_dict = {}
                # for x in datasets[d][4].T.itertuples():
                #     zscore = x[1]
                #     pval = stats.norm.cdf(zscore)
                #     p_vals_dict[x[0]] = pval
                # pvals_df = pd.DataFrame.from_dict(p_vals_dict, orient='index', columns=["P-value"])
                # padj = sm.stats.multipletests(pvals_df["P-value"], 0.05, method=m)[1]
                # pvals_df["P-adjust"] = padj
                # padj_sorted = pvals_df.sort_values(by=["P-adjust"])
                # print(DA_metabolites)
                # ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                # res_list.append([d, i, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])
            else:
                t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                t_test_res["ranked_padj"] = stats.rankdata(t_test_res["P-adjust"], method='min')
                t_test_res = t_test_res.sort_values(by="ranked_padj")
                # t_test_res.to_csv("Labbe_ttest_res" + m + ".csv")

                max_rank = t_test_res["ranked_padj"].max()
                print(max_rank)
                for i in range(1, len(t_test_res["Metabolite"])+1):
                    metabolite_rank = t_test_res.iloc[i-1, 3]
                    print("Metabolite:", i, "rank:", metabolite_rank)
                    DA_metabolites = t_test_res[t_test_res["ranked_padj"] <= metabolite_rank]["Metabolite"].tolist()
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_list.append([d, i, m, len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())])

    res_df = pd.DataFrame(res_list, columns=["Dataset", "Cutoff_P", "Multiple_correction_method", "n_p_less_01"])

    # add annotations
    cutoffs = [0.005, 0.05, 0.1]
    annotations = []
    for d in [dataset]:
        for c in cutoffs:
            for m in multiple_test_options:
                if d.startswith("Fuhrer"):
                    p_vals_dict = {}
                    # for x in datasets[d][4].T.itertuples():
                    #     zscore = x[1]
                    #     pval = stats.norm.cdf(zscore)
                    #     p_vals_dict[x[0]] = pval
                    # pvals_df = pd.DataFrame.from_dict(p_vals_dict, orient='index', columns=["P-value"])
                    # padj = sm.stats.multipletests(pvals_df["P-value"], 0.05, method=m)[1]
                    # pvals_df["P-adjust"] = padj
                    # DA_metabolites = pvals_df[pvals_df["P-adjust"] <= c].index.tolist()
                    # percentage = (len(DA_metabolites) / len(datasets[d][1])) * 100
                    # ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    # res_df.loc[len(res_df), :] = [d, percentage, m,
                    #                               len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())]
                    #
                    # annotations.append(
                    #     [c, percentage, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist()), d, m])
                else:
                    t_test_res = utils.t_tests(datasets[d][4].iloc[:, :-1], datasets[d][4]["Group"], m)
                    DA_metabolites = t_test_res[t_test_res["P-adjust"] <= c]["Metabolite"].tolist()
                    n_DA = len(DA_metabolites)
                    ora_res = utils.over_representation_analysis(DA_metabolites, datasets[d][1], datasets[d][2])
                    res_df.loc[len(res_df), :] = [d, n_DA, m,
                                                  len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist())]

                    annotations.append(
                        [c, n_DA, len(ora_res[ora_res["P-value"] <= 0.1]["P-value"].tolist()), d, m])

    cutoffs_df = pd.DataFrame(annotations)
    cutoffs_bonferroni = cutoffs_df[cutoffs_df[4] == "bonferroni"]
    cutoffs_FDR_BH = cutoffs_df[cutoffs_df[4] == "fdr_bh"]

    res_df = res_df.sort_values(by='Cutoff_P')
    res_df_bonferroni = res_df[res_df["Multiple_correction_method"] == "bonferroni"]
    res_df_FDR_BH = res_df[res_df["Multiple_correction_method"] == "fdr_bh"]
    print(res_df)
    with plt.style.context('seaborn-darkgrid'):
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        cols = sns.color_palette("muted", 8)
        for num, i in enumerate(datasets.keys()):
            ax1.plot(res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_bonferroni[res_df_bonferroni["Dataset"] == i]['n_p_less_01'].tolist(), 'o',
                     label=i + " Bonferroni", linestyle='-', color=cols[num], markersize=3)
            ax1.plot(res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['Cutoff_P'].to_list(),
                     res_df_FDR_BH[res_df_FDR_BH["Dataset"] == i]['n_p_less_01'].tolist(), 'o', label=i + " FDR BH",
                     linestyle='dotted', color=cols[0], markersize=3)
        ax1.legend(fontsize=11)
        ax1.set_ylabel("Number of pathways significant at P ≤ 0.1", fontsize=13)
        ax1.set_xlabel("Top N metabolites (ranked by adjusted P or Q-value)", fontsize=13)
        for num, row in enumerate(cutoffs_bonferroni.itertuples()):
            plt.text(row[2], row[3], row[1], bbox=dict(boxstyle="round",
                                                       ec=cols[0],
                                                       fc='white'), ha="right", va="baseline")
        for num, row in enumerate(cutoffs_FDR_BH.itertuples()):
            plt.text(row[2], row[3], row[1], bbox=dict(boxstyle="round",
                                                       ec=cols[0],
                                                       fc='white',
                                                       linestyle="dotted"), ha="right", va="baseline")
        plt.tight_layout()
        plt.savefig("../Figures/vary_input_metabolites_Labbe.png", dpi=600)
        plt.show()


# vary_dam_size_single("Labbé")