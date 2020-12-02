import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils

KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)

def stevens_data():
    md_raw = pd.read_csv("../Stevens/MTBLS136_compressed_files/s_MTBLS136.txt", sep="\t")
    metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'],
                             md_raw['Factor Value[AgeAtBloodDraw]'],
                             ['Over 75' if val not in ['<=55', '56-60', '61-65', '66-70', '71-75'] else 'Under 75' for val in md_raw['Factor Value[AgeAtBloodDraw]']]))
    metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))
    sample_status_dict = dict(zip(md_raw['Sample Name'].values, md_raw['Factor Value[CurrentPMH]']))

    replicate_samples = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', 'E+P']]
    nonusers = [k for k, v in metadata_dict.items() if v[0] not in [np.nan, 'E-only', 'E+P']]
    estrogen_only = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', np.nan, 'E+P']]
    estrogen_progesterone = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', np.nan]]

    # Get abundance matrix, transpose to n-samples by m-metabolites
    mat = pd.read_csv("Stevens_matrix_named_compounds_only.csv", index_col=0)
    mat_nonusers_estrogen = mat.drop((replicate_samples + estrogen_progesterone), axis=1)
    stevens_matrix_proc = utils.data_processing(mat_nonusers_estrogen.T, 8, 0)
    stevens_matrix_proc["Group"] = stevens_matrix_proc.index.map(sample_status_dict)
    ttest_res = utils.t_tests(stevens_matrix_proc.iloc[:,:-1], stevens_matrix_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()

    # Add back the metadata columns
    mat.columns = mat.columns.str.replace(' ', '')
    metadata_cols = ['KEGG', 'SampleHMDB_ID']
    stevens_matrix_proc_annotated = stevens_matrix_proc.T.join(mat[metadata_cols])
    background_list = stevens_matrix_proc.columns.tolist()
    DEM_KEGG_id = mat[mat.index.isin(DEM)]['KEGG'].tolist()
    DEM_KEGG_id = [i for i in DEM_KEGG_id if str(i) != 'nan']
    stevens_background_list = stevens_matrix_proc_annotated['KEGG'].dropna().tolist()

    # Calculate the classes of metabolites not mapped to KEGG
    # classes = pd.read_csv("Stevens_metabolite_classes.csv", index_col=1)
    #
    # mat_classes = stevens_matrix_proc_annotated.merge(classes, left_on="SampleHMDB_ID", right_on="HMDB")
    #
    # metabolites_not_in_KEGG = mat_classes.loc[mat_classes['KEGG'].isna(), 'SampleHMDB_ID'].tolist()
    # metabolites_not_in_KEGG_classes = mat_classes[mat_classes['SampleHMDB_ID'].isin(metabolites_not_in_KEGG)][
    #     'HMDB_class']
    # sns.countplot(y=metabolites_not_in_KEGG_classes, data=mat_classes)
    # plt.title("Classes of metabolites not mapped to KEGG pathways")
    # plt.savefig("Stevens_not_in_KEGG_all.png")
    # plt.show()
    # DE_metabolites_not_in_KEGG = mat_classes.loc[mat_classes['KEGG'].isin(DEM_KEGG_id), 'SampleHMDB_ID'].tolist()
    # print(len(DE_metabolites_not_in_KEGG))
    # DE_metabolites_not_in_KEGG_classes = mat_classes[mat_classes['SampleHMDB_ID'].isin(DE_metabolites_not_in_KEGG)]['HMDB_class']
    # sns.countplot(y=DE_metabolites_not_in_KEGG_classes, data=mat_classes)
    # plt.title("Classes of DA metabolites not mapped to KEGG pathways")
    # # plt.savefig("Stevens_not_in_KEGG_all.png")
    # plt.show()

    return DEM_KEGG_id, stevens_background_list

DEM, background = stevens_data()
print(DEM)
print(background)



def plot_PCA(matrix, metadata, title):
    '''
    PCA plot
    :param matrix:
    :return: PCA plot
    '''

    pca = PCA(n_components=100)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    group = np.array([metadata[i][3] for i in samples])
    uniq_sample = np.unique(group)
    print(uniq_sample)
    cmap = ['tab:green', 'tab:orange', 'tab:blue', 'red', 'cyan', 'magenta', 'brown', 'yellow']
    cdict = {samp: cmap[num] for num, samp in enumerate(uniq_sample)}

    plt.style.use("ggplot")
    fig, ax = plt.subplots()
    for g in np.unique(group):
        ix = np.where(group == g)
        ax.scatter(scatter_x[ix], scatter_y[ix], c=cdict[g], label=g, s=15)
    ax.legend()

    plt.title("PCA for " + title)
    plt.xlabel("Component 1:  " + str(round(pca.explained_variance_ratio_[0]*100, 2)) + "%")
    plt.ylabel("Component 2: " + str(round(pca.explained_variance_ratio_[1]*100, 2)) + "%")
    plt.savefig(title + ".png")
    plt.show()

def linear_regression(matrix):
    # Add new columns based on metadata from dict
    matrix['PMH_status'] = matrix.index.map(lambda x: metadata_dict[x][0])
    matrix['Target'] = pd.factorize(matrix['PMH_status'])[0]
    coefs = []
    pvals = []
    metabolites = matrix.columns.tolist()[:-2]

    for metabolite in metabolites:
        X = matrix[str(metabolite)]
        y = matrix['Target']
        model = sm.OLS(y, X.astype(float)).fit()
        coefs.append(round(model.params[0], 4))
        pvals.append(round(model.pvalues[0], 4))
    padj = sm.stats.multipletests(pvals, 0.05, method="bonferroni")
    print("Corrected alpha:", padj[3])
    results = pd.DataFrame(zip(metabolites, coefs, pvals, padj[1]), columns=["Metabolite", "Fold-change", "P-value", "P-adjust"])
    return results


def hierarchical_clustering(matrix):
    h = sns.clustermap(matrix, metric="correlation", cmap="vlag")
    plt.show()
    return h

def over_representation_analysis(DEM_list, background_list, pathways_df):
    # analyse each pathway
    pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:], inplace=True)
    pathways = pathways_df.index.tolist()
    pathway_names = pathways_df["Pathway_name"].tolist()
    pathways_df.drop('Pathway_name', axis=1, inplace=True)

    pvalues = []
    for pathway in pathways:
        pathway_compounds = pathways_df.loc[pathway, :].tolist()
        pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
        if not pathway_compounds:
            break
        else:
            # Create 2 by 2 contingency table
            DEM_in_pathway = len(set(DEM_list) & set(pathway_compounds))
            DEM_not_in_pathway = len(np.setdiff1d(DEM_list, pathway_compounds))
            compound_in_pathway_not_DEM = len(set(background_list) & set(pathway_compounds))
            compound_not_in_pathway_not_DEM = len(np.setdiff1d(background_list, pathway_compounds))
            contingency_table = np.array([[DEM_in_pathway, compound_in_pathway_not_DEM],
                                          [DEM_not_in_pathway, compound_not_in_pathway_not_DEM]])
            # Run right tailed Fisher's exact test
            oddsratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
            pvalues.append(pvalue)

    padj = sm.stats.multipletests(pvalues, 0.05, method="fdr_bh")
    results = pd.DataFrame(zip(pathways, pathway_names, pvalues, padj[1]),
                           columns=["Pathway_ID", "Pathway_name", "P-value", "P-adjust"])
    return results

# Calculate the classes of metabolites not mapped to KEGG
# classes = pd.read_csv("Stevens_metabolite_classes.csv", index_col=1)
#
# mat_classes = stevens_matrix_proc_annotated.merge(classes, left_on="SampleHMDB_ID", right_on="HMDB")
#
# metabolites_not_in_KEGG = mat_classes.loc[mat_classes['KEGG'].isna(), 'SampleHMDB_ID'].tolist()
# metabolites_not_in_KEGG_classes = mat_classes[mat_classes['SampleHMDB_ID'].isin(metabolites_not_in_KEGG)]['HMDB_class']
# sns.countplot(y=metabolites_not_in_KEGG_classes, data=mat_classes)
# plt.title("Classes of metabolites not mapped to KEGG pathways")
# plt.savefig("Stevens_not_in_KEGG_all.png")
# plt.show()
# DE_metabolites_not_in_KEGG = mat_classes.loc[mat_classes['KEGG'].isin(DEM_KEGG_id), 'SampleHMDB_ID'].tolist()
# print(len(DE_metabolites_not_in_KEGG))
# DE_metabolites_not_in_KEGG_classes = mat_classes[mat_classes['SampleHMDB_ID'].isin(DE_metabolites_not_in_KEGG)]['HMDB_class']
# sns.countplot(y=DE_metabolites_not_in_KEGG_classes, data=mat_classes)
# plt.title("Classes of DA metabolites not mapped to KEGG pathways")
# # plt.savefig("Stevens_not_in_KEGG_all.png")
# plt.show()


ORA_res = utils.over_representation_analysis(DEM, background, KEGG_pathways)
print(len(ORA_res[ORA_res["P-value"] < 0.1]["P-value"].tolist()))
print(len(ORA_res[ORA_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
# ORA_res.to_csv("../Stevens/ORA_Stevens.csv")

all_KEGG_compounds = list(set([x for x in KEGG_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
print(len(all_KEGG_compounds))

ora_whole_KEGG_bg = utils.over_representation_analysis(DEM, all_KEGG_compounds, KEGG_pathways)
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-value"] < 0.1]["P-value"].tolist()))
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-adjust"] < 0.1]["P-adjust"].tolist()))
#
# print(DEM_KEGG_id)
# differential_metabolites.to_csv("differential_metabolites_nonuser_v_estrogen.csv")
# name_map = pd.read_csv("../Stevens/name_map_metaboanalyst.csv", dtype=str)
# Reactome_pathways = pd.read_csv("Reactome_pathway_set.csv", dtype=str, index_col=0)
# Reactome_human_ids = [i for i in Reactome_pathways.index if i.startswith("R-HSA")]
# Reactome_human = Reactome_pathways[Reactome_pathways.index.isin(Reactome_human_ids)]
# background_list_chEBI = name_map[name_map["Query"].isin(background_list)]['ChEBI'].dropna().tolist()

# DA_metabolites_chEBI = name_map[name_map["Query"].isin(DEM)]['ChEBI'].dropna().tolist()
# print(len(background_list_chEBI), len(DA_metabolites_chEBI))
# ORA_reactome = utils.over_representation_analysis(DA_metabolites_chEBI, background_list_chEBI, Reactome_human)
# print(len(ORA_reactome[ORA_reactome["P-adjust"] < 0.1]["P-adjust"].tolist()))
# ORA_reactome.to_csv("Stevens_ORA_Reactome.csv")
# mat_nonusers_e_p = stevens_matrix_proc.drop((replicate_samples + estrogen_only), axis=0)
# # mat_estrogen_vs_e_p = stevens_matrix_proc.drop((replicate_samples + nonusers), axis=1)
# mat_all = stevens_matrix_proc.drop(replicate_samples, axis=0)
#
# plot_PCA(mat_all, metadata_dict, "all groups")
# plot_PCA(mat_nonusers_e_p, metadata_dict, "nonusers vs. estrogen+progestin")

