# Tools for pathway analysis in metabolomics

# TODO: make requirements.txt

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm

pd.options.mode.chained_assignment = None  # default='warn'

def data_processing(raw_matrix, firstrow, firstcol):
    '''
    Filtering low abundance metabolites, data cleaning and imputation using minimum value.
    :param raw_matrix: raw abundance matrix with n-samples and m-metabolites
    :param firstrow: First row containing values (integer)
    :return: imputed, log-transformed and standardised matrix
    '''

    # Remove commas and convert to numeric
    processed_matrix = raw_matrix.iloc[firstrow:, firstcol:]
    processed_matrix = processed_matrix.replace(',', '', regex=True)
    processed_matrix = processed_matrix.apply(pd.to_numeric)
    processed_matrix.replace(0, np.nan, inplace=True)

    processed_matrix = processed_matrix.loc[:, processed_matrix.isnull().mean() < 0.9]
    # Remove metabolites not present in > 90% of samples
    # by indexing df by rows with metabolites present in more than 90% of samples

    # Missing value imputation using minimum value/2
    imputed_matrix = processed_matrix.replace(np.nan, processed_matrix.min(axis=0)/2)
    # Log2 transformation
    log2_matrix = np.log2(imputed_matrix)
    # Standardisation by mean centering and scaling to unit variance
    standarised_mat = StandardScaler().fit_transform(log2_matrix)
    log2_matrix.loc[:, :] = standarised_mat
    return log2_matrix

def plot_PCA(matrix, metadata, title, n_comp=100):
    '''
    PCA plot
    :param matrix: processed data matrix
    :param metadata: pandas series with group names
    :param title: plot title
    :param n_comp: number of PCA components
    :return: PCA plot
    '''

    pca = PCA(n_components=n_comp)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    group = metadata
    uniq_sample = np.unique(group)
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

def linear_regression(matrix, metadatadict):
    # Add new columns based on metadata from dict
    matrix['PMH_status'] = matrix.index.map(lambda x: metadatadict[x][0])
    matrix['Target'] = pd.factorize(matrix['PMH_status'])[0]
    coefs = []
    pvals = []
    metabolites = matrix.columns.tolist()[:-2]
    for metabolite in metabolites:
        X = matrix[str(metabolite)]
        y = matrix['Target']
        model = sm.OLS(y, X.astype(float)).fit()
        coefs.append(model.params[0])
        pvals.append(model.pvalues[0])
    padj = sm.stats.multipletests(pvals, 0.05, method="bonferroni")
    print("Corrected alpha:", padj[3])
    results = pd.DataFrame(zip(metabolites, coefs, pvals, padj[1]), columns=["Metabolite", "Fold-change", "P-value", "P-adjust"])
    return results

def t_tests(matrix, classes, multiple_correction_method):

    matrix['Target'] = pd.factorize(classes)[0]
    metabolites = matrix.columns.tolist()[:-1]

    pvalues = []
    for metabolite in metabolites:
        group1 = matrix[matrix['Target'] == 0][metabolite].tolist()
        group2 = matrix[matrix['Target'] == 1][metabolite].tolist()
        stat, pval = stats.ttest_ind(group1, group2)
        pvalues.append(pval)
    padj = sm.stats.multipletests(pvalues, 0.05, method=multiple_correction_method)
    results = pd.DataFrame(zip(metabolites, pvalues, padj[1]),
                           columns=["Metabolite", "P-value", "P-adjust"])
    return results

def over_representation_analysis(DEM_list, background_list, pathways_df):
    # analyse each pathway
    KEGG_pathways = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])

    pathways = KEGG_pathways.index.tolist()
    pathway_names = KEGG_pathways["Pathway_name"].tolist()
    pathway_dict = dict(zip(pathways, pathway_names))
    KEGG_pathways.drop('Pathway_name', axis=1, inplace=True)

    pathways_with_compounds = []
    pathway_names_with_compounds = []
    pvalues = []
    pathway_ratio = []
    for pathway in pathways:
        pathway_compounds = KEGG_pathways.loc[pathway, :].tolist()
        pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
        if not pathway_compounds or len(pathway_compounds) < 3:
            continue
        else:
            DEM_in_pathway = len(set(DEM_list) & set(pathway_compounds))
            DEM_not_in_pathway = len(np.setdiff1d(DEM_list, pathway_compounds))
            compound_in_pathway_not_DEM = len(set(background_list) & set(pathway_compounds))
            compound_not_in_pathway_not_DEM = len(np.setdiff1d(background_list, pathway_compounds))

            if (DEM_in_pathway and compound_in_pathway_not_DEM) == 0:
                continue
            else:
                # Create 2 by 2 contingency table
                pathway_ratio.append(str(DEM_in_pathway) + "/" + str(len(pathway_compounds)))
                pathways_with_compounds.append(pathway)
                pathway_names_with_compounds.append(pathway_dict[pathway])
                contingency_table = np.array([[DEM_in_pathway, compound_in_pathway_not_DEM],
                                              [DEM_not_in_pathway, compound_not_in_pathway_not_DEM]])
                # Run right tailed Fisher's exact test
                oddsratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
                pvalues.append(pvalue)
    # padj = sm.stats.multipletests(pvalues, 0.05, method="fdr_bh")
    try:
        padj = sm.stats.multipletests(pvalues, 0.05, method="fdr_bh")
        results = pd.DataFrame(zip(pathways_with_compounds, pathway_names_with_compounds, pathway_ratio, pvalues, padj[1]),
                               columns=["Pathway_ID", "Pathway_name", "Hits", "P-value", "P-adjust"])
    except ZeroDivisionError:
        padj = [1] * len(pvalues)
        results = pd.DataFrame(zip(pathways_with_compounds, pathway_names_with_compounds, pathway_ratio, pvalues, padj),
                           columns=["Pathway_ID", "Pathway_name", "Hits", "P-value", "P-adjust"])
    return results

    # results = pd.DataFrame(zip(pathways_with_compounds, pathway_names_with_compounds, pvalues, padj[1]),
    #                        columns=["Pathway_ID", "Pathway_name", "P-value", "P-adjust"])
    # return results

def reduce_background_list_ora(background_list, percentage, DEM_list, pathways_df):
    '''
    Reduces size of background list by random removal of compounds
    :param background_list: background list of compound names/IDs
    :param percentage: percentage reduction of list desired
    :return: reduced background list
    '''
    list_size = int(len(background_list)*((100-percentage)/100))
    print(list_size)

    p_vals = []
    q_vals = []
    for i in range(0, 100):
        bg_list_reduced = np.random.choice(background_list, list_size, replace=False)
        ora_res = over_representation_analysis(DEM_list, bg_list_reduced, pathways_df)
        p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
        q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    sd_p_signficant_paths = np.std(p_vals)
    sd_q_signficant_paths = np.std(q_vals)
    return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths]

