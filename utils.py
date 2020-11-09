# Tools for pathway analysis in metabolomics

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm

def data_processing(raw_matrix, firstrow):
    '''
    Filtering low abundance metabolites, data cleaning and imputation using minimum value.
    :param raw_matrix: raw abundance matrix with n-samples and m-metabolites
    :param firstrow: First row containing values (integer)
    :return: imputed, log-transformed and standardised matrix
    '''

    processed_matrix = raw_matrix.loc[:, raw_matrix.isnull().mean() < 0.9]
    # Remove metabolites not present in > 90% of samples
    # by indexing df by rows with metabolites present in more than 90% of samples
    processed_matrix = processed_matrix.iloc[firstrow:, :]
    # Remove commas and convert to numeric
    # processed_matrix = processed_matrix.replace(np.nan, 0)
    processed_matrix = processed_matrix.replace(',','', regex=True)
    processed_matrix = processed_matrix.apply(pd.to_numeric)
    # Missing value imputation using minimum value/2
    imputed_matrix = processed_matrix.replace(np.nan, processed_matrix.min(axis=0)/2)
    # Log2 transformation
    log2_matrix = np.log2(imputed_matrix)
    # Standardisation by mean centering and scaling to unit variance
    standarised_mat = StandardScaler().fit_transform(log2_matrix)
    log2_matrix.loc[:, :] = standarised_mat
    return log2_matrix

def plot_PCA(matrix, metadata, title, labelpos):
    '''
    PCA plot
    :param matrix: processed data matrix
    :param labelpos: position in metadata list of sample labels of interest
    :return: PCA plot
    '''

    pca = PCA(n_components=100)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    group = np.array([metadata[i][labelpos] for i in samples])
    uniq_sample = np.unique(group)
    cmap = ['tab:green', 'tab:orange', 'tab:blue']
    cdict = {samp: cmap[num] for num, samp in enumerate(uniq_sample)}

    fig, ax = plt.subplots()
    for g in np.unique(group):
        ix = np.where(group == g)
        ax.scatter(scatter_x[ix], scatter_y[ix], c=cdict[g], label=g, s=10)
    ax.legend()
    plt.style.use("ggplot")
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

def over_representation_analysis(DEM_list, background_list, pathways_df):
    # analyse each pathway
    pathways = pathways_df.index.tolist()
    pvalues = []
    for pathway in pathways:
        pathway_compounds = pathways_df.loc[pathway, :].tolist()
        pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
        if not pathway_compounds:
            print("Pathway contains no compounds ")
            pvalues.append(np.nan)
        else:
            # Create 2 by 2 contingency table
            DEM_in_pathway = len(set(DEM_list) & set(pathway_compounds))
            DEM_not_in_pathway = len(np.setdiff1d(pathway_compounds, DEM_list))
            compound_in_pathway_not_DEM = len(set(background_list) & set(pathway_compounds))
            compound_not_in_pathway_not_DEM = len(np.setdiff1d(background_list, pathway_compounds))
            contingency_table = np.array([[compound_in_pathway_not_DEM, DEM_in_pathway], [compound_not_in_pathway_not_DEM, DEM_not_in_pathway]])
            # Run right tailed Fisher's exact test
            oddsratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
            pvalues.append(pvalue)
    padj = sm.stats.multipletests(pvalues, 0.05, method="bonferroni")
    results = pd.DataFrame(zip(pathways, pvalues, padj[1]),
                           columns=["Pathway_ID", "P-value", "P-adjust"])
    return results