# Tools for pathway analysis in metabolomics

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import seaborn as sns

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

def plot_PCA(matrix, metadata, title):
    '''
    PCA plot
    :param matrix: processed data matrix
    :return: PCA plot
    '''

    pca = PCA(n_components=100)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    group = np.array([metadata[i][0] for i in samples])
    cdict = {'Nonuser': 'red', 'E+P': 'blue', 'E-only': 'green'}
    fig, ax = plt.subplots()
    for g in np.unique(group):
        ix = np.where(group == g)
        ax.scatter(scatter_x[ix], scatter_y[ix], c=cdict[g], label=g, s=10)
    ax.legend()
    plt.style.use("ggplot")
    plt.title("PCA for " + title)
    plt.xlabel("Component 1:  " + str(round(pca.explained_variance_ratio_[0], 2)) + "%")
    plt.ylabel("Component 2: " + str(round(pca.explained_variance_ratio_[1], 2)) + "%")
    plt.savefig(title + ".png")
    plt.show()