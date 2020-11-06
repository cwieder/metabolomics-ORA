import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import seaborn as sns

# Get metadata
md_raw = pd.read_csv("../Stevens/MTBLS136_compressed_files/s_MTBLS136.txt", sep="\t")

metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'], md_raw['Factor Value[AgeAtBloodDraw]']))
metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))

replicate_samples = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', 'E+P']]
nonusers = [k for k, v in metadata_dict.items() if v[0] not in [np.nan, 'E-only', 'E+P']]
estrogen_only = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', np.nan, 'E+P']]
estrogen_progesterone = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', np.nan]]

# Get abundance matrix, transpose to n-samples by m-metabolites
mat = pd.read_csv("Stevens_matrix_named_compounds_only.csv", index_col=0).T

def data_processing(raw_matrix, firstrow):
    '''
    Filtering low abundance metabolites, data cleaning and imputation using minimum value.
    :param raw_matrix: raw abundance matrix
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
    :param matrix:
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

stevens_matrix_proc = data_processing(mat, 8)

mat_nonusers_estrogen = stevens_matrix_proc.drop((replicate_samples + estrogen_progesterone), axis=0)
mat_nonusers_e_p = stevens_matrix_proc.drop((replicate_samples + estrogen_only), axis=0)
# mat_estrogen_vs_e_p = stevens_matrix_proc.drop((replicate_samples + nonusers), axis=1)
mat_all = stevens_matrix_proc.drop(replicate_samples, axis=0)

plot_PCA(mat_all, metadata_dict, "all groups")
plot_PCA(mat_nonusers_estrogen, metadata_dict, "non-users vs. estrogen")
# plot_PCA(mat_estrogen_vs_e_p, metadata_dict)
plot_PCA(mat_nonusers_e_p, metadata_dict, "nonusers vs. estrogen+progestin")