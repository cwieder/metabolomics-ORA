import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats

# Get metadata
md_raw = pd.read_csv("../Stevens/MTBLS136_compressed_files/s_MTBLS136.txt", sep="\t")

metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'], md_raw['Factor Value[AgeAtBloodDraw]']))
metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))

replicate_samples = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', 'E+P']]
nonusers = [k for k, v in metadata_dict.items() if v[0] not in [np.nan, 'E-only', 'E+P']]
estrogen_only = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', np.nan, 'E+P']]
estrogen_progesterone = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', np.nan]]

# Get abundance matrix, transpose to n-samples by m-metabolites
mat = pd.read_csv("Stevens_matrix_named_compounds_only.csv", index_col=0)

def data_processing(raw_matrix, firstrow):
    '''
    Filtering low abundance metabolites, data cleaning and imputation using minimum value.
    :param raw_matrix: raw abundance matrix
    :param firstrow: First row containing values (integer)
    :param transpose: transpose input matrix
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
    uniq_sample = np.unique(group)
    cmap = ['tab:green', 'tab:orange', 'tab:blue']
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
    #plt.savefig(title + ".png")
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



stevens_matrix_proc = data_processing(mat.T, 8)
mat_nonusers_estrogen = stevens_matrix_proc.drop((replicate_samples + estrogen_progesterone), axis=0)
#plot_PCA(mat_nonusers_estrogen, metadata_dict, "non-users vs. estrogen")

differential_metabolites = linear_regression(mat_nonusers_estrogen)
DEM = differential_metabolites[differential_metabolites["P-adjust"] < 0.05]["Metabolite"].tolist()
mat_nonusers_estrogen_DEM = mat_nonusers_estrogen[DEM]
# plot_PCA(mat_nonusers_estrogen_DEM, metadata_dict, "Non-user vs estrogen (differentially expressed metabolites)")

KEGG_pathways = pd.read_csv("KEGG_reference_pathways_compounds.csv", dtype=str, index_col=0)
DEM_KEGG_id = mat[mat.index.isin(DEM)]['KEGG'].tolist()
DEM_KEGG_id = [i for i in DEM_KEGG_id if str(i) != 'nan']
print(DEM_KEGG_id)

with open("DEM_Stevens.txt", "a") as outfile:
    for i in DEM_KEGG_id:
        outfile.write(i + "\n")

stevens_background_list = mat['KEGG'].dropna().tolist()
print(len(stevens_background_list))

with open("background_list_Stevens.txt", "a") as outfile_bglist:
    for i in stevens_background_list:
        outfile_bglist.write(i + "\n")

print(len(DEM_KEGG_id))
ORA_res = over_representation_analysis(DEM_KEGG_id, stevens_background_list, KEGG_pathways)
ORA_res.to_csv("ORA_Stevens.csv")
# differential_metabolites.to_csv("differential_metabolites_nonuser_v_estrogen.csv")


# mat_nonusers_e_p = stevens_matrix_proc.drop((replicate_samples + estrogen_only), axis=0)
# mat_estrogen_vs_e_p = stevens_matrix_proc.drop((replicate_samples + nonusers), axis=1)
# mat_all = stevens_matrix_proc.drop(replicate_samples, axis=0)

# plot_PCA(mat_all, metadata_dict, "all groups")
# plot_PCA(mat_nonusers_e_p, metadata_dict, "nonusers vs. estrogen+progestin")