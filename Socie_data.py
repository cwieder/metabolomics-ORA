import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils

def plot_PCA(matrix, metadata, title, n_comp=100):
    '''
    PCA plot
    :param matrix: processed data matrix
    :param labelpos: position in metadata list of sample labels of interest
    :return: PCA plot
    '''

    pca = PCA(n_components=n_comp)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    print(samples)
    group = np.array([metadata[i] for i in samples])
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
    plt.savefig(title + ".png")
    plt.show()


# Get metadata
md_raw = pd.read_csv("../Socie/s_MTBLS204.txt", sep="\t")

metadata = dict(zip(md_raw["Source Name"], md_raw["Factor Value[Class]"]))

socie_mat = pd.read_excel("../Socie/abundance_mat.xlsx", index_col=0)
# remove_cols = [i if metadata[i] not in ["acute GVHD", "no GVHD"] else None for i in socie_mat.columns[6:]]
# remove_cols = [i for i in remove_cols if i != None]
# socie_mat = socie_mat.drop(remove_cols, axis=1)

socie_mat_proc = utils.data_processing(socie_mat.T, 6, 1)

# socie_mat["Group"] = socie_mat.iloc[:,1].map(metadata)
# print(socie_mat["Group"])

# plot_PCA(socie_mat_proc, metadata, "Donor vs Recipient", n_comp=10)
socie_mat_proc["Group"] = socie_mat_proc.index.map(metadata)

ttest_res = utils.t_tests(socie_mat_proc.iloc[:,:-1], socie_mat_proc["Group"], "fdr_bh")
# ttest_res.to_csv("DEM_Stevens_ttest.csv")

# DEM = differential_metabolites[differential_metabolites["P-adjust"] < 0.05]["Metabolite"].tolist()
DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
print(DEM)
# with open("DEM_Stevens_ttest.txt", "a") as infile:
#     for i in DEM:
#         infile.write(i+"\n")
print(len(DEM), "Differential metabolites")

KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)
socie_mat.columns = socie_mat.columns.str.replace(' ', '')
metadata_cols = ['KEGG', 'GroupHMDB_ID']
socie_mat_proc_annotated = socie_mat_proc.T.join(socie_mat[metadata_cols])

DEM_KEGG_id = socie_mat[socie_mat.index.isin(DEM)]['KEGG'].tolist()
DEM_KEGG_id = [i for i in DEM_KEGG_id if str(i) != 'nan']
print(DEM_KEGG_id)
socie_mat[socie_mat.index.isin(DEM)]['KEGG'].dropna().to_csv("Socie_DEM_KEGG.csv", index=False)
print(len(DEM_KEGG_id), "Differential metabolites mapping to KEGG pathways")
background_list = socie_mat_proc_annotated['KEGG'].dropna().tolist()
print(len(background_list), "background in KEGG")
print(socie_mat_proc.shape)

ORA_res = utils.over_representation_analysis(DEM_KEGG_id, background_list, KEGG_pathways)
print(ORA_res[ORA_res['P-adjust'] < 0.2].count())
print(ORA_res)
ORA_res.to_csv("../Socie/ORA_Socie_ttest_FDR.csv")