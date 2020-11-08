import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import utils

mat = pd.read_csv("miyamoto_data.csv", index_col=0).T
print(mat.head)

md_raw = pd.read_csv("Miyamoto_metadata.csv")

metadata_list = list(zip(md_raw['Organ'], md_raw['Cancer status'], md_raw['Smoker'], md_raw['Gender']))
metadata_dict = dict(zip(md_raw['Sample name'].values, metadata_list))

matrix_proc = utils.data_processing(mat, 0)



utils.plot_PCA(matrix_proc, metadata_dict, "Miyamoto (organ)", 0)
