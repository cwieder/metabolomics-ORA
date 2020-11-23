# Gut microbiome structure and metabolic activity in inflammatory bowel disease
# Xavier et al

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats
import utils


xavier_data = pd.read_excel("../Xavier/Xavier_data.xlsx", index_col=0, header=1).T
xavier_data_proc = utils.data_processing(xavier_data, 0, 7)
utils.plot_PCA(xavier_data_proc, xavier_data["Diagnosis"], "Xavier")

