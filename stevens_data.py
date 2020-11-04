import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Get metadata
md_raw = pd.read_csv("../Stevens/MTBLS136_compressed_files/s_MTBLS136.txt", sep="\t")
print(md_raw.columns)
metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'], md_raw['Factor Value[AgeAtBloodDraw]']))
metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))
print(metadata_dict)

# Get abundance matrix
# mat = pd.read_excel("../Stevens/stevens_peak_table.xlsx", header=3)
# mat.to_csv("Stevens_matrix.csv", index=False)
# print(mat.head)
mat = pd.read_csv("Stevens_matrix.csv")

mat_values = mat.iloc[:,9:]
print(mat_values.head)
