# Repeat for positive

import pandas as pd
raw_data_neg = pd.read_csv("../Zamboni/rawdata_pos_all.tsv", sep="\t", index_col=0)

strains = pd.read_excel("../Zamboni/sample_id.xls", header=None)
strains_list = strains.iloc[:,0].tolist()
col_list = []
for strain in strains_list:
    col_list.append([i for i in raw_data_neg.columns if i.startswith(strain)])

average_series = []
for strain in col_list:
    avg_strain = raw_data_neg.loc[:, strain].mean(axis=1)
    avg_strain.name = strain[0]
    average_series.append(avg_strain)

avg_df = pd.concat(average_series, axis=1)
# check for outliers per row
# Q1 = avg_df.quantile(0.25)
# Q3 = avg_df.quantile(0.75)
# IQR = Q3 - Q1
# print(avg_df.shape)
# avg_df_out = avg_df[~((avg_df < (Q1 - 1.5 * IQR)) |(avg_df > (Q3 + 1.5 * IQR))).any(axis=1)]
# print(avg_df_out.shape)

# calculate modified z-score
# Z = (x – median(x) ) / sd(x)

mod_zscore_df = avg_df.apply(lambda x: round((x - x.median(axis=0))/x.mad(axis=0),5), axis=1)
print(mod_zscore_df.head)
mod_zscore_df.to_csv("../Zamboni/mod_zscore_pos_CW.csv")

mod_zscore_df = mod_zscore_df[mod_zscore_df.index.notnull()]
print(mod_zscore_df.shape)