import pandas as pd
import utils

data = pd.read_excel("../Yamada/Yamada.xlsx", index_col=0, header=0).T
data.columns = [col[0:6] for col in data.columns]
CRC_or_healthy = ["Healthy" if i in ["Healthy"] else "CRC" if i in ["Stage_I_II", "Stage_III_IV"] else "Null" for i in data["Group"]]
BMI = ["Healthy" if i < 24 else "Overweight" for i in data["BMI"]]
age = ["Over60" if i > 60 else "Under60" for i in data["Age"]]
data.insert(1, "binary", CRC_or_healthy)

data = data.iloc[:, ~data.columns.duplicated()]
df = data[data.binary != "Null"]
data_proc = utils.data_processing(df, 0, 5)
ttest_res = utils.t_tests(data_proc, data[data.binary != "Null"]["binary"], "fdr_bh")
DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
print(DEM)
KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)

print(len(DEM), "Differential metabolites")
background = data_proc.columns.tolist()
print(len(background), "in background list")
ora_res = utils.over_representation_analysis(DEM, background, KEGG_pathways)
print(ora_res.head)
print(len(ora_res[ora_res["P-adjust"] < 0.2]["P-adjust"].tolist()))
# ora_res.to_csv("../Yamada/Yamada_ORA_human.csv", index=False)
