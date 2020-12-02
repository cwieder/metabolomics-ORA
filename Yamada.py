import pandas as pd
import numpy as np
import utils

KEGG_pathways = pd.read_csv("KEGG_human_pathways_compounds.csv", dtype=str, index_col=0)

def yamada_data():
    data = pd.read_excel("../Yamada/Yamada.xlsx", index_col=0, header=0).T
    data.columns = [col[0:6] for col in data.columns]
    CRC_or_healthy = ["Healthy" if i in ["Healthy"] else "CRC" if i in ["Stage_I_II", "Stage_III_IV"] else "Null" for i in data["Group"]]
    # BMI = ["Healthy" if i < 24 else "Overweight" for i in data["BMI"]]
    # age = ["Over60" if i > 60 else "Under60" for i in data["Age"]]
    data.insert(1, "binary", CRC_or_healthy)

    data = data.iloc[:, ~data.columns.duplicated()]
    df = data[data.binary != "Null"]
    data_proc = utils.data_processing(df, 0, 5)
    ttest_res = utils.t_tests(data_proc, data[data.binary != "Null"]["binary"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()

    background = data_proc.columns.tolist()
    return DEM, background

DEM, background = yamada_data()
# ora_res = utils.over_representation_analysis(DEM, background, KEGG_pathways)
# print(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
# print(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
# # ora_res.to_csv("../Yamada/Yamada_ORA_human.csv", index=False)
#
# all_KEGG_compounds = list(set([x for x in KEGG_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
# print(len(all_KEGG_compounds))
#
# ora_whole_KEGG_bg = utils.over_representation_analysis(DEM, all_KEGG_compounds, KEGG_pathways)
# print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-value"] < 0.1]["P-value"].tolist()))
# print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-adjust"] < 0.1]["P-adjust"].tolist()))

ora1 = utils.reduce_background_list_ora(background, 70, DEM, KEGG_pathways)
print(ora1)