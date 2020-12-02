import pandas as pd
import numpy as np
import utils

KEGG_pathways = pd.read_csv("KEGG_mouse_pathways_compounds.csv", dtype=str, index_col=0)

def brown_data():
    mat = pd.read_excel("../Brown_mouse_diet/abundance.xlsx", index_col=0, header=1).T

    metadata = pd.read_csv("../Brown_mouse_diet/s_metabolon.txt", sep="\t")
    sample_name = [i[0:10] for i in metadata["Sample Name"]]
    diet = metadata["Factor Value[Genotype]"].tolist()
    metadata_dict = dict(zip(sample_name, diet))
    mat_proc = utils.data_processing(mat, firstrow=6, firstcol=1)


    mat_proc["Group"] = mat_proc.index.map(metadata_dict)
    utils.plot_PCA(mat_proc.iloc[:, :-1], mat_proc["Group"], "Mouse diet", 10)

    ttest_res = utils.t_tests(mat_proc.iloc[:, :-1], mat_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    mat = mat.T

    DEM_KEGG_id = mat[mat.index.isin(DEM)]['KEGG'].tolist()
    background_KEGG = mat['KEGG'].tolist()
    return DEM_KEGG_id, background_KEGG
DEM_KEGG_id, background_KEGG = brown_data()
print(len(DEM_KEGG_id))
print(len(background_KEGG))

ORA_res = utils.over_representation_analysis(DEM_KEGG_id, background_KEGG, KEGG_pathways)
# ORA_res.to_csv("../Brown_mouse_diet/ORA_KEGG.csv")

print(len(ORA_res[ORA_res["P-value"] < 0.1]["P-value"].tolist()))
print(len(ORA_res[ORA_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
# ora_res.to_csv("../Yamada/Yamada_ORA_human.csv", index=False)

all_KEGG_compounds = list(set([x for x in KEGG_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
print(len(all_KEGG_compounds))

ora_whole_KEGG_bg = utils.over_representation_analysis(DEM_KEGG_id, all_KEGG_compounds, KEGG_pathways)
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-value"] < 0.1]["P-value"].tolist()))
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-adjust"] < 0.1]["P-adjust"].tolist()))