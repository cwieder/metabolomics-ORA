import pandas as pd
import numpy as np
import utils

pd.options.mode.chained_assignment = None  # default='warn'

n_zscore = pd.read_csv("../Zamboni/mod_zscore_neg_CW.csv", index_col=0)
p_zscore = pd.read_csv("../Zamboni/mod_zscore_pos_CW.csv", index_col=0)

# remove unannotated
n_zscore = n_zscore[n_zscore.index.notnull()]
p_zscore = p_zscore[p_zscore.index.notnull()]

# get all possible annotations for each ion
putative_annotations = pd.read_excel("../Zamboni/annotations_EV1.xlsx", sheet_name="Table EV1B")

neg_putative_annotations = putative_annotations[putative_annotations['Ionization Mode'] == "neg"]
pos_putative_annotations = putative_annotations[putative_annotations['Ionization Mode'] == "pos"]
annotations_neg = dict.fromkeys(neg_putative_annotations["Ion Index"])
annotations_pos = dict.fromkeys(pos_putative_annotations["Ion Index"])
for row in neg_putative_annotations.itertuples():
    annotation_cols = row[4:]
    ion_index = row[2]
    annotations = []
    for i in annotation_cols:
        if not pd.isna(i):
            annos = i.split()
            annotations.append(annos[-2])
    annotations_neg[ion_index] = annotations

for row in pos_putative_annotations.itertuples():
    annotation_cols = row[4:]
    ion_index = row[2]
    annotations = []
    for i in annotation_cols:
        if not pd.isna(i):
            annos = i.split()
            annotations.append(annos[-2])
    annotations_pos[ion_index] = annotations

# remove empty annotations
# empty_neg =
# empty_pos =

# filter z score
strain_DA_compounds = dict.fromkeys(n_zscore.columns)

for strain in strain_DA_compounds.keys():
    cur_col = n_zscore.loc[:, strain]
    DA_metabolites = []
    for items in cur_col.iteritems():
        if items[1] > 6 or items[1] < -6:
            DA_metabolites.append(annotations_neg[items[0]])
    DA_KEGG = [j for i in DA_metabolites for j in i]
    strain_DA_compounds[strain] = DA_KEGG

for strain in strain_DA_compounds.keys():
    cur_col = p_zscore.loc[:, strain]
    DA_metabolites = []
    for items in cur_col.iteritems():
        if items[1] > 6 or items[1] < -6:
            DA_metabolites.append(annotations_pos[items[0]])
    DA_KEGG = [j for i in DA_metabolites for j in i]
    strain_DA_compounds[strain] = DA_KEGG

# for k, v in strain_DA_compounds.items():
#     print(len(v))

# ORA

background_list_all_annotations = list(set(sum(annotations_neg.values(), []) + sum(annotations_pos.values(), [])))
print(len(background_list_all_annotations))
KEGG_pathways = pd.read_csv("KEGG_ecoMG1655_pathways_compounds.csv", dtype=str, index_col=0)
#
# n_significant = []
# for strain in strain_DA_compounds.keys():
#     DA_met_unique = list(set(strain_DA_compounds[strain]))
#     if DA_met_unique:
#         ora_res = utils.over_representation_analysis(DA_met_unique, background_list_all_annotations, KEGG_pathways)
#         if len(ora_res[ora_res["P-adjust"] < 0.05]["P-adjust"].tolist()) > 5:
#             print(strain, len(DA_met_unique), len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
#         n_significant.append(len(ora_res[ora_res["P-adjust"] < 0.05]["P-adjust"].tolist()))
#     else:
#         continue
# print(n_significant)


print(len(list(set(strain_DA_compounds["arcB"]))))
ora_res = utils.over_representation_analysis(list(set(strain_DA_compounds["arcB"])), background_list_all_annotations, KEGG_pathways)

# ora_res.to_csv("../Zamboni/creB_ORA.csv", index=False)
print(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
print(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
# ora_res.to_csv("../Yamada/Yamada_ORA_human.csv", index=False)

all_KEGG_compounds = list(set([x for x in KEGG_pathways.iloc[:, 1:].values.flatten() if x is not np.nan]))
print(len(all_KEGG_compounds))

ora_whole_KEGG_bg = utils.over_representation_analysis(list(set(strain_DA_compounds["arcB"])), all_KEGG_compounds, KEGG_pathways)
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-value"] < 0.1]["P-value"].tolist()))
print(len(ora_whole_KEGG_bg[ora_whole_KEGG_bg["P-adjust"] < 0.1]["P-adjust"].tolist()))
