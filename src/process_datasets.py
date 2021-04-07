import pandas as pd
import numpy as np
import utils
import re
from bioservices import *
import pickle

kegg_db = KEGG(verbose=False)

# with open('../data/MetaCyc_compound_mapping.pickle', 'rb') as handle:
#     metacyc_mapping = pickle.load(handle)
# IMPORT DATASETS, PRE-PROCESS THEM AND RUN T-TESTS TO OBTAIN LIST OF DIFFERENTIALLY ABUNDANT METABOLITES

def yamada_data(db="KEGG"):
    data = pd.read_csv("../example_data/yachida_abundance.csv", index_col=0, header=0).T
    data = data.rename(columns={'Group': 'disease'})
    sample_disease_dict = dict(zip(data.index, data['disease']))
    data.columns = data.columns[0:4].tolist() + [col[0:6] for col in data.columns[4:]]

    removecols = []
    for i in data.columns.tolist():
        matchObj = re.search("^[C]\d{5}$", i)
        if not matchObj:
            removecols.append(i)

    data = data.drop(removecols[4:], axis=1)
    # TODO: Try and map these compounds

    CRC_or_healthy_dict = dict.fromkeys(data.index.tolist())
    for k, v in sample_disease_dict.items():
        if v in ['Healthy']:
            CRC_or_healthy_dict[k] = "Healthy"
        elif v in ["Stage_I_II", "Stage_III_IV"]:
            CRC_or_healthy_dict[k] = "CRC"
        else:
            CRC_or_healthy_dict[k] = "Null"
    CRC_or_healthy = ["Healthy" if i in ["Healthy"] else "CRC" if i in ["Stage_I_II", "Stage_III_IV"] else "Null" for i in data["disease"]]

    data.insert(1, "Group", CRC_or_healthy)

    data = data.iloc[:, ~data.columns.duplicated()]
    df = data[data.disease != "Null"]
    data_proc = utils.data_processing(df, 0, 5)
    data_proc["Group"] = data_proc.index.map(CRC_or_healthy_dict)
    data_proc = data_proc[data_proc.Group != "Null"]

    if db == "Reactome":
        map_kegg_chebi = kegg_db.conv("chebi", "compound")
        KEGG2Reactome = dict.fromkeys(data_proc.columns.tolist()[:-1])
        for cpd in KEGG2Reactome.keys():
            chebiID = map_kegg_chebi['cpd:' + cpd]
            KEGG2Reactome[cpd] = chebiID[6:]
        data_proc = data_proc.rename(columns=KEGG2Reactome)

    if db == "Cyc":
        mapping = pd.read_csv("../yamada2metacyc.txt", sep="\t")
        KEGG2BioCyc = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
        data_proc = data_proc.rename(columns=KEGG2BioCyc)
        data_proc = data_proc[data_proc.columns.dropna()]

    ttest_res = utils.t_tests(data_proc.iloc[:, :-1], data_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    background = data_proc.iloc[:, :-1].columns.tolist()
    return DEM, background, data_proc

def brown_data(db="KEGG"):
    # mat = pd.read_excel("../example_data/Labbe_abundance.xlsx", index_col=0, header=1, sheet_name="OrigScale(tissue)", engine='openpyxl').T
    mat = pd.read_csv("../example_data/labbe_abundance.csv", index_col=0, header=1).T
    mapping = dict(zip(mat.columns.tolist(), mat.loc["KEGG", :].tolist()))
    mat = mat.rename(columns=mapping)
    mat = mat.loc[:, mat.columns.notnull()]
    mat = mat.loc[:, ~mat.columns.duplicated()]
    metadata = pd.read_csv("../example_data/Labbe_metadata.txt", sep="\t")
    sample_name = [i[0:10] for i in metadata["Sample Name"]]
    diet = metadata["Factor Value[Genotype]"].tolist()
    metadata_dict = dict(zip(sample_name, diet))
    mat_proc = utils.data_processing(mat, firstrow=6, firstcol=1)
    mat_proc["Group"] = mat_proc.index.map(metadata_dict)
    mat_proc = mat_proc.iloc[:, ~mat_proc.columns.duplicated()]

    if db == "Reactome":
        map_kegg_chebi = kegg_db.conv("chebi", "compound")
        KEGG2Reactome = dict.fromkeys(mat_proc.columns.tolist()[:-1])
        for cpd in KEGG2Reactome.keys():
            try:
                chebiID = map_kegg_chebi['cpd:' + cpd]
                KEGG2Reactome[cpd] = chebiID[6:]
            except KeyError:
                KEGG2Reactome[cpd] = np.nan
        mat_proc = mat_proc.rename(columns=KEGG2Reactome)
        mat_proc = mat_proc.loc[:, mat_proc.columns.notnull()]

    if db == "Cyc":
        mapping = pd.read_csv("../brown2metacyc.txt", sep="\t")
        KEGG2BioCyc = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
        mat_proc = mat_proc.rename(columns=KEGG2BioCyc)
        mat_proc = mat_proc[mat_proc.columns.dropna()]

    ttest_res = utils.t_tests(mat_proc.iloc[:, :-1], mat_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].dropna().tolist()
    background = mat_proc.columns.tolist()[:-1]
    return DEM, background, mat_proc

def stevens_data(db="KEGG"):
    md_raw = pd.read_csv("../example_data/Stevens_metadata.txt", sep="\t")
    metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'],
                             md_raw['Factor Value[AgeAtBloodDraw]'],
                             ['Over 75' if val not in ['<=55', '56-60', '61-65', '66-70', '71-75'] else 'Under 75' for val in md_raw['Factor Value[AgeAtBloodDraw]']]))
    metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))
    sample_status_dict = dict(zip(md_raw['Sample Name'].values, md_raw['Factor Value[CurrentPMH]']))

    replicate_samples = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', 'E+P']]
    nonusers = [k for k, v in metadata_dict.items() if v[0] not in [np.nan, 'E-only', 'E+P']]
    estrogen_only = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', np.nan, 'E+P']]
    estrogen_progesterone = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', np.nan]]

    # Get abundance matrix, transpose to n-samples by m-metabolites
    mat = pd.read_csv("../example_data/Stevens_matrix_named_compounds_only.csv", index_col=0)
    mat_nonusers_estrogen = mat.drop((replicate_samples + estrogen_progesterone), axis=1)
    stevens_matrix_proc = utils.data_processing(mat_nonusers_estrogen.T, 8, 0)
    stevens_matrix_proc["Group"] = stevens_matrix_proc.index.map(sample_status_dict)
    mapping_dict = dict(zip(mat.index, mat['KEGG']))
    stevens_matrix_proc = stevens_matrix_proc.rename(columns=mapping_dict)

    if db == "Reactome":
        map_kegg_chebi = kegg_db.conv("chebi", "compound")
        KEGG2Reactome = dict.fromkeys(stevens_matrix_proc.columns.tolist()[:-1])
        for cpd in KEGG2Reactome.keys():
            try:
                chebiID = map_kegg_chebi['cpd:' + str(cpd)]
                KEGG2Reactome[cpd] = chebiID[6:]
            except KeyError:
                KEGG2Reactome[cpd] = np.nan
        stevens_matrix_proc = stevens_matrix_proc.rename(columns=KEGG2Reactome)
        stevens_matrix_proc = stevens_matrix_proc.loc[:, stevens_matrix_proc.columns.notnull()]

    stevens_matrix_proc = stevens_matrix_proc.loc[:, stevens_matrix_proc.columns.notnull()]
    stevens_matrix_proc = stevens_matrix_proc.loc[:, ~stevens_matrix_proc.columns.duplicated()]

    if db == "Cyc":
        mapping = pd.read_csv("../stevens2metacyc.txt", sep="\t")
        KEGG2BioCyc = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
        stevens_matrix_proc = stevens_matrix_proc.rename(columns=KEGG2BioCyc)
        stevens_matrix_proc = stevens_matrix_proc[stevens_matrix_proc.columns.dropna()]

    ttest_res = utils.t_tests(stevens_matrix_proc.iloc[:, :-1], stevens_matrix_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    background = stevens_matrix_proc.iloc[:, :-1].columns.tolist()
    return DEM, background, stevens_matrix_proc


def zamboni_data(knockout, db="KEGG"):
    # import modified z-scores
    n_zscore = pd.read_csv("../example_data/Fuhrer_mod_zscore_neg_CW.csv.zip", index_col=0)
    p_zscore = pd.read_csv("../example_data/Fuhrer_mod_zscore_pos_CW.csv.zip", index_col=0)

    # remove unannotated
    n_zscore = n_zscore[n_zscore.index.notnull()]
    p_zscore = p_zscore[p_zscore.index.notnull()]

    # get all possible annotations for each ion
    # putative_annotations = pd.read_excel("../Zamboni/annotations_EV1.xlsx", sheet_name="Table EV1B")
    #
    # neg_putative_annotations = putative_annotations[putative_annotations['Ionization Mode'] == "neg"]
    # pos_putative_annotations = putative_annotations[putative_annotations['Ionization Mode'] == "pos"]
    # annotations_neg = dict.fromkeys(neg_putative_annotations["Ion Index"])
    # annotations_pos = dict.fromkeys(pos_putative_annotations["Ion Index"])
    # for row in neg_putative_annotations.itertuples():
    #     annotation_cols = row[4:]
    #     ion_index = row[2]
    #     annotations = []
    #     for i in annotation_cols:
    #         if not pd.isna(i):
    #             annos = i.split()
    #             annotations.append(annos[-2])
    #     annotations_neg[ion_index] = annotations
    #
    # for row in pos_putative_annotations.itertuples():
    #     annotation_cols = row[4:]
    #     ion_index = row[2]
    #     annotations = []
    #     for i in annotation_cols:
    #         if not pd.isna(i):
    #             annos = i.split()
    #             annotations.append(annos[-2])
    #     annotations_pos[ion_index] = annotations

    with open('../data/zamboni_pos_annotation_dict2.pickle', 'rb') as handle:
        annotations_pos = pickle.load(handle)

    with open('../data/zamboni_neg_annotation_dict2.pickle', 'rb') as handle:
        annotations_neg = pickle.load(handle)

    # convert to CHEBI for Reactome
    if db == "Reactome":
        # kegg_id_neg = sorted({x for v in annotations_neg.values() for x in v})
        # kegg_id_pos = sorted({x for v in annotations_pos.values() for x in v})
        # kegg_id = list(set(kegg_id_pos + kegg_id_neg))
        # mapping_dict = dict.fromkeys(kegg_id)
        #
        # for num, cpd in enumerate(mapping_dict.keys()):
        #     map_kegg_chebi = kegg_db.conv("chebi", "compound")
        #     try:
        #         chebiID = map_kegg_chebi['cpd:' + cpd]
        #         mapping_dict[cpd] = chebiID
        #     except KeyError:
        #         mapping_dict[cpd] = np.nan
        # with open('zamboni_CHEBI_mapping.pickle', 'wb') as handle:
        #     pickle.dump(mapping_dict, handle, protocol=4)

        with open('../data/zamboni_CHEBI_mapping.pickle', 'rb') as handle:
            mapping_dict = pickle.load(handle)
            annotations_neg = {k: [str(mapping_dict[i])[6:] for i in v] for k, v in annotations_neg.items()}
            annotations_pos = {k: [str(mapping_dict[i])[6:] for i in v] for k, v in annotations_pos.items()}

    if db == "Cyc":
        mapping = pd.read_csv("../zamboni2metacyc.txt", sep="\t")
        mapping_dict = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
        annotations_neg = {k: [mapping_dict[i] if i in mapping_dict.keys() else "" for i in v] for k, v in annotations_neg.items()}
        annotations_pos = {k: [mapping_dict[i] if i in mapping_dict.keys() else "" for i in v] for k, v in
                           annotations_pos.items()}

        # annotations_neg = {k: ([str(mapping_dict[i]) for i in v] if mapping_dict[i] else "Null") for k, v in annotations_neg.items() else "Null"}
        # annotations_pos = {k: [str(mapping_dict[i]) for i in v] for k, v in annotations_pos.items()}
    strain_DA_compounds = dict.fromkeys(n_zscore.columns)

    for strain in strain_DA_compounds.keys():
        cur_col = n_zscore.loc[:, strain]
        DA_metabolites = []
        for items in cur_col.iteritems():
            if items[1] > 6 or items[1] < -6:
                DA_metabolites.append(annotations_neg[items[0]])
        DA_KEGG = [j for i in DA_metabolites for j in i] # filter out only annotated compounds
        strain_DA_compounds[strain] = DA_KEGG

    for strain in strain_DA_compounds.keys():
        cur_col = p_zscore.loc[:, strain]
        DA_metabolites = []
        for items in cur_col.iteritems():
            if items[1] > 6 or items[1] < -6:
                DA_metabolites.append(annotations_pos[items[0]])
        DA_KEGG = [j for i in DA_metabolites for j in i] # filter out only annotated compounds
        strain_DA_compounds[strain] += DA_KEGG
    background_list_all_annotations = list(filter(None, list(set(sum(annotations_neg.values(), []) + sum(annotations_pos.values(), [])))))
    DEM = list(filter(None, list(set(strain_DA_compounds[knockout]))))

    # return strain matrix with all annotations
    empty_anno_pos = [k for k, v in annotations_pos.items() if not v]

    strain_pos_all = p_zscore.drop(empty_anno_pos).loc[:, knockout]
    pos_rows = []
    for ion, val in strain_pos_all.items():
        ion_annos = annotations_pos[ion]
        for i in ion_annos:
            row = ("pos"+str(ion), i, val)
            pos_rows.append(row)
    pos_annos_df = pd.DataFrame(pos_rows, columns=['Ion', 'Annotation', 'Z-score'])

    empty_anno_neg = [k for k, v in annotations_neg.items() if not v]
    strain_neg_all = n_zscore.drop(empty_anno_neg).loc[:, knockout]
    neg_rows = []
    for ion, val in strain_neg_all.items():
        ion_annos = annotations_neg[ion]
        for i in ion_annos:
            row = ("neg"+str(ion), i, val)
            neg_rows.append(row)
    neg_annos_df = pd.DataFrame(neg_rows, columns=['Ion', 'Annotation', 'Z-score'])
    all_annotations_df = pd.concat([pos_annos_df, neg_annos_df])
    all_annotations_df = all_annotations_df.set_index(all_annotations_df["Annotation"])
    mat = all_annotations_df.iloc[:, 2:].T
    # dem_count = []
    # for x in all_annotations_df.itertuples():
    #     if x[3] > 6 or x[3] < -6:
    #         dem_count.append(x[2])
    # mat_pval = mat
    # mat_pval.iloc[0, :] = [stats.norm.cdf(i) for i in mat_pval.iloc[0, :]]
    # mat_pval.T.to_csv("Zamboni_dcuS_for_IPA.csv")
    return DEM, background_list_all_annotations, mat


def auwerx_data(db="KEGG"):
    mat = pd.read_csv("../example_data/quiros_abundance.csv", index_col=6).T
    mat = mat.iloc[6:, :]
    mat = mat.loc[:, ~mat.columns.duplicated(keep='first')]
    mat = mat.loc[:, mat.columns.notnull()]
    groups = [i.split(".", 1)[0] for i in mat.index.tolist()]
    mat['Group'] = groups
    mat_selected_groups = mat.loc[mat['Group'].isin(['Acti', 'Dox'])] # Acti vs Dox
    matrix_proc = utils.data_processing(mat_selected_groups.iloc[:, :-1], 0, 0)
    matrix_proc_copy = matrix_proc.copy()
    matrix_proc_copy['Group'] = mat_selected_groups['Group']

    if db == "Reactome":
        map_kegg_chebi = kegg_db.conv("chebi", "compound")
        KEGG2Reactome = dict.fromkeys(matrix_proc_copy.columns.tolist()[:-1])
        for cpd in KEGG2Reactome.keys():
            try:
                chebiID = map_kegg_chebi['cpd:' + str(cpd)]
                KEGG2Reactome[cpd] = chebiID[6:]
            except KeyError:
                KEGG2Reactome[cpd] = np.nan
        matrix_proc_copy = matrix_proc_copy.rename(columns=KEGG2Reactome)
        matrix_proc_copy = matrix_proc_copy.loc[:, matrix_proc_copy.columns.notnull()]

    if db == "Cyc":
        mapping = pd.read_csv("../auwerx2metacyc.txt", sep="\t")
        KEGG2BioCyc = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
        matrix_proc_copy = matrix_proc_copy.rename(columns=KEGG2BioCyc)
        matrix_proc_copy = matrix_proc_copy[matrix_proc_copy.columns.dropna()]

    ttest_res = utils.t_tests(matrix_proc_copy.iloc[:, :-1], matrix_proc_copy["Group"], "fdr_bh")
    DA_metabolites = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    background = matrix_proc_copy.iloc[:, :-1].columns.tolist()
    return DA_metabolites, background, matrix_proc_copy
