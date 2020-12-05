import pandas as pd
import numpy as np
import utils

# IMPORT DATASETS, PRE-PROCESS THEM AND RUN T-TESTS TO OBTAIN LIST OF DIFFERENTIALLY ABUNDANT METABOLITES

def yamada_data():
    data = pd.read_excel("../Yamada/Yamada.xlsx", index_col=0, header=0).T
    data = data.rename(columns={'Group': 'disease'})
    sample_disease_dict = dict(zip(data.index, data['disease']))
    data.columns = data.columns[0:4].tolist() + [col[0:6] for col in data.columns[4:]]
    remove_not_in_KEGG = [i for i in data.columns[4:] if i[0] != "C"]
    data = data.drop(remove_not_in_KEGG, axis=1)
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
    ttest_res = utils.t_tests(data_proc.iloc[:, :-1], data_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    background = data_proc.iloc[:, :-1].columns.tolist()
    return DEM, background, data_proc

def brown_data():
    mat = pd.read_excel("../Brown_mouse_diet/abundance.xlsx", index_col=0, header=1).T

    metadata = pd.read_csv("../Brown_mouse_diet/s_metabolon.txt", sep="\t")
    sample_name = [i[0:10] for i in metadata["Sample Name"]]
    diet = metadata["Factor Value[Genotype]"].tolist()
    metadata_dict = dict(zip(sample_name, diet))
    mat_proc = utils.data_processing(mat, firstrow=6, firstcol=1)
    mat_proc["Group"] = mat_proc.index.map(metadata_dict)

    ttest_res = utils.t_tests(mat_proc.iloc[:, :-1], mat_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
    mat = mat.T
    DEM_KEGG_id = mat[mat.index.isin(DEM)]['KEGG'].tolist()
    background_KEGG = mat['KEGG'].dropna().tolist()

    # Return KEGG matrix
    KEGG_dict = dict(zip(mat.index.tolist(), mat['KEGG'].tolist()))
    mat_KEGG = mat_proc.rename(columns=KEGG_dict)
    mat_KEGG = mat_KEGG.loc[:, mat_KEGG.columns.notnull()]

    return DEM_KEGG_id, background_KEGG, mat_KEGG

def stevens_data():
    md_raw = pd.read_csv("../Stevens/MTBLS136_compressed_files/s_MTBLS136.txt", sep="\t")
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
    mat = pd.read_csv("Stevens_matrix_named_compounds_only.csv", index_col=0)
    mat_nonusers_estrogen = mat.drop((replicate_samples + estrogen_progesterone), axis=1)
    stevens_matrix_proc = utils.data_processing(mat_nonusers_estrogen.T, 8, 0)
    stevens_matrix_proc["Group"] = stevens_matrix_proc.index.map(sample_status_dict)
    ttest_res = utils.t_tests(stevens_matrix_proc.iloc[:,:-1], stevens_matrix_proc["Group"], "fdr_bh")
    DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()

    # Add back the metadata columns
    mat.columns = mat.columns.str.replace(' ', '')
    metadata_cols = ['KEGG', 'SampleHMDB_ID']
    stevens_matrix_proc_annotated = stevens_matrix_proc.T.join(mat[metadata_cols])
    background_list = stevens_matrix_proc.columns.tolist()
    DEM_KEGG_id = mat[mat.index.isin(DEM)]['KEGG'].tolist()
    DEM_KEGG_id = [i for i in DEM_KEGG_id if str(i) != 'nan']
    stevens_background_list = stevens_matrix_proc_annotated['KEGG'].dropna().tolist()

    # Return KEGG matrix
    KEGG_dict = dict(zip(mat.index.tolist(), mat['KEGG'].tolist()))
    mat_KEGG = stevens_matrix_proc.rename(columns=KEGG_dict)
    mat_KEGG = mat_KEGG.loc[:, mat_KEGG.columns.notnull()]

    return DEM_KEGG_id, stevens_background_list, mat_KEGG

def zamboni_data(knockout):
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

    background_list_all_annotations = list(set(sum(annotations_neg.values(), []) + sum(annotations_pos.values(), [])))

    DEM = list(set(strain_DA_compounds[knockout]))

    return DEM, background_list_all_annotations