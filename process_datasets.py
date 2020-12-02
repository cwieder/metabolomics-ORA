import pandas as pd
import numpy as np
import utils

# IMPORT DATASETS, PRE-PROCESS THEM AND RUN T-TESTS TO OBTAIN LIST OF DIFFERENTIALLY ABUNDANT METABOLITES

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
    background_KEGG = mat['KEGG'].tolist()
    return DEM_KEGG_id, background_KEGG

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

    return DEM_KEGG_id, stevens_background_list