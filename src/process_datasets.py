import pandas as pd
import numpy as np
import utils
import re
from bioservices import *
import pickle

class Dataset:
    """
    Class for databset processing
    """
    def __init__(self, name, dbname="KEGG"):
        self.name = name
        self.dbname = dbname
        self.dem = []
        self.background = []
        self.proc_mat = []

        # Each datasets requires a unique processing method due to dataset differences
        if not self.name.startswith(("Yachida", "Labbe", "Stevens", "Fuhrer", "Quiros")):
            raise NotImplementedError(
                "Dataset not supported. Add to this class to include and pre-process custom datasets.")
        if self.dbname not in ["KEGG", "reactome", "biocyc"]:
            raise ValueError("Database not supported.")

        # Process individual datasets
        if self.name == "Yachida":
            self.process_yachida()
        elif self.name == "Labbe":
            self.process_labbe()
        elif self.name == "Stevens":
            self.process_stevens()
        elif self.name.startswith("Fuhrer"):
            self.process_fuhrer()
        elif self.name == "Quiros":
            self.process_quiros()

        if not self.name.startswith("Fuhrer"):
            # Change pathway database
            if dbname != "KEGG":
                self.pathway_db()
            # Get DA metabolites and background list
            self.get_dem()

    def get_dem(self):
        """
        Runs two-sided t-tests to get differentially abundant metabolites.
        Also assigns background list.
        """
        ttest_res = utils.t_tests(self.proc_mat.iloc[:, :-1], self.proc_mat["Group"], "fdr_bh")
        DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
        background = self.proc_mat.iloc[:, :-1].columns.tolist()
        self.dem = DEM
        self.background = background

    def pathway_db(self):
        """
        Change pathway database from KEGG (default) to BioCyc or Reactome
        :return: new processed data matrix with pathway database-specific IDs
        """

        kegg_db = KEGG(verbose=False)

        if self.dbname == "KEGG":
            pass

        elif self.dbname == "reactome":
            map_kegg_chebi = kegg_db.conv("chebi", "compound")
            KEGG2Reactome = dict.fromkeys(self.proc_mat.columns.tolist()[:-1])
            for cpd in KEGG2Reactome.keys():
                try:
                    chebiID = map_kegg_chebi['cpd:' + cpd]
                    KEGG2Reactome[cpd] = chebiID[6:]
                except KeyError:
                    KEGG2Reactome[cpd] = np.nan
            data_proc = self.proc_mat.rename(columns=KEGG2Reactome)
            data_proc = data_proc.loc[:, data_proc.columns.notnull()]
            self.proc_mat = data_proc

        elif self.dbname == "biocyc":
            biocyc_mapping = {"Yachida": "../yamada2metacyc.txt",
                              "Labbe": "../brown2metacyc.txt",
                              "Stevens": "../stevens2metacyc.txt",
                              "Quiros": "../auwerx2metacyc.txt",
                              "Fuhrer_dcuS": "../zamboni2metacyc.txt",
                              "Fuhrer_yfgM": "../zamboni2metacyc.txt"}
            mapping = pd.read_csv(biocyc_mapping[self.name], sep="\t")
            KEGG2BioCyc = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
            data_proc = self.proc_mat.rename(columns=KEGG2BioCyc)
            data_proc = data_proc[data_proc.columns.dropna()]
            self.proc_mat = data_proc

        # Implement other databases here

    def process_yachida(self):
        data = pd.read_csv("../example_data/yachida_abundance.csv", index_col=0, header=0).T
        data = data.rename(columns={'Group': 'disease'})
        data = data.dropna(0)
        sample_disease_dict = dict(zip(data.index, data['disease']))
        data.columns = data.columns[0:4].tolist() + [col[0:6] for col in data.columns[4:]]

        removecols = []
        for i in data.columns.tolist():
            matchObj = re.search("^[C]\d{5}$", i)
            if not matchObj:
                removecols.append(i)

        data = data.drop(removecols[4:], axis=1)

        CRC_or_healthy_dict = dict.fromkeys(data.index.tolist())
        for k, v in sample_disease_dict.items():
            if v in ['Healthy']:
                CRC_or_healthy_dict[k] = "Healthy"
            elif v in ["Stage_I_II", "Stage_III_IV"]:
                CRC_or_healthy_dict[k] = "CRC"
            else:
                CRC_or_healthy_dict[k] = "Null"
        CRC_or_healthy = ["Healthy" if i in ["Healthy"] else "CRC" if i in ["Stage_I_II", "Stage_III_IV"] else "Null"
                          for i in data["disease"]]

        data.insert(1, "Group", CRC_or_healthy)

        data = data.iloc[:, ~data.columns.duplicated()]
        df = data[data.disease != "Null"]
        data_proc = utils.data_processing(df, 0, 5)
        data_proc["Group"] = data_proc.index.map(CRC_or_healthy_dict)
        data_proc = data_proc[data_proc.Group != "Null"]
        self.proc_mat = data_proc

    def process_labbe(self):
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
        self.proc_mat = mat_proc

    def process_stevens(self):
        md_raw = pd.read_csv("../example_data/Stevens_metadata.txt", sep="\t")
        metadata_list = list(zip(md_raw['Factor Value[CurrentPMH]'], md_raw['Factor Value[Gender]'],
                                 md_raw['Factor Value[AgeAtBloodDraw]'],
                                 ['Over 75' if val not in ['<=55', '56-60', '61-65', '66-70', '71-75'] else 'Under 75'
                                  for val in md_raw['Factor Value[AgeAtBloodDraw]']]))
        metadata_dict = dict(zip(md_raw['Sample Name'].values, metadata_list))
        sample_status_dict = dict(zip(md_raw['Sample Name'].values, md_raw['Factor Value[CurrentPMH]']))

        replicate_samples = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', 'E+P']]
        nonusers = [k for k, v in metadata_dict.items() if v[0] not in [np.nan, 'E-only', 'E+P']]
        estrogen_only = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', np.nan, 'E+P']]
        estrogen_progesterone = [k for k, v in metadata_dict.items() if v[0] not in ['Nonuser', 'E-only', np.nan]]

        # Get abundance matrix, transpose to n-samples by m-metabolites
        mat = pd.read_csv("../example_data/Stevens_matrix_named_compounds_only.csv", index_col=0, dtype=object)
        mat_nonusers_estrogen = mat.drop((replicate_samples + estrogen_progesterone), axis=1)
        stevens_matrix_proc = utils.data_processing(mat_nonusers_estrogen.T, 8, 0)
        stevens_matrix_proc["Group"] = stevens_matrix_proc.index.map(sample_status_dict)
        mapping_dict = dict(zip(mat.index, mat['KEGG']))
        stevens_matrix_proc = stevens_matrix_proc.rename(columns=mapping_dict)
        stevens_matrix_proc = stevens_matrix_proc.loc[:, stevens_matrix_proc.columns.notnull()]
        stevens_matrix_proc = stevens_matrix_proc.loc[:, ~stevens_matrix_proc.columns.duplicated()]
        self.proc_mat = stevens_matrix_proc

    def process_fuhrer(self):
        knockout = self.name[7:]
        # import modified z-scores
        n_zscore = pd.read_csv("../example_data/Fuhrer_mod_zscore_neg_CW.csv.zip", index_col=0)
        p_zscore = pd.read_csv("../example_data/Fuhrer_mod_zscore_pos_CW.csv.zip", index_col=0)

        # remove unannotated
        n_zscore = n_zscore[n_zscore.index.notnull()]
        p_zscore = p_zscore[p_zscore.index.notnull()]

        with open('../data/zamboni_pos_annotation_dict2.pickle', 'rb') as handle:
            annotations_pos = pickle.load(handle)

        with open('../data/zamboni_neg_annotation_dict2.pickle', 'rb') as handle:
            annotations_neg = pickle.load(handle)

        # convert to CHEBI for Reactome
        if self.dbname == "reactome":
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

        if self.dbname == "biocyc":
            mapping = pd.read_csv("../zamboni2metacyc.txt", sep="\t")
            mapping_dict = dict(zip(mapping["Kegg"].tolist(), mapping["BioCyc"].tolist()))
            annotations_neg = {k: [mapping_dict[i] if i in mapping_dict.keys() else "" for i in v] for k, v in
                               annotations_neg.items()}
            annotations_pos = {k: [mapping_dict[i] if i in mapping_dict.keys() else "" for i in v] for k, v in
                               annotations_pos.items()}

        strain_DA_compounds = dict.fromkeys(n_zscore.columns)

        for strain in strain_DA_compounds.keys():
            cur_col = n_zscore.loc[:, strain]
            DA_metabolites = []
            for items in cur_col.iteritems():
                if items[1] > 6 or items[1] < -6:
                    DA_metabolites.append(annotations_neg[items[0]])
            DA_KEGG = [j for i in DA_metabolites for j in i]  # filter out only annotated compounds
            strain_DA_compounds[strain] = DA_KEGG

        for strain in strain_DA_compounds.keys():
            cur_col = p_zscore.loc[:, strain]
            DA_metabolites = []
            for items in cur_col.iteritems():
                if items[1] > 6 or items[1] < -6:
                    DA_metabolites.append(annotations_pos[items[0]])
            DA_KEGG = [j for i in DA_metabolites for j in i]  # filter out only annotated compounds
            strain_DA_compounds[strain] += DA_KEGG
        background_list_all_annotations = list(
            filter(None, list(set(sum(annotations_neg.values(), []) + sum(annotations_pos.values(), [])))))
        DEM = list(filter(None, list(set(strain_DA_compounds[knockout]))))

        # return strain matrix with all annotations
        empty_anno_pos = [k for k, v in annotations_pos.items() if not v]

        strain_pos_all = p_zscore.drop(empty_anno_pos).loc[:, knockout]
        pos_rows = []
        for ion, val in strain_pos_all.items():
            ion_annos = annotations_pos[ion]
            for i in ion_annos:
                row = ("pos" + str(ion), i, val)
                pos_rows.append(row)
        pos_annos_df = pd.DataFrame(pos_rows, columns=['Ion', 'Annotation', 'Z-score'])

        empty_anno_neg = [k for k, v in annotations_neg.items() if not v]
        strain_neg_all = n_zscore.drop(empty_anno_neg).loc[:, knockout]
        neg_rows = []
        for ion, val in strain_neg_all.items():
            ion_annos = annotations_neg[ion]
            for i in ion_annos:
                row = ("neg" + str(ion), i, val)
                neg_rows.append(row)
        neg_annos_df = pd.DataFrame(neg_rows, columns=['Ion', 'Annotation', 'Z-score'])
        all_annotations_df = pd.concat([pos_annos_df, neg_annos_df])
        all_annotations_df = all_annotations_df.set_index(all_annotations_df["Annotation"])
        mat = all_annotations_df.iloc[:, 2:].T

        self.proc_mat = mat
        self.dem = DEM
        self.background = background_list_all_annotations


    def process_quiros(self):
        mat = pd.read_csv("../example_data/quiros_abundance.csv", index_col=6).T
        mat = mat.iloc[6:, :]
        mat = mat.loc[:, ~mat.columns.duplicated(keep='first')]
        mat = mat.loc[:, mat.columns.notnull()]
        groups = [i.split(".", 1)[0] for i in mat.index.tolist()]
        mat['Group'] = groups
        mat_selected_groups = mat.loc[mat['Group'].isin(['Acti', 'Dox'])]  # Acti vs Dox
        matrix_proc = utils.data_processing(mat_selected_groups.iloc[:, :-1], 0, 0)
        matrix_proc_copy = matrix_proc.copy()
        matrix_proc_copy['Group'] = mat_selected_groups['Group']
        self.proc_mat = matrix_proc_copy

fuhrer_dcus_data_r = Dataset("Fuhrer_dcuS", dbname='reactome')
print(fuhrer_dcus_data_r.dem)
