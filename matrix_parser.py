import numpy as np
import pandas as pd

def parse_matrix(filename):
    # File formats to support csv or excel
    print("Reading from:", filename)
    mat = pd.read_excel(filename)
    print(mat.head())
    return mat

peak_area_matrix = parse_matrix("../Polonis_EVA_Data/EVA_aligned filtered data_positive_negative.xls")

# peak_area_matrix = pd.concat([peak_area_matrix,pd.DataFrame(columns=['Name', 'KEGG', 'HMDB','LipidMaps', 'MetLin', 'InChIKey', 'SMILES'])])
peak_area_matrix = peak_area_matrix.reindex(peak_area_matrix.columns.tolist() + ['Name', 'KEGG', 'HMDB','LipidMaps', 'MetLin', 'InChIKey', 'SMILES'], axis=1)

# Define sample classes
putative_characterisation = pd.read_excel("../Polonis_EVA_Data/EVA_MS_putative characetrization.xls")
print(putative_characterisation.head())

# Match up compound names to mass/RT
mass = peak_area_matrix["Mass"].tolist()
retention_time = peak_area_matrix["RT"].tolist()
mass_RT = tuple(zip(mass, retention_time))


for i, pair in enumerate(mass_RT):
    for idx, row in putative_characterisation.iterrows():
        if pair[0] == row[0] and pair[1] == row[1]:
            peak_area_matrix.loc[i,"Name"] = row[2]
            peak_area_matrix.loc[i, "KEGG"] = row[5]
            peak_area_matrix.loc[i, "HMDB"] = row[6]
            peak_area_matrix.loc[i, "LipidMaps"] = row[7]
            peak_area_matrix.loc[i, "MetLin"] = row[8]
            peak_area_matrix.loc[i, "InChIKey"] = row[9]
            peak_area_matrix.loc[i, "SMILES"] = row[10]

peak_area_matrix.to_csv("Polonis_peak_area_matrix_with_annotations.csv")