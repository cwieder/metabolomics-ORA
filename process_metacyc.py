import pandas as pd
from itertools import groupby

raw = pd.read_csv("All_pathways_of_H._sapiens.txt", sep="\t", index_col=0)

print(raw.columns)
pathway_compounds = dict.fromkeys(raw.index.tolist())

for i in raw.itertuples():
    cpds = i[2].split(" // ")
    pathway_compounds[i[0]] = cpds

# print(len(pathway_compounds))

alt_compounds = {}
with open("../tier1/24.1/data/compounds.dat", "r", encoding='latin-1', newline="\n") as infile:
    # infile = infile.read().splitlines()
    all_lines = []
    for line in infile:
        line = line.strip()
        if not line.startswith("#"):
            all_lines.append(line)
    split_at = "//"
    all_entries = [list(g) for k, g in groupby(all_lines, lambda x: x != split_at) if k]
    for i in all_entries:
        cpd = str(i[0]).replace("UNIQUE-ID - ", "")
        alt_compounds[cpd] = []
        for field in i:
            if field.startswith("DBLINKS"):
                field = field.split()
                alt_compounds[cpd].append(field[3].replace('"', ''))

# print(pathway_compounds)
pathway_df = pd.DataFrame.from_dict(pathway_compounds, orient="index")
print(pathway_df)