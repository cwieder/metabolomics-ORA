import pandas as pd

infile = pd.read_csv("../ChEBI2Reactome.txt", sep="\t", header=None)
infile = infile.astype(str)
name_dict = dict(zip(infile[1], infile[3]))

groups = infile.groupby([1])[0].apply(list).to_dict()
df = pd.DataFrame.from_dict(groups, orient='index')
pathway_ids = df.index.tolist()

df['Pathway_name'] = df.index.map(name_dict)
col = df.pop("Pathway_name")
df.insert(0, "Pathway_name", col)
print(df.head)
# df.to_csv("Reactome_pathway_set.csv")