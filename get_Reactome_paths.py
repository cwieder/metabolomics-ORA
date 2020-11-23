import pandas as pd

infile = pd.read_csv("../ChEBI2Reactome.txt", sep="\t", header=None)
infile = infile.astype(str)



groups = infile.groupby([1])[0].apply(list).to_dict()
df = pd.DataFrame.from_dict(groups, orient='index')
pathway_ids = df.index.tolist()
# pathway_names = infile[infile[1].isin(groups[1])][3].tolist()

df.insert(loc=0, column='Pathway_name', value=['' for i in range(df.shape[0])])

df.to_csv("Reactome_pathway_set.csv")