
import requests
import re
import pandas as pd
import pickle

url = 'http://rest.kegg.jp/list/compound'
data = requests.get(url)
compounds = data.text

compounds = compounds.split("\n")
compounds = filter(None, compounds)

compound_dict = dict()

for compound in compounds:
    compound = compound.split("\t")
    cpd_id = re.search(r"cpd:(.*)", compound[0]).group(1)
    cpd_names = compound[1].split(";")
    compound_dict[cpd_id] = cpd_names

print(compound_dict)
with open('KEGG_compounds.pickle', 'wb') as handle:
    pickle.dump(compound_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("Done.")