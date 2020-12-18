import requests
import re
import pandas as pd

print("Beginning KEGG download...")
# get all pathways
url = 'http://rest.kegg.jp/list/compound'
# change organism name
data = requests.get(url)
compounds = data.text
compounds = compounds.split("\n")
compound_dict = dict()

for entry in compounds:
    if entry:
        cpd = entry.split("\t")
        cpd_id = cpd[0].replace("cpd:", "")
        compound_dict[cpd_id] = ""

for entry in compound_dict.keys():
    cpd_url = 'http://rest.kegg.jp/get/' + entry
    data = requests.get(cpd_url)
    cpd_info = data.text
    cpd_info = cpd_info.split("\n")
    for i in cpd_info:
        if i.startswith("FORMULA"):
            formula = i.split()
            compound_dict[entry] = formula[1]

df = pd.DataFrame.from_dict(compound_dict, orient="index", columns=["molecular_formula"])
df.to_csv("KEGG_compound_formulae.csv")
print("Completed KEGG formula download!")


