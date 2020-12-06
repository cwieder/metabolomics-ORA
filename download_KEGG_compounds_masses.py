import requests
import re
import pandas as pd

print("Beginning KEGG download...")
# get all compounds
dfs = []
intervals = ['0-100', '100-200', '200-300', '300-400', '400-500',
             '500-600', '600-700', '700-800', '800-900', '900-1000', '1000-1100', '1100-5000']
for i in intervals:
    url = 'http://rest.kegg.jp/find/compound/'+ str(i) +'/exact_mass'
    # change organism name
    data = requests.get(url)
    compounds = data.text
    compounds = compounds.split("\n")
    compound_dict = {}

    for entry in compounds:
        if entry:
            cpd = entry.split("\t")
            cpd_id = cpd[0].replace("cpd:", "")
            compound_dict[cpd_id] = cpd[1]

    df = pd.DataFrame.from_dict(compound_dict, orient="index", columns=["molecular_weight"])
    dfs.append(df)
compounds_masses = pd.concat(dfs)
print(compounds_masses.shape)
compounds_masses.to_csv("KEGG_compounds_exact_mass.csv")
