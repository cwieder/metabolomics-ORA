# Get latest KEGG pathways using rest API

import requests
import re
import pandas as pd

print("Beginning KEGG download...")
# get all pathways
url = 'http://rest.kegg.jp/list/pathway'
data = requests.get(url)
pathways = data.text
pathways = pathways.split("\n")
pathways = filter(None, pathways)
pathway_dict = dict()

for path in pathways:
    path = path.split("\t")
    name = path[1]
    pathid = re.search(r"path:(.*)", path[0]).group(1)
    pathway_dict[pathid] = name

# get compounds for each pathway
base_url = 'http://rest.kegg.jp/link/cpd/'
pathway_ids = [*pathway_dict]
pathway_names = list(pathway_dict.values())

pathway_compound_mapping = dict()
for i in pathway_ids:
    current_url = base_url + i
    compounds = requests.get(current_url)
    c_list = compounds.text.split("\n")
    c_list = filter(None, c_list)
    complist = []
    for comp in c_list:
        comp = comp.split("\t")
        comp_id = re.search(r"cpd:(.*)", comp[1]).group(1)
        complist.append(comp_id)
    pathway_compound_mapping[i] = complist


df = pd.DataFrame.from_dict(pathway_compound_mapping, orient='index')
df.insert(0, 'Pathway_name', pathway_names)

df.to_csv("KEGG_reference_pathways_compounds.csv")
print("Complete!")



