# Get latest KEGG pathways using rest API

import requests
import re
import pandas as pd

print("Beginning KEGG download...")
# get all pathways
url = 'http://rest.kegg.jp/list/pathway/mmu'
# change organism name
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
base_url = 'http://rest.kegg.jp/get/'
#eco00010

pathway_ids = [*pathway_dict]
pathway_names = list(pathway_dict.values())
pathway_compound_mapping = dict()
for i in pathway_ids:
    print(i)
    current_url = base_url + i + "/conf"
    page = requests.get(current_url)
    lines = page.text.split("\n")
    complist = []
    for l in lines:
        if l.startswith("circ"):
            l = l.split()
            complist.append(l[4])
    pathway_compound_mapping[i] = complist

df = pd.DataFrame.from_dict(pathway_compound_mapping, orient='index')
df.insert(0, 'Pathway_name', pathway_names)
df.to_csv("KEGG_mouse_pathways_compounds.csv")
print("Complete!")



