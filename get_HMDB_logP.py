from bs4 import BeautifulSoup as bs
import re
import pickle

content = []
# Read the XML file
# with open("csf_metabolites.xml", "r") as file:
with open("/rds/general/user/cw2019/home/metabo_pathways/hmdb_metabolites.xml", "r") as file:
    # Read each line in the file, readlines() returns a list of lines
    content = file.readlines()
    # Combine the lines in the list into a string
    content = "".join(content)
    bs_content = bs(content, "lxml")

logP_dict = {}

all_metabolites = bs_content.find_all("metabolite")
print(len(all_metabolites))
for metabolite in all_metabolites:
    metabolite = str(metabolite)
    accessions = []
    logp = []
    kegg_id = ""
    metabolite_list = metabolite.split(">\n<")
    print(len(metabolite_list))
    for num, entry in enumerate(metabolite_list):
        if entry.startswith("accession"):
            acc = re.search("(?<=>).*?(?=<)", entry)
            accessions.append(acc.group(0))
        if entry == "kind>logp</kind":
            lp = re.search("(?<=>).*?(?=<)", metabolite_list[num+1])
            logp.append(lp.group(0))
        if entry.startswith("kegg_id"):
            ki = re.search("(?<=>).*?(?=<)", entry)
            kegg_id = ki.group(0)
    logP_dict[accessions[0]] = {"accessions": accessions, "kegg_ID": kegg_id, "logP": logp}

print(len(logP_dict))

with open('HMDB_logP_dict_all_metabolites.pickle', 'wb') as handle:
    pickle.dump(logP_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)