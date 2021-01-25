import pickle

with open('HMDB_logP_dict_all_metabolites.pickle', 'rb') as handle:
    logp_dict = pickle.load(handle)

print(len(logp_dict))
print(logp_dict)