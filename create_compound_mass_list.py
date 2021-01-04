import pandas as pd
from molmass import Formula
import numpy as np

kegg_cpd = pd.read_csv("KEGG_compound_formulae.csv", index_col=0)
kegg_mass = pd.read_csv("KEGG_compounds_exact_mass.csv", index_col=0)

kegg_cpd['mass'] = kegg_mass[kegg_mass.index.isin(kegg_cpd.index.tolist())]["molecular_weight"]

print(kegg_cpd['mass'].isna().sum())
mapping_dict = {}
for row in kegg_cpd.itertuples():
    formula = str(row[1])
    mass = str(row[2])
    if formula != "nan" and mass == "nan":
        try:
            f = Formula(formula.replace("n", "").replace("R", ""))
            mapping_dict[formula] = f.mass
        except:
            print(formula)


kegg_cpd['mass'] = kegg_cpd['molecular_formula'].map(mapping_dict).fillna(kegg_cpd['mass'])
print(kegg_cpd['mass'].isna().sum())

kegg_cpd.to_csv("KEGG_compounds_masses_estimated.csv")