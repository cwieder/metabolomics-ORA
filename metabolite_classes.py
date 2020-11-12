# Get metabolite classes from HMDB

from getDB import hmdb
import pandas as pd
import numpy as np

hmdb_list = pd.read_csv("Stevens_HMDB_list.txt", header=None).iloc[:,0].tolist()
metabo_class = []

print("Starting...")
for item in hmdb_list:
    hmdb_id = item.split("HMDB")
    print(item)

    temp = hmdb.getMetabolite(ID = "HMDB" + "00" + hmdb_id[1])
    if temp != None:
        ##check the contents of temp
        taxonomy = temp["taxonomy"]
        if not isinstance(taxonomy, str):
            taxonomy = taxonomy.rename(columns=taxonomy.iloc[0]).drop(taxonomy.index[0])
            # print(temp["name"], taxonomy["class_new"].values)
            metabo_class.append(taxonomy["class_new"].values[0])
        else:
            metabo_class.append("Unknown")
    else:
        print(item, "Not found")
        metabo_class.append("Unknown")

results = pd.DataFrame(zip(hmdb_list, metabo_class), columns=["HMDB", "HMDB_class"])
results.to_csv("Stevens_metabolite_classes.csv")
print("Complete!")

# temp = hmdb.getMetabolite(ID = "HMDB0000327")

