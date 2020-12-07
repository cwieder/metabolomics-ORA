import pandas as pd
import numpy as np

anno_df = pd.read_excel("../Zamboni/annotations_EV1.xlsx", sheet_name="Table EV1B")
anno_df["anno_split"] = anno_df['Annotations (Compound Name, Modification, Sum-Formula, KEGG Code, HMDB Code)'].str.split(' ')

annotations = (anno_df["anno_split"])
kegg_id = [i[-2] if i is not np.nan else 'nan' for i in annotations]

anno_df["KEGG"] = kegg_id

anno_df.to_csv("../Zamboni/Zamboni_KEGG_annotation.csv", index=False)
