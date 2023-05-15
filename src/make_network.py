###############################################################################
# This script converts the table to a network and prepares the files to the   #
# gephi format                                                                #
#                                                                             #
# 10.05.2022: reated by ana.trianahoyos@aalto.fi                              #
###############################################################################

import pandas as pd

savepath = './results'

df = pd.read_excel("./results/results_organized.xlsx")

cols = ['blood_pressure','comm_patterns','EMA','glucose','heart_rate','hrv','mobility_patterns','images',
        'physical_activity_sensor','resp_rate_sensor','screen_use','sleep_sensor','temperature','app_use',
        'dti', 'smri', 'fmri', 'mr_angiography']
df = df[cols]
df.rename(columns={"blood_pressure": "blood pressure", "comm_patterns": "communication patterns",
                     "EMA": "EMA", "glucose":"glucose", "heart_rate": "heart rate", "hrv": "heart rate variability",
                     "mobility_patterns":"mobility patterns", "images":"images", 
                     "physical_activity_sensor":"physical activity", "resp_rate_sensor":"respiration rate",
                     "screen_use":"screen use", "sleep_sensor":"sleep", "temperature":"temperature",
                     "app_use":"app use", "dti":"DWI", "smri":"T1/T2-weighted MRI", "fmri":"fMRI", 
                     "mr_angiography":"angiography"}, inplace=True)

graph = df.T.dot(df)
graph.to_csv(f'{savepath}/graph.csv')

#create an edge list too
nodes = graph.columns

# Open a file to write the edge list
with open(f'{savepath}/network.txt', 'w') as f:
    # Iterate through the upper triangle of the adjacency matrix
    for i in range(graph.shape[0]):
        for j in range(i+1, graph.shape[1]):
            node_i = nodes[i]
            node_j = nodes[j]
            value = graph.iloc[i, j]

            if value > 0:
                f.write(f"{node_i}\t{node_j}\t{value}\n")
