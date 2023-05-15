###############################################################################
# This script execute an automated search of all possible combinations of two #
# sets of terms. The search is done in pubmed. The script retrieves all pubmed#
# ids and store them in text files. The search can be updated by changing the #
# "term" variable.                                                            #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
# Date: 29.10.2020                                                            #
###############################################################################

import json
import os
import pandas as pd
import numpy as np
from Bio import Entrez
       
#Organize the terms
Entrez.email = 'your_email@provider.com'
path = "/results/pubmeddata/"
first_terms = pd.read_csv("/src/keys1.csv",header=None)
second_terms = pd.read_csv("/src/keys2.csv",header=None)
first_terms[0] = first_terms[0].str.replace("+", " ")
second_terms[0] = second_terms[0].str.replace("+", " ")
first_terms[0] = first_terms[0].str.replace("%22", '')
second_terms[0] = second_terms[0].str.replace("%22", '')
first_terms['key'] = 1
second_terms['key'] = 1
combi = pd.merge(first_terms,second_terms,on='key').drop('key',axis=1)

#Request results and save them      
for i,row in combi.iterrows():
    term1 = row['0_x']
    term2 = row['0_y']
    title = path+term1.replace(' ','_')+"-"+term2.replace(' ','_')+".txt" 
    term = ('((((("%s"[Title/Abstract]) AND ("%s"[Title/Abstract])) AND ("english"[Language])) AND (("journal article"[Publication Type]) OR ("congress"[Publication Type]) OR ("clinical study"[Publication Type]))))' % (term1, term2))
    handle = Entrez.esearch(db="pubmed", term=term, RetMax='100000')
    record = Entrez.read(handle)
    pubmedids = record["IdList"]
    with open(title, 'w') as f:
        for item in pubmedids:
            f.write("%s\n" % item)
    print(term)
