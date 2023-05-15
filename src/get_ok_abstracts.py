###############################################################################
# This script This function extracts the pubmed id of the papers that are     #
# selected based on abstracts reviews.                                        #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
# Date: 30.08.2021: updat                                                     #
###############################################################################

import pandas as pd

###                         Start with Pubmed                               ###
df = pd.read_excel('../results/Literature Review.xlsx', sheet_name="PubMed")
abs_ok = df.loc[df['Agreement']==1]
abs_ok_list = set(abs_ok["pubmed_id"].tolist())

#Check those pubmed that have been checked
done = pd.read_excel('../results/Literature Review.xlsx', sheet_name="Full-text include")
done_list = set(done["pubmed_id"].tolist())

if not done_list.issubset(abs_ok_list):
    print("WARNING: mismatch between pubmed ids. The abstracts that have been checked are not a subset of the abstracts to include.")

to_check = pd.concat([abs_ok, done])
to_check.drop_duplicates(subset="pubmed_id", keep=False, inplace=True)

to_check.to_excel('./results/Abstract_FullText.xlsx',sheet_name="PubMed",index=False)


###                         Continue with Scopus                            ###
df = pd.read_excel('./results/Literature Review.xlsx', sheet_name="Scopus")
abs_ok = df.loc[df['Agreement']==1]
abs_ok_list = set(abs_ok["scopus_id"].tolist())

abs_ok.to_excel('./results/Abstract_FullText.xlsx',sheet_name="Scopus",index=False)
