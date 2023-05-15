###############################################################################
# This script reads pubmed ids and retrieves the title, abstract, and publica #
# tion date for each of them. The script also removes duplicates and saves the#
# results in an excel file.                                                   #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
# Date: 29.10.2020                                                            #
###############################################################################

from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import os

Entrez.email = 'your_email@provider.com'
path = '/results/pubmeddata/'
savepath = '/results/literature.xlsx'

files = [os.path.join(path,x) for x in os.listdir(path)]

pmids=[]
for each_file in files:
    with open(each_file) as f:
        lines = f.readlines()
        lines = [x.replace('\n', '') for x in lines]
        pmids.extend(lines)
print('Files read')
  
handle = Entrez.efetch(db="pubmed", id=','.join(map(str, pmids)),
                       rettype="xml", retmode="text")
records = Entrez.read(handle)
print('handle done')

abstract_dict = {}
without_abstract = {}

for pubmed_article in records['PubmedArticle']:
    pmid = int(str(pubmed_article['MedlineCitation']['PMID']))
    article = pubmed_article['MedlineCitation']['Article']
    if 'Abstract' in article:
        abstract = article['Abstract']['AbstractText'][0]
        title = article['ArticleTitle']
        if article['ArticleDate']:
            date = article['ArticleDate'][0]
            date = '.'.join(str(x) for x in date.values())
        else:
            date = ' '
        abstract_dict[pmid] = [abstract, title, date]
    else:
        noabstract = 'no abstract'
        noabs_title = article['ArticleTitle']
        if article['ArticleDate']:
            noabs_date = article['ArticleDate'][0]
            noabs_date = '.'.join(str(x) for x in noabs_date.values())
        else:
            noabs_date = 'no date'
        without_abstract[pmid] = [noabstract, noabs_title, noabs_date]
       
print('abstracts read')
       
df = pd.DataFrame(data=abstract_dict)
df = df.drop_duplicates()
df = (df.T)
df = df.rename(columns={0: "abstract", 1: "title", 2:"date"})
df.index.names = ['pubmed_id']
df['date'] = df['date'].str.replace('.', '-', regex=True)
df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d", errors='coerce')
df = df[(df['date'] > '2000-01-01')] #Filter in the last 20 years
df['date'] = df['date'].dt.strftime('%Y-%m-%d')

no_df = pd.DataFrame(data=without_abstract)
no_df = no_df.drop_duplicates()
no_df = (no_df.T)
no_df = no_df.rename(columns={0: "abstract", 1: "title", 2:"date"})
no_df.index.names = ['pubmed_id']
if len(no_df)>0:
    no_df['date'] = no_df['date'].str.replace('.', '-', regex=True)
    no_df['date'] = pd.to_datetime(no_df['date'], format="%Y-%m-%d", errors='coerce')
    no_df = no_df[(no_df['date'] > '2000-01-01')] #Filter in the last 20 years
    no_df['date'] = no_df['date'].dt.strftime('%Y-%m-%d')

print('writing in excel')

with pd.ExcelWriter(savepath) as writer:  
    df.to_excel(writer, sheet_name='Abstracts')
    no_df.to_excel(writer, sheet_name='No_Abstracts')
    
print('done :)')
