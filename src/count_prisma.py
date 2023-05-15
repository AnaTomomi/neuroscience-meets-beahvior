##############################################################################
#This file makes the count for generating the PRISMA numbers automatically   #
#                                                                            #
#Created: 29.03.2023                                                         #
#Author: ana.trianahoyos@aalto.fi                                            #
##############################################################################

from Bio import Entrez
import pandas as pd
import json
import xml.etree.ElementTree as ET
import os
import numpy as np

Entrez.email = 'your_email@provider.com'
pubmed_update_path = './results/pubmeddata/'
scopus_update_path = './results/scopus/'

master = './results/Literature Review.xlsx'

###############################################################################
#Utility functions do not touch
def organize_pmids(path):
    files = [os.path.join(path,x) for x in os.listdir(path)]

    pmids=[]
    for each_file in files:
        with open(each_file) as f:
            lines = f.readlines()
            lines = [x.replace('\n', '') for x in lines]
            pmids.extend(lines)
    return pmids

def retrieved_pmids(pmids):
    handle = Entrez.efetch(db="pubmed", id=','.join(map(str, pmids)),
                       rettype="xml", retmode="text")
    records = Entrez.read(handle)
    
    abstract_dict={}
    count=0
    for pubmed_article in records['PubmedArticle']:
        pmid = int(str(pubmed_article['MedlineCitation']['PMID']))
        count = count+1
        article = pubmed_article['MedlineCitation']['Article']
        if article['ArticleDate']:
            date = article['ArticleDate'][0]
            date = '.'.join(str(x) for x in date.values())
        else:
            date = ' '
        abstract_dict[count] = [date, pmid]
    df = pd.DataFrame(data=abstract_dict)
    df = (df.T)
    df = df.rename(columns={0: "date", 1: "pmid"})
    df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d", errors='coerce')
    df = df[(df['date'] > '2000-01-01')] #Filter in the last 20 years
    return df

def discard_pmids_other(pmids):
    handle = Entrez.efetch(db="pubmed", id=','.join(map(str, pmids)),
                       rettype="xml", retmode="text")
    records = Entrez.read(handle)
    
    abstract_dict={}
    count=0
    for pubmed_article in records['PubmedBookArticle']:
        pmid = pubmed_article['BookDocument']['PMID']
        date = pubmed_article['BookDocument']['Book']['PubDate']
        date = '.'.join(str(x) for x in date.values())
        abstract_dict[count] = [date, pmid]
    df = pd.DataFrame(data=abstract_dict)
    df = (df.T)
    df = df.rename(columns={0: "date", 1: "pmid"})
    df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d", errors='coerce')
    df = df[(df['date'] > '2000-01-01')] #Filter in the last 20 years
    return df   

def organize_scopusids(path):
    files = [os.path.join(path,x) for x in os.listdir(path)]
    scopusids = []
    for i,file in enumerate(files):
        with open(file) as json_file:
            data = json.load(json_file)
        if int(data['search-results']['opensearch:totalResults'])>0: #does it have results?
            for j,item in enumerate(data['search-results']['entry']):
                 scopusids.append(item['dc:identifier'].strip('SCOPUS_ID:'))
    return scopusids
                 
def retrieved_scopus(path):
    files = [os.path.join(path,x) for x in os.listdir(path)]
    abstract_dict = {}
    count=0
    for i,file in enumerate(files):
        with open(file) as json_file:
            data = json.load(json_file)
        if int(data['search-results']['opensearch:totalResults'])>0: #does it have results?
            for j,item in enumerate(data['search-results']['entry']):
                 count = count+1
                 scopus_id = item['dc:identifier'].strip('SCOPUS_ID:')
                 if 'prism:coverDate' in item.keys():
                     date = item['prism:coverDate']
                 else:
                    date = "no date found"
                 abstract_dict[count] = [date,scopus_id]
    df = pd.DataFrame(data=abstract_dict)
    df = (df.T)
    df = df.rename(columns={0: "date", 1: "scopus_id"})
    df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d", errors='coerce')
    df = df[(df['date'] > '2000-01-01')] #Filter in the last 20 years
    return df

####### Pubmed ########

pmids_ud = organize_pmids(pubmed_update_path)
fetched_pmids_ud = retrieved_pmids(pmids_ud)
discard_pmids_ud = discard_pmids_other(pmids_ud)

####### Scopus ########
sids_ud = organize_scopusids(scopus_update_path)
fetched_sids_ud = retrieved_scopus(scopus_update_path)

print(f'Records identified total: {len(fetched_pmids_ud)+len(fetched_sids_ud)+len(discard_pmids_ud)}')
print(f'Records identified from PubMed: {len(fetched_pmids_ud)+len(discard_pmids_ud)}')
print(f'Records identified from Scopus: {len(fetched_sids_ud)}')
print(f'Records discarded from PubMed for other reasons: {len(discard_pmids_ud)}')

##############################################################################
pubmed = pd.read_excel(master, sheet_name="PubMed")
scopus = pd.read_excel(master, sheet_name="Scopus")

print(f'Records excluded due to duplicates total: {(len(fetched_pmids_ud)-len(pubmed))+(len(fetched_sids_ud)-len(scopus))}')
print(f'Records excluded due to duplicates PubMed: {len(fetched_pmids_ud)-len(pubmed)}')
print(f'Records excluded due to duplicates Scopus: {len(fetched_sids_ud)-len(scopus)}')

pubmed_discard = pubmed[pubmed['Agreement']==-1]
pubmed_include = pd.read_excel(master, sheet_name="Full-text include Pubmed")
assert len(pubmed_discard)+len(pubmed_include)==len(pubmed)

scopus_discard = scopus[scopus['Agreement']==-1]
scopus_include = pd.read_excel(master, sheet_name="Full-text include Scopus")
assert len(scopus_discard)+len(scopus_include)==len(scopus)

print(f'Records screened pubmed: {len(pubmed)}')
print(f'Records screened scopus: {len(scopus)}')
print(f'Records screened total: {len(pubmed)+len(scopus)}')
print(f'Records discarded pubmed: {len(pubmed_discard)}')
print(f'Records discarded scopus: {len(scopus_discard)}')
print(f'Records discarded total: {len(pubmed_discard)+len(scopus_discard)}')
print(f'Records sought for retrieval pubmed: {len(pubmed_include)}')
print(f'Records sought for retrieval scopus: {len(scopus_include)}')
print(f'Records sought for retrieval total: {len(pubmed_include)+len(scopus_include)}')

pubmed_real_include = pubmed_include[pubmed_include['include']==1]
pubmed_real_exclude = pubmed_include[pubmed_include['include']==-1]
assert len(pubmed_real_include)+len(pubmed_real_exclude)==len(pubmed_include)

scopus_real_include = scopus_include[scopus_include['include']==1]
scopus_real_exclude = scopus_include[scopus_include['include']==-1]
assert len(scopus_real_include)+len(scopus_real_exclude)==len(scopus_include)

print(f'Finally included pubmed: {len(pubmed_real_include)}')
print(f'Finally included scopus: {len(scopus_real_include)}')
print(f'Finally included total: {len(pubmed_real_include)+len(scopus_real_include)}')

print("Exclusion by reason Pubmed")
pubmed_real_exclude['exlusion reason prisma'].value_counts()
print("Exclusion by reason Scopus")
scopus_real_exclude['exlusion reason prisma'].value_counts()
print("Exclusion by reason Total")
total = pd.concat([pubmed_real_exclude,scopus_real_exclude])
total['exlusion reason prisma'].value_counts()
