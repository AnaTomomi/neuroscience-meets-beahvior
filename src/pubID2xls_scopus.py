###############################################################################
# This script reads scopus ids and retrieves the title, abstract, and publica #
# tion date for each of them. The script also compares the results with pubmed#
# results (that were previously organized using pubID2xls.py) and removes du- #
# plicates and saves the results in an excel file.                            #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
# Date: 29.10.2020                                                            #
###############################################################################

from elsapy.elsclient import ElsClient
import json
import os
import pandas as pd
import numpy as np


#Load main file of results
rpath="./results/literature.xlsx"
savepath = './results/literature.xlsx'
savepath2 = './results/literature2.xlsx'
abstracts = pd.read_excel(rpath, sheet_name="Abstracts",header=0)
no_abstracts = pd.read_excel(rpath, sheet_name="No_Abstracts",header=0)

abstracts["title"] = abstracts["title"].str.strip(".")
no_abstracts["title"] = no_abstracts["title"].str.strip(".")

## Initialize client
#client = ElsClient(config['apikey'])

#Read files 
path = './results/scopus/'
files = [os.path.join(path,x) for x in os.listdir(path)]

#Quick summary of results
total = 0
for i,file in enumerate(files):
    with open(file) as json_file:
        data = json.load(json_file)
        total = total + int(data['search-results']['opensearch:totalResults'])
print("total files %s" % str(total))

def check_pubmedid(abstracts, no_abstracts, item):
    # Checks if the result already exists in the pubmed id database search by pubmedid. 
    
    pubmed_id = np.int64(item['pubmed-id'])
    if not(pubmed_id in abstracts["pubmed_id"].tolist() or pubmed_id in no_abstracts["pubmed_id"].tolist()): #Is this a new paper?
        scopus_id = item['dc:identifier'].strip('SCOPUS_ID:')
        title = item['dc:title']
        if 'prism:doi' in item.keys():
            doi= item['prism:doi']
        else:
            doi="no doi found"
        if 'dc:description' in item.keys():
            abstract = item['dc:description']
        else:
            abstract = "no abstract found"
        if 'prism:coverDate' in item.keys():
            date = item['prism:coverDate']
        else:
            date = "no date found"
    else:
        pubmed_id = "nothing"
        scopus_id = "nothing"
        title = "nothing"
        doi = "nothing"
        abstract = "nothing"
        date = "nothing"
    return pubmed_id, scopus_id, title, doi, abstract, date

def check_title(abstracts, no_abstracts, item):
    # Checks if the result already exists in the pubmed id database search by title comparison. 
    
    if 'title' in item.keys():
        title = item["title"]
    elif "dc:title" in item.keys():
        title = item["dc:title"]
        
    if not(title in abstracts["title"].tolist() or title in no_abstracts["title"].tolist()):
        scopus_id = item['dc:identifier'].strip('SCOPUS_ID:')
        pubmed_id = "no pubmed id"
        if 'prism:doi' in item.keys():
            doi= item['prism:doi']
        else:
            doi="no doi found"
        if 'dc:description' in item.keys():
            abstract = item['dc:description']
        else:
            abstract = "no abstract found"
        if 'prism:coverDate' in item.keys():
            date = item['prism:coverDate']
        else:
            date = "no date found"
    else:
        pubmed_id = "nothing"
        scopus_id = "nothing"
        title = "nothing"
        doi = "nothing"
        abstract = "nothing"
        date = "nothing"
    return pubmed_id, scopus_id, title, doi, abstract, date

#Organize results
abstract_dict = {}
for i,file in enumerate(files):
    with open(file) as json_file:
        data = json.load(json_file)
    if int(data['search-results']['opensearch:totalResults'])>0: #does it have results?
        for j,item in enumerate(data['search-results']['entry']):
            if 'pubmed-id' in item.keys(): #does it have pubmed id?
                pubmed_id, scopus_id, title, doi, abstract, date  = check_pubmedid(abstracts, no_abstracts, item)
                abstract_dict[scopus_id] = [pubmed_id, title, doi, abstract, date]
            elif ('title' in item.keys() or "dc:title" in item.keys()): #if not pubmedid, does it have a title?
                pubmed_id, scopus_id, title, doi, abstract, date  = check_title(abstracts, no_abstracts, item)
                abstract_dict[scopus_id] = [pubmed_id, title, doi, abstract, date]
            else:
                print("file %s and item %s not found \n" % str(i) % str(j) )
    else:
        del data

df = pd.DataFrame(data=abstract_dict)
df = (df.T)
df = df.drop_duplicates()
df = df.rename(columns={0: "pubmed_id", 1: "title", 2:"doi", 3:"abstract", 4:"date"})
df.index.names = ['scopus_id']
df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d", errors='coerce')
df = df[(df['date'] > '2000-01-01')] #Filter in the last 20 years
df['date'] = df['date'].dt.strftime('%Y-%m-%d')

#Write in excel
with pd.ExcelWriter(savepath2) as writer:  
    df.to_excel(writer, sheet_name='Scopus')
