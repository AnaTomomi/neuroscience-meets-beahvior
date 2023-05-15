###############################################################################
# This script execute an automated search of all possible combinations of two #
# sets of terms. The search is done in scopus. The script retrieves all scopus#
# results and store them in json files. The search can be updated by changing #
# the "request_string" variable.                                              #
# An API key is needed to run this code. Make sure to request the API key from#
# Scopus and store it in the file ./src/apikey.txt. Do not reveal this key    #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
# Date: 29.10.2020                                                            #
###############################################################################

import requests
import pandas as pd

def save_request(path,term1,term2,counter,r):
    if r.status_code==200: #if the request is ok, save the results
        file = open(path+term1+"_"+term2+"_"+str(counter), "w")
        file.write(r.text)
        file.close()
    else:
        print("Bad output")
        
#Request the results
path = "./results/scopus"
first_terms = pd.read_csv("./src/keys1.csv",header=None)
second_terms = pd.read_csv("./src/keys2.csv",header=None)
first_terms[0] = first_terms[0].str.replace("+", "%20")
second_terms[0] = second_terms[0].str.replace("+", "%20")
first_terms['key'] = 1
second_terms['key'] = 1
combi = pd.merge(first_terms,second_terms,on='key').drop('key',axis=1)

#open the API key
con_file = open("./src/apikey.txt", "r")
apikey = con_file.read()[:-1]
con_file.close()
        
for i,row in combi.iterrows():
    counter = 0
    term1 = row['0_x']
    term2 = row['0_y']
    request_string = "https://api.elsevier.com/content/search/scopus?count=25&start="+str(counter)+"&query=(TITLE("+term1+"%20AND%20"+term2+")%20OR%20ABS("+term1+"%20AND%20"+term2+")%20AND%20LANGUAGE%20(english))%20AND%20DOCTYPE(ar%20OR%20cp)&view=COMPLETE&apiKey="+apikey
    r = requests.get(request_string)
    save_request(path,term1,term2,counter,r)
    print(request_string)
    if int(r.json()['search-results']['opensearch:totalResults'])>25:
        repetitions = int(r.json()['search-results']['opensearch:totalResults'])//25
        if int(r.json()['search-results']['opensearch:totalResults'])%25>0: #if we still have reminders, add one more iteration
            repetitions = repetitions+1
        repetitions = repetitions -1 #because we already have the first request
        while repetitions > 0:
            counter=counter+25
            request_string = "https://api.elsevier.com/content/search/scopus?count=25&start="+str(counter)+"&query=(TITLE("+term1+"%20AND%20"+term2+")%20OR%20ABS("+term1+"%20AND%20"+term2+")%20AND%20LANGUAGE%20(english))%20AND%20DOCTYPE(ar%20OR%20cp)&view=COMPLETE&apiKey="+apikey
            r = requests.get(request_string)
            save_request(path,term1,term2,counter,r)
            print(request_string)
            repetitions = repetitions-1
