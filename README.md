# Neuroscience meets behavior: a systematic literature review on Magnetic Resonance Imaging of the brain combined with real-world digital phenotyping
Authors: Ana María Triana, Jari Saramäki, Enrico Glerean, and Nicholas Mark Edward Alexander Hayward.

This repository contains the code used to produce this systematic literature review findings.

There are two main folders in the repository: src and results. Src contains the code and results contains the figures and supplementary material. Please modify the paths accordingly. 
Due to copyright, we cannot post all the files online. However, we hope this could help the readers in understanding and replicating if needed.

## Retrieving the abstracts automatically
The first step is to retrieve the abstracts from the databases (PubMed and Scopus). To run systematic queries to pubmed, you need to:
1. Edit the keywords files (keys1.csv has the "brain" keywords, keys2.csv has the "behaviour" keywords)
2. Run the script ``. /src/pubmed_search.py``
3. All outputs are stored into a subfolder ``/results/pubmeddata``
4. Retrieve abstracts using ``. /src/pubId2xls.py`` script to make an excel file with the abstracts

To run the systematic queries to Scopus, you need to:
1. Edit the keywords (keys1.csv has the "brain" keywords, keys2.csv has the "behaviour" keywords)
2. Create a file ``apikey.txt`` and store your Elsevier API key there. You need to request this key in advance. Do not reveal the key publicly.
3. Run the script ``/src/scopus_search.py``
4. All output is stored into a subfolder ``/results/scopus``
5. Retrieve abstracts using ``pubId2xls_scopus.py`` script to make an excel file with the abstracts

This process will create an excel file with several tabs, one for each database. If needed, the keywords that yielded the abstract can be added. To do so, run the scripts get_pubmed_keywords.sh and get_scopus_keywords.sh.
Copy paste the keywords in the same order to the excel files. If needed, the get_ok_abstracts.py script retrieves the PubMed IDs for those papers found in Scopus that were not found in PubMed initially.

## Extracting the information
All abstracts are ready to be assessed for inclusion or exclusion. This is manual work. When done, simply move the included abstracts into a new spreadsheet. Then, it is time to start reading. 
We employed Atlas.ti to tag and extract the information with ease. However, we cannot post this file online due to copyright. Briefly, codes can be created in Atlas.ti. These codes will contain information for binary cateogries (e.g. fMRI; 1 includes, 0 excludes). Other qualitative data can be input in the comment section of each paper in Atlas.ti.
Once all the papers are read and all information is tagged, a summary document can be exported to excel. The document will have one column for the paper ID, one column for the comments, and one column for the codes.

## Organizing the information in a table
Once the excel file from Atlas.ti is obtained, create a table organizing the information. Use the script ``full_text.py`` to do so. If needed, the script ``make_final_tables.py`` will help formatting everything as well.

## Graph data analysis
Use the script ``make_network.py`` to create the co-occurence network. This script will create the network as an adjacency matrix and as an edge list. 
Use the edge list to create the link communities according to Yong-Yeol Ahn, James P. Bagrow, and Sune Lehmann, Link communities reveal multiscale complexity in networks, Nature 466, 761 (2010). (https://github.com/bagrow/linkcomm)

# Brain regions analysis
Use the scripts in the folder ./src/brainroisynth to map and extract the number of brain regions used by each MRI-PAD combination. The scripts are labelled in steps, so follow the steps. 

## Visualization
Finally, run the following scripts to create the visualizations.
1. To compute the numbers for the PRISMA dataflow chart, use ``count_prisma.py``
2. To create Tables 3, 4, and 5, and Figures 2, 3, and 5, run the jupyter notebook ``Final_plots.ipynb``
3. Network visualizations are created in gephi. Use the adjacency matrix generated from the Graph analysis and the link communities files as data sources for gephi (https://gephi.org/).
4. To visualize the areas in the brain, use the script ``plot_brainmasks.py``

