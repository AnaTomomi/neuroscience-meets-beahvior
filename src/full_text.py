import pandas as pd

#define paths
data_location = "./results/literature-review-results.xlsx"
master_location = "./results/Literature Review.xlsx"
save_location = "./results/results_organized.xlsx"

data = pd.read_excel(data_location) #Read the data

# Define those items that were markers for humans, but useless to computers
items_removed = ['sample_size', 'min_age', 'max_age', 'median_age', 'mean_age', 'no_diff_diag', 
                 'patient_sample', 'no_samples_discarded_mri', 'no_samples_discarded_dev', 
                 'dev_discard_reason', 'no_samples_discarded_other', 'dur_study', 'dev_freq', 
                 'dev_dur', 'mri_sample', 'mri_freq', 'method', 'amb_aut_data', 'mri_data_ana', 
                 'quo', 'brain_area', 'main_finding', 'conclusion']

dfout = pd.DataFrame()

# Let's put the tags and comments into a column format
for index, row in data.iterrows():
    # for each row, separate the data with the strange separator from AtlasTI
    separator="\u2029"
    article_comments = row['Comment'].split(separator)
    
    if not pd.isna(row['Codes']):
        article_codes_temp = row['Codes'].split('\n')
    article_codes = list(set(article_codes_temp) - set(items_removed))
    
    # prepare the dictionary to store the data in the current row
    tempdict={}
    tempdict['ID'] = row['ID']
    tempdict['Document'] = row['Document']
    
    for item in article_codes:
        tempdict[item] = 1 #if we are here, the tag is 
    
    for item in article_comments:
        # for each item that was tag, separate the key and the value
        temp=item.split(':')
        if(len(temp) <= 0):
            print("This should not happen")
            print(item)
        
        if(len(temp)==1):
            tempdict[temp[0]]="" #if we are here, the key has no value
        else:
            tempdict[temp[0]]=": ".join(temp[1:]) # note the 1: to include also those cases where a ":" is in the value part
         
    tempdf = pd.DataFrame(tempdict, index=[0]) # make a dataframe out of the dictionary
    
    #dfout = dfout.append(tempdf, ignore_index = True)
    dfout = pd.concat([dfout, tempdf]) # append the dataframe to the output df

#Drop columns with no data (markers for humans)    
dfout.drop(columns=['sample description', 'data quality', 'duration of the study', 
                    'duration of the study', 'duration of the study', 'frequency sensor data', 
                    'sample frequency mri', 'duration sensor use', 'analysis methods', 
                    'question_of_interest', 'findings', 'conclusion'], inplace=True)

#Format the data into correct datatypes
cols = ['sample_size', 'min_age', 'max_age', 'median_age', 'mean_age', 
        'no_diff_diag', 'patient_sample', 'no_samples_discarded_mri', 
        'no_samples_discarded_dev', 'no_samples_discarded_other', 'mri_sample']
dfout[cols] = dfout[cols].apply(pd.to_numeric, errors='coerce')
dfout.fillna(0, inplace=True)

#Merge the dti and diff-mri information because it describes the same thing
dfout['dti'] = dfout['dti'] + dfout['diff-mri']
dfout.drop(columns=['diff-mri', 'fractional_anisotropy'], inplace=True)

#Now let's merge this with information about the date, title, pubmed or scopus ID,
#and keywords, as taken from the master document

pubmed = pd.read_excel(master_location, sheet_name="Full-text include Pubmed")
scopus = pd.read_excel(master_location, sheet_name="Full-text include Scopus")

master = pd.concat([pubmed, scopus])
master = master[master["include"]==1]
master = master[['Atlas.ti code', 'pubmed_id', 'doi', 'title', 'date', 'keywords', 'scopus_id', 'abstract']]
master.rename(columns={"Atlas.ti code": "Document"}, inplace=True)

dfout = dfout.merge(master, how="outer", on="Document")

#Re-order the column names
dfout = dfout[['ID', 'Document', 'pubmed_id', 'scopus_id', 'title', 'doi', 'date', 'keywords', 
               'cohort', 'controlled-trial', 'cross-sectional', 'intervention', 'longitudinal',
               'independent_study', 'sample_size', 'min_age', 'max_age', 'median_age', 'mean_age',
               'clinical_setting', 'patient_sample', 'no_diff_diag', 'dur_study',
               'dti', 'smri', 'fmri', 'mr_angiography', 'mri_sample', 'mri_freq', 
               'rs-fmri', 'task-mri', 'blocked', 'event-related', 'mixed-design', 'naturalistic',
               'blood_pressure', 'comm_patterns', 'EMA', 'glucose', 'heart_rate', 'hrv', 'mobility_patterns', 'images',
               'physical_activity_sensor', 'resp_rate_sensor', 'screen_use', 'sleep_sensor', 'temperature', 'app_use',
               'actigraph', 'beeper', 'fitness_tracker', 'smartphone', 'step_counter', 'wristwatch', 'multisensor',
               'accelerometer', 'ECG', 'inclinometer', 'spo2_sensors', 
               'dev_dur', 'dev_freq', 
               'no_samples_discarded_mri', 'no_samples_discarded_dev', 'no_samples_discarded_other', 'dev_discard_reason', 
               'cognitive_decline', 'depression','anorexia', 'psychosis', 'schizophrenia', 'fibromyalgia', 'obesity', 
               'stroke', 'alzheimer', 'dementia', 'pediatric_sickle_cell_disease', 'chronic_fatigue_syndrome',
               'cognitive_impairment', 'knee_osteoarthritis', 'multiple_sclerosis', 'psichoaffective_dissorder', 
               'hypertension', 'anxiety', 'down_syndrome', 'chronic_obstructive_pulmonary_disease', 'bipolar',
               'method', 'amb_aut_data', 'mri_data_ana',
               'brain_area', 'main_finding', 'quo', 'abstract', '']]

#Save to excel
dfout.to_excel(save_location)