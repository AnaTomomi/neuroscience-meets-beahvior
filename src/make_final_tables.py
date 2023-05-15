import pandas as pd
import numpy as np

path = "./results"

df = pd.read_excel(f'{path}/results_organized.xlsx')

#drop the columns that we will not use for sure
df.drop(columns=['Unnamed: 0', 'ID','abstract'], inplace=True)

#rename the columns
df.rename(columns={"blood_pressure": "blood pressure", "comm_patterns": "communication patterns",
                     "EMA": "EMA", "glucose":"glucose", "heart_rate": "heart rate", 
                     "hrv": "heart rate variability", "mobility_patterns":"mobility patterns", 
                     "images":"images", "physical_activity_sensor":"physical activity", 
                     "resp_rate_sensor":"respiration rate", "screen_use":"screen use", 
                     "sleep_sensor":"sleep", "temperature":"temperature", "app_use":"app use", 
                     "dti":"DWI", "smri":"T1/T2-weighted MRI", "fmri":"fMRI", 
                     "mr_angiography":"angiography", "pubmed_id":"PubMed ID", 
                     "scopus_id": "Scopus ID", "doi":"DOI", "sample_size":"sample size",
                     "patient_sample":"patient sample", "no_diff_diag":"number of diagnostics",
                     "mri_sample":"MRI sample", "rs-fmri":"rs-fMRI", "task-mri":"task fMRI",
                     "no_samples_discarded_mri":"MRI", "no_samples_discarded_dev":"PAD",
                     "no_samples_discarded_other":"other", "dev_discard_reason":"discard reason",
                     "alzheimer":"Alzheimer's", "pediatric_sickle_cell_disease": "pediatric sickel cell disease",
                     "chronic_fatigue_syndrome":"chronic fatigue syndrome", "cognitive_impairment":"cognitive impairement",
                     "knee_osteoarthritis":"knee osteoarthritis", "multiple_sclerosis":"multiple sclerosis",
                     "psichoaffective_dissorder":"psichoaffective disorder", "down_syndrome": "Down syndrome",
                     "chronic_obstructive_pulmonary_disease":"chronic obstructive pulmonary disease",
                     "amb_aut_data":"PAD data analysis", "mri_data_ana":"MRI data analysis", 
                     "brain_area":"brain area", "main_finding":"main finding", "quo":"question of interest",
                     "min_age":"minimum", "max_age":"maximum", "median_age":"median", 
                     "mean_age":"mean", "step_counter":"step counter"}, inplace=True)

#add some columns
df['wristwatch'] = df['wristwatch'] + df['fitness_tracker']

#define the columns to aggregate
study_type_cols = ['cohort', 'controlled-trial', 'cross-sectional', 'intervention', 'longitudinal']
mri_technique_cols = ['DWI', 'T1/T2-weighted MRI', 'fMRI', 'angiography']
fmri_technique_cols = ['rs-fMRI', 'task fMRI']
fmri_stimuli_cols = ['blocked', 'event-related', 'mixed-design', 'naturalistic']
device_function_cols = ['blood pressure', 'communication patterns', 'EMA',
                        'glucose', 'heart rate', 'heart rate variability', 'mobility patterns',
                        'images', 'physical activity', 'respiration rate', 'screen use',
                        'sleep', 'temperature', 'app use']
device_cols = ['actigraph', 'beeper', 'smartphone', 'step counter', 'wristwatch']
illness_cols = ['cognitive_decline', 'depression', 'anorexia', 'psychosis',
                'schizophrenia', 'fibromyalgia', 'obesity', 'stroke', "Alzheimer's",
                'dementia', 'pediatric sickel cell disease', 'chronic fatigue syndrome',
                'cognitive impairement', 'knee osteoarthritis', 'multiple sclerosis',
                'psichoaffective disorder', 'hypertension', 'anxiety', 'Down syndrome',
                'chronic obstructive pulmonary disease', 'bipolar']
age_cols = ['mean', 'median', 'minimum', 'maximum']
discarded_cols =["MRI", "PAD", "other"]

#Define functions that will be needed
def get_age_range(row, cols):
    ages = {col: row[col] for col in cols if row[col] != 0}
    if not ages:
        return np.nan
    return ', '.join(f"{key}: {value}" for key, value in ages.items())

def get_feature(row, cols):
    features = [c for c in cols if row[c]==1]
    return ', '.join(features) if features else np.nan

def convert_to_days(dev_dur_str):
    if not isinstance(dev_dur_str, str) or len(dev_dur_str)<2:
        return None

    num = int(dev_dur_str[:-1])
    unit = dev_dur_str[-1]

    if unit == 'H':
        return num / 24
    elif unit == 'D':
        return num
    elif unit == 'W':
        return num * 7
    elif unit == 'M':
        return num * 30.44  # Approximate days in a month
    else:
        return None

#Start formatting the columns and get rid of the ones I do not need
df['study type'] = df.apply(get_feature, args=(study_type_cols,), axis=1)
df.drop(columns=study_type_cols, inplace=True)

df['MRI technique'] = df.apply(get_feature, args=(mri_technique_cols,), axis=1)
df.drop(columns=mri_technique_cols, inplace=True)

df['fMRI technique'] = df.apply(get_feature, args=(fmri_technique_cols,), axis=1)
df.drop(columns=fmri_technique_cols, inplace=True)

df['fMRI stimulus'] = df.apply(get_feature, args=(fmri_stimuli_cols,), axis=1)
df.drop(columns=fmri_stimuli_cols, inplace=True)

df['device (by function)'] = df.apply(get_feature, args=(device_function_cols,), axis=1)
df.drop(columns=device_function_cols, inplace=True)

df['device'] = df.apply(get_feature, args=(device_cols,), axis=1)
df.drop(columns=device_cols, inplace=True)

df['illness'] = df.apply(get_feature, args=(illness_cols,), axis=1)
df.drop(columns=illness_cols, inplace=True)

df['age range'] = df.apply(get_age_range, args=(age_cols,), axis=1)
df.drop(columns=age_cols, inplace=True)

df['data source'] = 'first-use data'
df.loc[df['independent_study'] == 0, 'data source'] = 're-used'

df['device use'] = df['dev_dur'].apply(convert_to_days)

df['number of discarded samples'] = df.apply(get_age_range, args=(discarded_cols,), axis=1)
df.drop(columns=discarded_cols, inplace=True)

df.drop(columns=['dev_dur','independent_study', 'clinical_setting', 'dur_study', 
                 'mri_freq', 'fitness_tracker', 'multisensor', 'ECG', 'inclinometer',
                 'spo2_sensors', 'dev_freq', 'accelerometer', ], inplace=True)

#Re-order the columns
cols = ['Document', 'PubMed ID', 'Scopus ID', 'title', 'DOI', 'date', 'keywords', 
        'study type', 'data source', 'sample size', 'age range', 'patient sample', 
        'number of diagnostics','illness',
        'MRI technique', 'fMRI technique','fMRI stimulus', 'MRI sample', 
        'device (by function)', 'device', 'device use', 'number of discarded samples',
        'discard reason', 'method', 'MRI data analysis', 'PAD data analysis',
        'question of interest', 'main finding', 'brain area']
df = df[cols]

#save the dataframe
df.to_excel(f'{path}/supplementary_table1.xlsx', index=False)
