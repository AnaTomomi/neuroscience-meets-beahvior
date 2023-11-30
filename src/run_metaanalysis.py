###############################################################################
# This script runs meta-analysis using NIMARE.                                #
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
###############################################################################

import pandas as pd
import json

import nibabel as nib

from nimare.dataset import Dataset
from nimare.workflows.cbma import CBMAWorkflow

# Define changes
path = '../results'
df = pd.read_excel(f'{path}/metaanalysis.xlsx')
df.rename(columns={"Document":"study_id","sample_sizes":"sample_size"},inplace=True)

# select only those studies with coordinates
df = df[df['x'].notna()]

# check the combinations of MRI - PAD
mri = ['dti', 'smri', 'fmri', 'mr_angiography']
beh = ['blood_pressure', 'comm_patterns', 'EMA', 'glucose', 'heart_rate',
'hrv', 'mobility_patterns', 'images', 'physical_activity_sensor',
'resp_rate_sensor', 'screen_use', 'sleep_sensor', 'temperature',
'app_use']

selected_rows = {}
for mri_col in mri:
    for beh_col in beh:
        key = f"{mri_col}-{beh_col}"
        temp_df = df[(df[mri_col] == 1) & (df[beh_col] == 1)]
        if len(temp_df)>0:
            #only compute for combinations that have more than 5 research papers
            if len(temp_df['study_id'].unique())>=4: 
                selected_rows[key] = temp_df
# add for the MRI alone
for mri_col in mri:
    key = f"{mri_col}-all"
    temp_df = df[(df[mri_col] == 1)]
    if len(temp_df)>0:
        #only compute for modalities that have more than 5 research papers
        if len(temp_df['study_id'].unique())>=4: 
            selected_rows[key] = temp_df

for key in selected_rows.keys():
    print(key)
    #create the file names
    save_img = f'{path}/metaanalysis/img_{key}_all.nii'
    cluster_path = f'{path}/metaanalysis/cluster_{key}_all.xlsx'
    contribution_path = f'{path}/metaanalysis/contribution_{key}_all.xlsx'
    
    # create the dataset into json files
    dset_dict = {}
    coords_df = selected_rows[key][['study_id', 'x', 'y', 'z', 'space', 'sample_size']]
    
    study_ids = list(coords_df['study_id'].unique())
    
    for study in study_ids:
        this_study_coords = coords_df[coords_df['study_id'] == study]
        contrast = {"coords": { "space": this_study_coords['space'].unique()[0],
                    "x": list(this_study_coords['x']),
                    "y": list(this_study_coords['y']),
                    "z": list(this_study_coords['z'])},
                    "images":{},
                    "metadata": {"sample_sizes":[float(this_study_coords['sample_size'].unique()[0])]}
                    }
        dset_dict[study] = {"contrasts": {"1": contrast }}
    with open(f'{path}/metaanalysis/{key}.json', 'w') as fp:
        json.dump(dset_dict, fp)

    # now run the meta analysis
    dset_file = f'{path}/metaanalysis/{key}.json'
    dset = Dataset(dset_file)

    workflow = CBMAWorkflow(corrector="fdr")
    result = workflow.fit(dset)

    # summarize the results
    img = result.get_map("z_corr-FDR_method-indep")
    nib.save(img, save_img)
    
    if len(result.tables["z_corr-FDR_method-indep_tab-clust"])>0:
        cluster = result.tables["z_corr-FDR_method-indep_tab-clust"]
        contribution = result.tables["z_corr-FDR_method-indep_diag-Jackknife_tab-counts_tail-positive"]
        cluster.to_excel(cluster_path)
        contribution.to_excel(contribution_path)