clear all
close all
mkdir tempMAout
mkdir MAout
data=readtable('../../results/results_organized.xlsx');
%data=readtable('../../../PMmeta/AAL3/myAAL3/results_organized.xlsx');

dtistudies={
'P18' % there would be also study P13 but we keep it since it has fMRI results
'P19'
'P20'
'P21'
'P22'
'P40'
'P48'
'P49'
'P50'
'P59'
'P71'
'S03'
'S04'
'S18'
'P62'
'P69'
'P73'
'P88'
'P68'
'P70'
'P87'};
filters={'smri', 'physical_activity_sensor'
'smri', 'sleep_sensor'
'smri', 'blood_pressure'
'smri', 'heart_rate'
'smri', 'EMA'
'fmri', 'images'
'fmri', 'app_use'
'fmri', 'temperature'
'fmri', 'glucose'
'fmri', 'physical_activity_sensor'
'fmri', 'sleep_sensor'
'fmri', 'screen_use'
'fmri', 'EMA'
'fmri', 'heart_rate'
'fmri', 'hrv'
'fmri', 'mobility_patterns'
'fmri', 'comm_patterns'
'fmri', ' '
'smri', ' '};

for l=1:size(filters,1)
    l1=filters{l,1};
    l2=filters{l,2};
    if(l2 == ' ' )
       docs=data.ID(find((data.(l1))));
        fileout=[l1 '-ALL-']; 
    else
        docs=data.ID(find((data.(l1).*data.(l2))));
        fileout=[l1 '-' l2 '-' ];
    end
    outcommandB='';
    cc=0;
    for d=1:length(docs)
        if(ismember(docs{d},dtistudies))
            disp(['Skipping dti study ' docs{d}])
            continue
        end
        if(isfile(['output/' docs{d} '.nii.gz']))
            cc=cc+1;
            outcommandB = [outcommandB 'output/' docs{d} '.nii.gz '];
        else
            error(['Missing output/' docs{d} '.nii.gz'])
        end
    end
    fileout=[fileout num2str(cc)];
    outcommandA=['fslmerge -t tempMAout/' fileout '.nii.gz ' ];
    outcommand=[outcommandA outcommandB];
    disp(outcommand)
    system(outcommand)

    disp(['fslmaths tempMAout/' fileout '.nii.gz -Tmean -mul ' num2str(cc) ' MAout/' fileout '.nii.gz'])
    e= system(['fslmaths tempMAout/' fileout '.nii.gz -Tmean -mul ' num2str(cc) ' MAout/' fileout '.nii.gz'])
    e=0;
    if(e==0)
        disp(' ')
    else
        error('stop')
    end
end