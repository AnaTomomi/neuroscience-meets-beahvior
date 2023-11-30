clear all
close all
%%
% check the miss
% process the strings to get AAL in front
% run fslmerge from here

% list of dti stuies to exclude
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

map=readtable('map.xlsx');
% manual fixes
% hyppocampus -> hippocampus
% cerebral vermis -> cerebellar vermis

%data2=readtable('datav2.xlsx');
opts = detectImportOptions('../../results/results_organized.xlsx');
opts.SelectedVariableNames = opts.SelectedVariableNames([2, 93]); % columns 'ID' and 'brain_area'

data = readtable('../../results/results_organized.xlsx',opts);
% for d=1:height(data)
%     out=strcmp(strtrim(data{d,2}{1}),strtrim(data2{d,2}{1}));
%     if(out == 0)
%         d
%         disp([data{d,2} ' ' data2{d,2}])
%     end
% end
% error('stop')
mkdir tempoutput
mkdir output
for m = 1:height(data)
    brain_areas = split(data{m,2},',');
    studylabel=data{m,1}{1};
    if(ismember(studylabel,dtistudies))
        disp(['Skipping dti study ' studylabel])
        continue
    end
    out='';
    %disp(studylabel)
    if length(brain_areas{1})==0
        %disp(['Study ' studylabel ' has no brain areas' ])
        disp(studylabel)
        out='NaN';
        
    else
        for ba = 1:length(brain_areas)
            curr=strtrim(brain_areas{ba});
            curr=strrep(curr,'.','');
            curr=strrep(curr,'  ',' ');
            if(length(curr)==0) continue; end
            if(length(curr)>3 && strcmp(curr(1:4),'and '))
                %disp(curr)
                curr=curr(5:end);
                %disp(curr)
            end
            id=find(strcmp(curr,map.ORIGINAL));
            if length(id)==0
                error(['miss: ' curr])
            else
                out=[out ' ' map.AAL{id}];
                %disp([studylabel ':' out])
            end
        end
    end
    out=strrep(out,'  ',' ');
    outarray=split(strtrim(out),' ');
    newout='';
    for c = 1:length(outarray)
        newout=[newout 'AAL/' outarray{c} '.nii.gz '];
    end
    disp(['fslmerge -t' ' tempoutput/' studylabel '.nii.gz ' newout])

    e=system(['fslmerge -t' ' tempoutput/' studylabel '.nii.gz ' newout])
    if(e ~= 0)
        out
        error('missing image ')
    end
    disp(['fslmaths tempoutput/' studylabel '.nii.gz -Tmean -bin output/' studylabel '.nii.gz' ])
    e=system(['fslmaths tempoutput/' studylabel '.nii.gz -Tmean -bin output/' studylabel '.nii.gz' ])
    if(e ~= 0)
        error('missing image ')
    end
end
