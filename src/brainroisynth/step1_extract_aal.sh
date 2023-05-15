#!/bin/bash

mkdir -p AAL
for n in $(cat AAL3v1.nii.txt |cut -d\  -f1); do 
	outname=$(grep "^$n " AAL3v1.nii.txt |cut -d\  -f2);
	echo fslmaths AAL3v1.nii.gz -thr $n -uthr $n  AAL/$outname".nii.gz"
	fslmaths AAL3v1.nii.gz -thr $n -uthr $n  AAL/$outname".nii.gz"
done

# let's make also an empty nii file
fslmaths AAL3v1.nii.gz -mul 0 AAL/NaN.nii.gz
