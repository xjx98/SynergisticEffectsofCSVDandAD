#!/bin/bash
cd $1

# Aligh FLAIR with T1 and perform intensity correction on FLAIR
if [ ! -e FLAIR/T2_FLAIR_unbiased.nii.gz ]; then
    flirt -in FLAIR/data.nii -ref T1/T1_unbiased.nii.gz -out FLAIR/FLAIR_in_T1.nii.gz -omat transform/FLAIR_to_T1.mat -interp nearestneighbour
    fslmaths FLAIR/FLAIR_in_T1.nii.gz -div T1/T1_fast/T1_brain_bias.nii.gz FLAIR/T2_FLAIR_unbiased.nii.gz
fi


if [ ! -d FLAIR/lesions ]; then
    mkdir FLAIR/lesions
fi

# Create an inclusion mask with T1 --> Used to remove GM from BIANCA results
if [ ! -e FLAIR/lesions/T1_unbiased_brain_mask.nii.gz ]; then
    make_bianca_mask T1/T1_unbiased.nii.gz T1/T1_fast/T1_brain_pve_0.nii.gz transform/T1_to_MNI_warp_coef_inv.nii.gz
    #TO CHECK: Move the inclusion mask to T2_FLAIR/lesions directory
    mv T1/T1_unbiased_bianca_mask.nii.gz T1/T1_unbiased_ventmask.nii.gz T1/T1_unbiased_brain_mask.nii.gz FLAIR/lesions/
fi

# Generate the configuration file to run Bianca
if [ ! -e FLAIR/lesions/conf_file.txt ]; then
    echo T1/T1_unbiased_brain.nii.gz FLAIR/T2_FLAIR_unbiased.nii.gz transform/T1_to_MNI_linear.mat > FLAIR/lesions/conf_file.txt;
fi

# Run BIANCA
if [ ! -e FLAIR/lesions/bianca_mask.nii.gz ]; then
    bianca --singlefile=FLAIR/lesions/conf_file.txt --querysubjectnum=1 --brainmaskfeaturenum=1 --loadclassifierdata=$2 --matfeaturenum=3 --featuresubset=1,2 -o FLAIR/lesions/bianca_mask
fi

# Apply the inclusion mask to BIANCA output to get the final thresholded mask
# Also computes the lesion volume
if [ ! -e FLAIR/lesions/volume.txt ]; then
    fslmaths FLAIR/lesions/bianca_mask -mul FLAIR/lesions/T1_unbiased_bianca_mask.nii.gz -thr 0.8 -bin FLAIR/lesions/final_mask
    #Get the volume of the lesions
    fslstats FLAIR/lesions/final_mask -V | awk '{print $1}' > FLAIR/lesions/volume.txt
fi
