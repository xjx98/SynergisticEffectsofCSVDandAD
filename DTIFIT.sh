#!/bin/bash
cd $1

# 1. Prepare file to run dti fit
# ditfit needs: data_ud nodif_brain_mask_ud 

# Extract b0
if [ ! -e DTI/nodif_ud.nii.gz ]; then
    fslroi DTI/data_ud.nii.gz DTI/nodif_ud.nii.gz 0 1
fi

# Skull stripping on nodif image and generates nodif brain mask for eddy (distorted)
if [ ! -e DTI/nodif_brain.nii.gz ]; then
    bet DTI/nodif_ud.nii.gz DTI/nodif_ud_brain -m -f 0.2 -o -R
fi

# 2. Fit the diffusion tensor model using dtifit
# Automatically creates FA, MD, MO maps from a diffusion MRI
if [ ! -d DTI/dtifit ]; then
    mkdir DTI/dtifit
fi

if [ ! -e DTI/dtifit/dti_FA.nii.gz ]; then
    dtifit -k DTI/data_ud -m DTI/nodif_ud_brain_mask -r DTI/bvecs -b DTI/bvals -o DTI/dtifit/dti
fi
