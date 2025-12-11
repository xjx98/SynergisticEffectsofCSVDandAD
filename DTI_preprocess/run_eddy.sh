#!/bin/bash
cd $1

# 1: working dir

# Run eddy current correction
if [ ! -d DTI/eddy ]; then
    mkdir DTI/eddy
fi

if [ ! -e DTI/eddy/data.nii.gz ]; then
    eddy --imain=DTI/raw/data.nii --mask=DTI/raw/nodif_brain_mask.nii.gz --acqp=synb0/INPUTS/acqparams.txt --index=DTI/raw/index.txt \
    --bvecs=DTI/raw/data.bvec --bvals=DTI/raw/data.bval --topup=synb0/OUTPUTS/topup --out=DTI/eddy/data --slm=linear --fwhm=2 --ff=5 \
    --very_verbose --repol --rms >> DTI/eddy/eddy_log.txt ###
fi

# Saves the corrected data
if [ ! -e DTI/bvec ]; then
    cp DTI/eddy/data.eddy_rotated_bvecs DTI/bvecs
fi

if [ ! -e DTI/bval ]; then
    cp DTI/raw/data.bval DTI/bvals
fi

if [ ! -e DTI/data_ud.nii.gz ]; then
    cp DTI/eddy/data.nii.gz DTI/data_ud.nii.gz
fi

if [ ! -e DTI/nodif_ud.nii.gz ]; then
    cp synb0/OUTPUTS/b0_u.nii.gz DTI/nodif_ud_syn.nii.gz
fi
