#!/bin/bash
cd $1

# 1: working dir
# 2: dir of freesurfer license, e.g. /media/lairai/Lairai10T/gulou/license.txt
# 3: dir of acqparams.txt e.g.  /media/lairai/Lairai10T/gulou/acqparams.txt 
# 4: dir of generate_index.sh, e.g. /media/lairai/Lairai10T/gulou/codes/shell/generate_index.sh
# 5: dir of volume number, e.g. 33


# Move original data into raw directory

if [ ! -d DTI/raw ]; then
    mkdir DTI/raw
    cd DTI
    find . -maxdepth 1 -not -name 'raw' -not -name '.' -exec mv {} raw/ \;
    cd ..
fi

# Extract b0
if [ ! -e DTI/raw/nodif.nii.gz ]; then
    fslroi DTI/raw/data.nii DTI/raw/nodif.nii.gz 0 1 ###
fi

# Skull stripping on nodif image and generates nodif brain mask for eddy (distorted)
if [ ! -e DTI/raw/nodif_brain.nii.gz ]; then
    bet2 DTI/raw/nodif DTI/raw/nodif_brain -m -f 0.2 ###
fi

# Prepare for synb0 procedure
if [ ! -d synb0 ]; then
    mkdir synb0
fi

if [ ! -d synb0/INPUTS ]; then
    mkdir synb0/INPUTS
fi

if [ ! -d synb0/OUTPUTS ]; then
    mkdir synb0/OUTPUTS
fi

if [ ! -e synb0/INPUTS/T1.nii.gz ]; then
    cp T1/T1_brain.nii.gz synb0/INPUTS/T1.nii.gz
fi

if [ ! -e synb0/INPUTS/b0.nii.gz ]; then
    cp DTI/raw/nodif.nii.gz synb0/INPUTS/b0.nii.gz
fi

if [ ! -e synb0/license.txt ]; then
    cp $2 synb0/license.txt
fi

if [ ! -e synb0/INPUTS/acqparams.txt ]; then
    cp $3 synb0/INPUTS/acqparams.txt
fi

# Generate index file for eddy correction
if [ ! -e DTI/raw/index.txt ]; then
    $4 $5 DTI/raw
fi
