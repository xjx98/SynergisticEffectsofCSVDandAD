#!/bin/bash
cd $1

# 1: working dir

# Run synb0
#if [ ! -e synb0/OUTPUTS/topup_fieldcoef.nii.gz ]; then
sudo -S docker run --rm -v $1/synb0/INPUTS/:/INPUTS/ -v $1/synb0/OUTPUTS:/OUTPUTS/ -v $1/synb0/license.txt:/extra/freesurfer/license.txt --user $(id -u):$(id -g) leonyichencai/synb0-disco:v3.0 --stripped << EOF
1
EOF
