#!/bin/bash
# Usage:
# 	run_freesurfer.sh <path_to_id_date>
# 
# Run freesurfer in align with UKB pipeline
# https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/blob/master/bb_FS_pipeline/bb_FS_run.sh
# Freesurfer will create a 'fs' dir within $1 for all the output
cd $1
export FREESURFER_HOME=/home/xjx/freesurfer/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
recon-all -all -s fs -i T1/T1_unbiased.nii.gz -sd $1 -FLAIR FLAIR/T2_FLAIR_unbiased.nii.gz -FLAIRpial
