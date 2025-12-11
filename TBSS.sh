#!/bin/bash
cd $1

# Prepare files
if [ ! -d DTI/TBSS ]; then
    mkdir DTI/TBSS
fi
cp DTI/dtifit/dti_FA.nii.gz DTI/TBSS
cd DTI/TBSS
# erode a little and zero end slices
tbss_1_preproc *.nii.gz 

# 2. Nonlinear registration, aligning all FA images to a 1x1x1mm standard space.
# Use FMRIB58_FA_1mm as target for nonlinear registration
tbss_2_reg -T

# creating valid mask and mean FA
tbss_3_postreg -T

cd stats
# Thresh set with reference to https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/blob/master/bb_diffusion_pipeline/bb_tbss/bb_tbss_general?ref_type=heads
fslmaths mean_FA_skeleton -thr 2000 -bin mean_FA_skeleton_mask

# Generates skeletonised FA
fslmaths all_FA -mas mean_FA_skeleton_mask all_FA_skeletonised

# Calculate averaged FA within 48 JHU tracts
fslstats -K "$FSLDIR"/data/atlases/JHU/JHU-ICBM-labels-1mm all_FA_skeletonised.nii.gz -M >JHUrois_FA.txt

# Generates average within-tract IDPs for other diffusion measures 
suffix="L1 L2 L3 MO MD"
for elem in $suffix ; do
    applywarp --rel -i $1/DTI/dtifit/dti_$elem -o all_"$elem" -r "$FSLDIR"/data/standard/FMRIB58_FA_1mm -w $1/DTI/TBSS/FA/dti_FA_FA_to_target_warp
    fslmaths all_"$elem" -mas mean_FA_skeleton_mask all_"$elem"_skeletonised
    fslstats -K "$FSLDIR"/data/atlases/JHU/JHU-ICBM-labels-1mm all_"$elem"_skeletonised.nii.gz -M >JHUrois_"$elem".txt
done

# TODO: Generates tract-wise WMHVs
# cd $1
# mask (not threshed) in FLAIR (in T1) -> MNI
#applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_1mm --in=FLAIR/lesions/bianca_mask --warp=transform/T1_to_MNI_warp --out=FLAIR/lesions/lesion_MNI