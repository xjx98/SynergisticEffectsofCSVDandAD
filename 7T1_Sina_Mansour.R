# ******************************************************************************************************************
# Generates T1 naive parcellations for each participant following UKB pipeline
# This is an adaption from https://github.com/sina-mansour/UKB-connectomics/blob/main/scripts/bash/UKB_connectivity_mapping_pipeline.sh
# Adapted for Gulou images
# Please cite 10.1016/j.neuroimage.2023.120407
# Created by Haojie Chen, Beijing Normal University, 20240401
# Requires: fsl, freesurfer, workbench, python3 (see above link for related modules)
# ******************************************************************************************************************

library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(glue)
#设置路径
path.nifti <- '/home/xjx/ADX/IMAGE' # Absolute dir for images for each participant
path.QC <- '/home/xjx/ADX/QC' # Absolute dir for quality control


# Constant variables and functions -----------------------

# Shell and python scripts
script_freesurfer <- '/home/xjx/ADX/run_freesurfer.sh'


# path_config <- '/home/shulab/GULOU/config/' # Dir for scripts and standard templates
#script_initialize <- '/home/shulab/GULOU/config/scripts/config_init.sh'
#script_step1 <- '' # Step 1: map surface atlases to native surface
#script_step2 <- '' # Step 2: map native surface atlases to native volumetric labels
#script_step3 <- '' # Step 3: map subcortical labels from standard to native space

# Functions
create_when_absent <<- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}
remove_images <- function(path.image, pattern, process.num = 10) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.image), function(ID){
    path_image = glue('{path.image}/{ID}/{pattern}')
    system(glue('rm {path_image}'))
  })
}

# Run freesurfer for all participants in the dir
# Reference runtime: 3H ~ 6H single participant (7.4.1)
fs_dir <- function(path.input, process.num = 8) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    system(str_c(script_freesurfer, ' ', path_id))
  })
}

# Using slicesdir for quality control
slices_dir_single <- function(path.input, image, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list <- str_c(path.input, '/', list.files(path.input), '/', image)
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir {image.list}'))
}


slices_dir_pair <- function(path.input, imageBase, imageLine, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list.A <- str_c(path.input, '/', list.files(path.input), '/', imageBase)
  image.list.B <- str_c(path.input, '/', list.files(path.input), '/', imageLine)
  image.list <- purrr::map2_chr(image.list.A, image.list.B, function(x, y) stringr::str_c(x, ' ', y))
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir -o {image.list}'))
}

# Run freesurfer and QC ----------------------------------------------
# Run freesurfer, which takes time
fs_dir(path.input = path.nifti, process.num = 10)

# Quality control
# Generates surf files of nifti format
future::plan(future::multisession, workers = 10)
furrr::future_walk(list.files(path.nifti), function(ID){
    path_id <- str_c(path.nifti, '/', ID)
    system(glue('cd {path_id}'))
    setwd(path_id)
    Sys.setenv(FREESURFER_HOME='/home/xjx/freesurfer/freesurfer')
    system('/home/xjx/freesurfer/freesurfer/bin/mri_convert fs/mri/T1.mgz T1_fs.nii.gz')
    system('/home/xjx/freesurfer/freesurfer/bin/mri_convert fs/mri/aparc+aseg.mgz aparc+aseg.nii.gz')
})

# QC
slices_dir_pair(path.input = path.nifti, imageBase = 'T1_fs.nii.gz', imageLine = 'aparc+aseg.nii.gz', path.report = str_c(path.QC, '/freesurfer'))

# Remove temp images after QC
remove_images(path.image = path.nifti, pattern = 'T1_fs.nii.gz', process.num = 10)
remove_images(path.image = path.nifti, pattern = 'aparc+aseg.nii.gz', process.num = 10)

# TODO: After running freesurfer, generate individualized T1 parcellations ------------------------------------------ 
# Initialize configuration files for individualized T1 parcellation, including copying related templates
# Step 1: map surface atlases to native surface


# Step 2: map native surface atlases to native volumetric labels

# Step 3: map subcortical labels from standard to native space