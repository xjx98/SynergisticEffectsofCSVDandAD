# Created by Haojie Chen, Beijing Normal University, 20240104
library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(glue)
#设置路径
path.nifti <- '/home/xjx/sharefolder/IMAGE' # Absolute path of nifti
path.QC <- '/home/xjx/sharefolder/QC' # Absolute path of quality control
path.IDP <- '/home/xjx/sharefolder/IDP'
path.shell.bianca <- '/home/xjx/sharefolder/run_bianca_single.sh'
path.bianca.train.data <- '/home/xjx/soft/bianca/bianca_class_data'

# Related functions -------------
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}
get_file_extension <- function(file_path) {
  # Split the file name by the dot and extract the last part
  parts <- strsplit(file_path, "\\.")[[1]]
  if (length(parts) > 1) {
    return(tail(parts, n=1))
  } else {
    return('')  # Return empty string if there is no file extension
  }
}

# Drop the file extension (e.g. X.nii.gz -> X)
drop_extension <- function(path) {
  if(!grepl(pattern = '\\..[^\\.]*$', path)) return(path)
  drop_extension(sub('\\..[^\\.]*$', '', path)) # Recursively removes all file extension (e.g. .nii.gz)
}
remove_images <- function(path.image, pattern, process.num = 10) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.image), function(ID){
    path_image = glue('{path.image}/{ID}/{pattern}')
    system(glue('rm {path_image}'))
  })
}

slices_dir_single <- function(path.input, image, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list <- str_c(path.input, '/', list.files(path.input), '/', image)
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir {image.list}'))
}

slices_dir_selected_id <- function(list_id_path, image, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list <- str_c(list_id_path, '/', image)
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir {image.list}'))
}

# QC with slices dir (paired)
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

# Functions for BIANCA
bianca_dir <- function(path.input, process.num = 8) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    system(str_c(path.shell.bianca, ' ', path_id, ' ', path.bianca.train.data))
  })
}

copy_IDP <- function(group_label, IDP_label, path.input, path.IDP.input, path.IDP.output, process.num = 10) {
  create_when_absent(glue('{path.IDP.output}/{IDP_label}'))
  create_when_absent(glue('{path.IDP.output}/{IDP_label}/{group_label}'))
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    IDP.format = get_file_extension(path.IDP.input)
    file.copy(
      from = glue('{path.input}/{ID}/{path.IDP.input}'),
      to = glue('{path.IDP.output}/{IDP_label}/{group_label}/{ID}.{IDP.format}')
    )
  })
}


# Run BIANCA -------
bianca_dir(path.input = path.nifti, process.num = 3)

slices_dir_pair(path.input = path.nifti, imageBase = 'FLAIR/T2_FLAIR_unbiased.nii.gz', imageLine = 'FLAIR/lesions/final_mask.nii.gz', path.report = str_c(path.QC, '/BIANCA_FINAL'))# QC

copy_IDP(group_label = 'AD_MCIX', path.IDP.input = 'FLAIR/lesions/volume.txt', IDP_label = 'WMHV_BIANCA', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical segmentation
