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
shell.dtifit <- '/home/xjx/soft/dtifit_tbss/DTIFIT.sh'
shell.tbss <- '/home/xjx/soft/dtifit_tbss/TBSS.sh'
# ToQC: 
  # DTI/nodif_ud_brain_overlay.nii.gz

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

# Functions for BIANCA#DTI预处理：提取b0，去头骨，生成FA、MD、MO等参数的图像
dtifit_dir <- function(path.input, process.num = 8) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    system(str_c(shell.dtifit, ' ', path_id))
  })
}

#TBSS分析：用分析软件tbss_1_preproc预处理，tbss_2_reg非线性配准，tbss_3_postreg后配准处理，stats下统计分析生成骨架FA，计算JHU48条束内的平均扩散参数
tbss_dir <- function(path.input, process.num = 8) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    system(str_c(shell.tbss, ' ', path_id))
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



# Check if all data ready--------------------------------------------------------------

path.id <- str_c(path.nifti, '/', list.files(path.nifti))

table_check_dti <- tibble(
  path.id = path.id
) %>% mutate(
  dti = map_int(path.id, ~if_else(file.exists( str_c(.,'/DTI/data_ud.nii.gz') ), 1, 0)),
  bvecs = map_int(path.id, ~if_else(file.exists( str_c(.,'/DTI/bvecs') ), 1, 0)),
  bvals = map_int(path.id, ~if_else(file.exists( str_c(.,'/DTI/bvals') ), 1, 0))  
)
# Run DTIFIT ------------------------
dtifit_dir(path.nifti, process.num = 10)

# Quality control of brain stripping on b0 images
slices_dir_single(path.input = path.nifti, image = 'DTI/nodif_ud_brain_overlay.nii.gz', path.report =  str_c(path.QC, '/nodif_ud_bet'))

# Quality control of FA images
slices_dir_single(path.input = path.nifti, image = 'DTI/dtifit/dti_FA.nii.gz', path.report =  str_c(path.QC, '/dtifit_fa'))

# Run TBSS ---------------------------
# Copy all FA images into TBSS
tbss_dir(path.nifti, process.num = 5)

# Quality control of FA images
slices_dir_single(path.input = path.nifti, image = 'DTI/TBSS/stats/all_FA.nii.gz', path.report =  str_c(path.QC, '/tbss_fa'))
# Quality control of MD images
slices_dir_single(path.input = path.nifti, image = 'DTI/TBSS/stats/all_MD.nii.gz', path.report =  str_c(path.QC, '/tbss_md'))

# Save mean FA within 48 tracts
copy_IDP(group_label = 'test', path.IDP.input = 'DTI/TBSS/stats/JHUrois_FA.txt', IDP_label = 'TBSS_FA', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical segmentation

# Save mean MD within 48 tracts
copy_IDP(group_label = 'test', path.IDP.input = 'DTI/TBSS/stats/JHUrois_MD.txt', IDP_label = 'TBSS_MD', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical segmentation


# Save mean MO within 48 tracts
copy_IDP(group_label = 'test', path.IDP.input = 'DTI/TBSS/stats/JHUrois_MO.txt', IDP_label = 'TBSS_MO', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical segmentation
