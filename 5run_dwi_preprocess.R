# Created by Haojie Chen, Beijing Normal University, 20240219
library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(glue)
#设置路径
path.nifti <- '/home/xjx/sharefolder/IMAGE' # Absolute path of nifti
path.QC <- '/home/xjx/sharefolder/QC' # Absolute path of quality control
path.IDP <- '/home/xjx/sharefolder/IDP'
path.freesurfer.lisence <- '/home/xjx/soft/DTI_preprocess/license.txt'
path.shell.file.prepare <- '/home/xjx/soft/DTI_preprocess/dwi_preprocess/prepare_file.sh'
path.shell.synb0 <- '/home/xjx/soft/DTI_preprocess/dwi_preprocess/run_synb0.sh'
path.shell.eddy <- '/home/xjx/soft/DTI_preprocess/dwi_preprocess/run_eddy.sh'
path.shell.generate.index <- '/home/xjx/soft/DTI_preprocess/generate_index.sh'

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

# Generate acq para file (single volume)参数文件，相位编码
acq_para_gen <- function(phase_encoding, read_out, acq_path) {
  pe = case_when(
    phase_encoding == 'i' ~ '1 0 0',
    phase_encoding == 'i-' ~ '-1 0 0',
    phase_encoding == 'j' ~ '0 1 0',
    phase_encoding == 'j-' ~ '0 -1 0',
    phase_encoding == 'k' ~ '0 0 1',
    TRUE ~ '0 0 -1' # k-
  )
  str_c(
    str_c(pe, ' ', read_out),
    str_c(pe, ' 0'),
    sep = '\n'
  )  %>% readr::write_file(acq_path)
}

# Prepare files for synb0 and eddy
file_prepare_dir <- function(path.input, process.num = 10) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    vol.num = readr::read_lines(str_c(path_id, '/DTI/data.bval')) %>% trimws() %>% strsplit(' ') %>% unlist() %>% length()
    para.parsed <- jsonlite::fromJSON(str_c(path_id, '/DTI/data.json'))
    pe = ifelse(!is_null(para.parsed$PhaseEncodingDirection), para.parsed$PhaseEncodingDirection, 'j-') # Direction of phase encoding
    rt = ifelse(!is_null(para.parsed$TotalReadoutTime), para.parsed$TotalReadoutTime, para.parsed$EstimatedTotalReadoutTime) # Readout time
    acq.path = str_c(path_id, '/acqparams_dti.txt')
    acq_para_gen(pe, rt, acq.path)
    system(str_c(path.shell.file.prepare, ' ', path_id, ' ', path.freesurfer.lisence, ' ', acq.path, ' ', path.shell.generate.index, ' ', vol.num))
  })
}



# Run synb0
synb0_dir <- function(path.input, process.num = 3) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    if(!file.exists(str_c(path_id, 'synb0/OUTPUTS/topup_fieldcoef.nii.gz'))){
      system(str_c(path.shell.synb0, ' ', path_id))
    }
  })
}

# Run eddy
eddy_dir <- function(path.input, process.num = 5) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    system(str_c(path.shell.eddy, ' ', path_id))
  })
}

# Run preprocess -------


# Prepare file for synb0 and eddy
file_prepare_dir(path.nifti, 5)

# Generate table for preprocess parameters (for inspection only)
table.dwi.param <- map_dfr(
  list.files(path.nifti),
  function(id) {
    vol.num = readr::read_lines(str_c(path.nifti, '/', id, '/DTI/raw/data.bval')) %>% trimws() %>% strsplit(' ') %>% unlist() %>% length()
    para.parsed <- jsonlite::fromJSON(str_c(path.nifti, '/', id, '/DTI/raw/data.json'))
    return(tibble(
      path.id = str_c(path.nifti, '/', id),
      vol.num, 
      pe = ifelse(!is_null(para.parsed$PhaseEncodingDirection), para.parsed$PhaseEncodingDirection, 'j-'), # Direction of phase encoding
      rt =  ifelse(!is_null(para.parsed$TotalReadoutTime), para.parsed$TotalReadoutTime, para.parsed$EstimatedTotalReadoutTime) # Readout time
    ))
  }
)

# Check b0 image and b0 mask
slices_dir_pair(path.input = path.nifti, imageBase = 'DTI/raw/nodif.nii.gz', imageLine = 'DTI/raw/nodif_brain_mask.nii.gz', path.report =  str_c(path.QC, '/nodif_d'))

# Chek T1 (stripped)
slices_dir_single(path.input = path.nifti, image = 'synb0/INPUTS/T1.nii.gz', path.report =  str_c(path.QC, '/T1_stripped'))

# Run synb0
synb0_dir(path.input = path.nifti, process.num = 3)
slices_dir_single(path.input = path.nifti, image = 'synb0/OUTPUTS/b0_u.nii.gz', path.report =  str_c(path.QC, '/synb0'))

# Run eddy
eddy_dir(path.input = path.nifti, process.num = 5)

# Check corrected diffusion image
slices_dir_single(path.input = path.nifti, image = 'DTI/data_ud.nii.gz', path.report =  str_c(path.QC, '/nodif_ud'))
slices_dir_single(path.input = path.nifti, image = 'DTI/nodif_ud_syn.nii.gz', path.report =  str_c(path.QC, '/nodif_ud_syn'))

# Find out participants failed the pipeline
preprocess.succeed <- map_dfr(
  list.files(path.nifti),
  function(id){
    tibble(
      id,
      done = file.exists(str_c(path.nifti, '/', id, '/DTI/data_ud.nii.gz'))
    )
  }
)


