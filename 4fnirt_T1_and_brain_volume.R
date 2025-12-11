# Created by Haojie Chen, Beijing Normal University, 20240123
library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(glue)
# Global variables ----------------------
path.nifti <<- '/home/xjx/sharefolder/IMAGE/' # Absolute path of nifti
path.QC <<- '/home/xjx/sharefolder/QC' # Absolute path of quality control
path.IDP <<- '/home/xjx/sharefolder/IDP' # Imaging-derivative phenotype

PATH_MNI152_1mm <<- '/home/xjx/soft/fsl_templates/fsl_templates/MNI152_T1_1mm.nii.gz'

PATH_FNIRT_MASK_1mm <<- '/home/xjx/soft/bb_data/MNI152_T1_1mm_brain_mask_dil_GD7.nii.gz'
PATH_FNIRT_CONFIG <<- '/home/xjx/soft/bb_data/bb_fnirt.cnf'

PATH_MNI152_2mm_STRUCSEG_PERIPH <<- '/home/xjx/soft/bb_data/MNI152_T1_2mm_strucseg_periph.nii.gz'
PATH_MNI152_T1_SEGVENT <<- '/home/xjx/soft/bb_data/T1_segvent.nii.gz'

PATH_ATLAS_HO <<- '/home/xjx/soft/atlases/UKB_HO_atlas.nii.gz'
PATH_ATLAS_MSA_S1 <<- '/home/xjx/soft/atlases/Tian_Subcortex_S1_3T_1mm.nii.gz'
PATH_ATLAS_MSA_S2 <<- '/home/xjx/soft/atlases/Tian_Subcortex_S2_3T_1mm.nii.gz'
PATH_ATLAS_MSA_S3 <<- '/home/xjx/soft/atlases/Tian_Subcortex_S3_3T_1mm.nii.gz'
PATH_ATLAS_MSA_S4 <<- '/home/xjx/soft/atlases/Tian_Subcortex_S4_3T_1mm.nii.gz'
# Related functions of file operations and QC-------------------------------------------------------
create_when_absent <<- function(path) {
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

remove_images <- function(path.image, pattern, dir = F, process.num = 10) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.image), function(ID){
    path_image = glue('{path.image}/{ID}/{pattern}')
    if(!dir) {
      system(glue('rm {path_image}'))
    } else (
      system(glue('rm -r {path_image}'))
    )
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

# QC with slices dir
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

slices_dir_pair_fixed_base <- function(path.input, imageBase, imageLine, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list.B <- str_c(path.input, '/', list.files(path.input), '/', imageLine)
  image.list <- purrr::map_chr(image.list.B, function(x) stringr::str_c(imageBase, ' ', x))
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir -o {image.list}'))
}

slices_dir_pair_fixed_line <- function(path.input, imageBase, imageLine, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list.B <- str_c(path.input, '/', list.files(path.input), '/', imageBase)
  image.list <- purrr::map_chr(image.list.B, function(x) stringr::str_c(x, ' ', imageLine))
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir -o {image.list}'))
}

# Functions of image processing ----------------------
MNI152_1mm = drop_extension(PATH_MNI152_1mm)
#非线性配准
t1_mni_non_linear_dir <- function(path.input, process.num = 3) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    if(file.exists('transform/T1_to_MNI_warp_coef_inv.nii.gz')) return()
    
    # 非线性配准First create warp: T1 to MNI (40' for each participant)
    if(!file.exists('transform/T1_to_MNI_warp_coef.nii.gz')) { 
      command.temp <- str_c(
        'fnirt --in=T1/T1_unbiased.nii.gz --ref=', MNI152_1mm, ' --aff=transform/T1_to_MNI_linear.mat --config=', PATH_FNIRT_CONFIG, 
        ' --refmask=', PATH_FNIRT_MASK_1mm,
        ' --cout=transform/T1_to_MNI_warp_coef --fout=transform/T1_to_MNI_warp --jout=transform/T1_to_MNI_warp_jac --iout=transform/T1_in_MNI.nii.gz --interp=spline'
      )
      system(command.temp)
    }
    # 配准需要用到的反向转换Invert the warp to create the warp from MNI to T1
    if(!file.exists('transform/T1_to_MNI_warp_coef_inv.nii.gz')) {
      system('invwarp --ref=T1/T1_unbiased.nii.gz -w transform/T1_to_MNI_warp_coef -o transform/T1_to_MNI_warp_coef_inv')
      # system('invwarp --ref=T1/T1_unbiased.nii.gz -w transform/T1_to_MNI_warp -o transform/T1_to_MNI_warp_inv') # Run when necessary
    }
  })
}

# 脑区体积分割Gets the volume for single segmentation
# Units: mm3
get_seg_volume <- function(path.image.seg) {
  S = system(str_c('fslstats ', path.image.seg, ' -m -v'), intern = T)
  values <- strsplit(S, ' ')[[1]]
  xa <- as.numeric(values[1])
  xb <- as.numeric(values[3])
  xa*xb
}

# Continuing with sienax, non-linear warps the FAST segmentation and calculates the volumes线性配准、非线性配准和计算体积
fast_warp_volume <- function(path.input, process.num = 5) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    if(file.exists('T1/T1_fast/volumes.csv')) return()
    
    # Transform peripheral GM
    if(!file.exists('T1/T1_fast/T1_segperiph.nii.gz')) {
	#非线性变换
      system(str_c('applywarp --rel --interp=trilinear --in=', PATH_MNI152_2mm_STRUCSEG_PERIPH, '  --ref=T1/T1_unbiased.nii.gz -w transform/T1_to_MNI_warp_coef_inv -o T1/T1_fast/T1_segperiph'))
      system('fslmaths T1/T1_fast/T1_segperiph -thr 0.5 T1/T1_fast/T1_segperiph')
    }

    # Transform ventricular CSF
    if(!file.exists('T1/T1_fast/T1_segvent.nii.gz')) {
      system(str_c('applywarp --rel --interp=nn --in=', drop_extension(PATH_MNI152_T1_SEGVENT), ' --ref=T1/T1_unbiased.nii.gz -w transform/T1_to_MNI_warp_coef_inv -o T1/T1_fast/T1_segvent'))
    }

    # Apply masks
    if(!file.exists('T1/T1_fast/T1_pve_1_segperiph.nii.gz')) {
      system('fslmaths T1/T1_fast/T1_brain_pve_1 -mas T1/T1_fast/T1_segperiph T1/T1_fast/T1_pve_1_segperiph -odt float') # Peripheral GM
    }
    
    if(!file.exists('T1/T1_fast/T1_pve_0_segvent.nii.gz')) {
      system('fslmaths T1/T1_fast/T1_brain_pve_0 -mas T1/T1_fast/T1_segvent T1/T1_fast/T1_pve_0_segvent -odt float') # Ventricular CSF
    }
    
    # Generates the volumes
    if(!file.exists('T1/T1_fast/volumes.csv')) {
      result = map2_dfr(
        list('T1/T1_fast/T1_pve_1_segperiph.nii.gz', 'T1/T1_fast/T1_pve_0_segvent.nii.gz', 'T1/T1_fast/T1_brain_pve_1.nii.gz', 'T1/T1_fast/T1_brain_pve_2.nii.gz', 'T1/T1_fast/T1_brain_pve_0.nii.gz'),
        list('GMVperi', 'CSFVvent', 'GMV', 'WMV', 'CSFV'),
        function(image_seg, image_label){
          image_seg = str_c(path_id, '/', image_seg)
          tibble(value = get_seg_volume(image_seg), label = image_label)
        }
      ) %>% 
        pivot_wider(names_from = 'label', values_from = 'value') %>% 
        mutate(GMVsub = GMV - GMVperi) %>% 
        mutate(BV = GMV + WMV) %>% 
        mutate(ICV = GMV + WMV + CSFV) %>% 
        select(GMV, WMV, CSFVvent, BV, GMVperi, GMVsub, CSFV, ICV)
      readr::write_csv(result, 'T1/T1_fast/volumes.csv')
    }
  })
}

# Warp atlas and gets the volume for brain parcellations用模板获取的脑区体积
# Units: mm3
warp_atlas_and_volume <- function(path.input, list.path.atlas, list.path.result, process.num = 8) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    create_when_absent('T1/atlases')
    walk2(
      list.path.atlas, list.path.result,
      function(path.atlas, path.atlas.to.t1) {
        #setwd(path_id)
        #system(str_c('cd ', path_id))
        path.atlas.to.t1 = drop_extension(path.atlas.to.t1)
        if(!file.exists(str_c(path.atlas.to.t1, '.nii.gz'))) {
          system(str_c('applywarp -i ', drop_extension(path.atlas), ' -o ', path.atlas.to.t1, ' -r T1/T1_unbiased -w transform/T1_to_MNI_warp_coef_inv.nii.gz --interp=nn'))
        }
        if(!file.exists(str_c(path.atlas.to.t1, '_raw.txt'))) {
          system(str_c('fslstats -K ', path.atlas.to.t1, '.nii.gz T1/T1_fast/T1_brain_pve_1.nii.gz -m -v >', path.atlas.to.t1, '_raw.txt')) # Computes volume within each region with GM segmentation as mask following UKB-pipeline
          result.raw <- read.table(str_c(path.atlas.to.t1, '_raw.txt')) %>% as_tibble()
          result <- result.raw %>% mutate(V = V1*V2)
          readr::write_lines(result %>% pull(V), str_c(path.atlas.to.t1, '.txt'))
        }
      }
    )
  })
}

# Non-linear registration (MNI & T1) --------------------------------
t1_mni_non_linear_dir(path.input = path.nifti, process.num = 3)

# QC of non-linear transformations of linear transformation
slices_dir_pair_fixed_line(path.input = path.nifti, 
                           imageBase = 'transform/T1_in_MNI.nii.gz', 
                           imageLine = PATH_MNI152_1mm, 
                           path.report = glue('{path.QC}/FNIRT_T1_to_MNI_non_linear'))

# Brain tissue volumes on FAST and SIENAX result-------------------- -------------------------------
fast_warp_volume(path.input = path.nifti, process.num = 3)
# QC of cortical GM
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/T1_fast/T1_pve_1_segperiph.nii.gz',
  path.report = str_c(path.QC, '/FAST_PERI_GM')
)

# QC of ventricular CSF
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/T1_fast/T1_pve_0_segvent.nii.gz',
  path.report = str_c(path.QC, '/FAST_VENT_CSF')
)


# Warp atlases and computes volumes ---------------------------------------
warp_atlas_and_volume(
  path.input = path.nifti, process.num = 3,
  list.path.atlas = list(PATH_ATLAS_HO, PATH_ATLAS_MSA_S1, PATH_ATLAS_MSA_S2, PATH_ATLAS_MSA_S3, PATH_ATLAS_MSA_S4), 
  list.path.result = list('T1/atlases/UKB_HO', 'T1/atlases/MSA_S1', 'T1/atlases/MSA_S2', 'T1/atlases/MSA_S3', 'T1/atlases/MSA_S4')
  )

# remove_images(path.nifti, pattern = 'T1/atlases', dir = T)
# QC of HO atlas
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/atlases/UKB_HO.nii.gz',
  path.report = str_c(path.QC, '/ATLAS_UKB_HO')
)

# QC of MSA S1
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/atlases/MSA_S1.nii.gz',
  path.report = str_c(path.QC, '/ATLAS_MSA_S1')
)

# QC of MSA S2
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/atlases/MSA_S2.nii.gz',
  path.report = str_c(path.QC, '/ATLAS_MSA_S2')
)

# QC of MSA S3
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/atlases/MSA_S3.nii.gz',
  path.report = str_c(path.QC, '/ATLAS_MSA_S3')
)

# QC of MSA S4
slices_dir_pair(
  path.input = path.nifti,
  imageBase = 'T1/T1_unbiased.nii.gz',
  imageLine = 'T1/atlases/MSA_S4.nii.gz',
  path.report = str_c(path.QC, '/ATLAS_MSA_S4')
)
# Copy Related IDPs ----------------
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/T1_fast/volumes.csv', IDP_label = 'volume_fast', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Volumetric scaling
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/atlases/UKB_HO.txt', IDP_label = 'volume_UKB_HO', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Harvard-Oxford parcellation
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/atlases/MSA_S1.txt', IDP_label = 'volume_MSA_S1', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical parcellation
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/atlases/MSA_S2.txt', IDP_label = 'volume_MSA_S2', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical parcellation
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/atlases/MSA_S3.txt', IDP_label = 'volume_MSA_S3', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical parcellation
copy_IDP(group_label = 'AD_CI', path.IDP.input = 'T1/atlases/MSA_S4.txt', IDP_label = 'volume_MSA_S4', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical parcellation
