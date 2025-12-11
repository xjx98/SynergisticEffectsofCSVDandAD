# Created by Haojie Chen, Beijing Normal University, 20240116
library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(glue)
# Global variables ----------------------设置路径
path.nifti <<- '/home/xjx/sharefolder/IMAGE' # Absolute path of nifti
path.QC <<- '/home/xjx/sharefolder/QC' # Absolute path of quality control
path.IDP <<- '/home/xjx/sharefolder/IDP' # Imaging-derivative phenotype
PATH_MNI152_2mm_BRAIN <<- '/home/xjx/soft/fsl_templates/fsl_templates/MNI152_T1_2mm_brain.nii.gz'
PATH_MNI152_2mm_SKULL <<- '/home/xjx/soft/fsl_templates/fsl_templates/MNI152_T1_2mm_skull.nii.gz'
PATH_MNI152_2mm <<- '/home/xjx/soft/fsl_templates/fsl_templates/MNI152_T1_2mm.nii.gz'
PATH_MNI152_1mm <<- '/home/xjx/soft/fsl_templates/fsl_templates/MNI152_T1_1mm.nii.gz'
# Related functions -------------------------------------------------------
#如果路径不存在则新建
create_when_absent <<- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}
#获取文件扩展名
get_file_extension <- function(file_path) {
  # Split the file name by the dot and extract the last part
  parts <- strsplit(file_path, "\\.")[[1]]
  if (length(parts) > 1) {
    return(tail(parts, n=1))
  } else {
    return('')  # Return empty string if there is no file extension
  }
}

# Drop the file extension (e.g. X.nii.gz -> X)移除文件扩展名
drop_extension <- function(path) {
  if(!grepl(pattern = '\\..[^\\.]*$', path)) return(path)
  drop_extension(sub('\\..[^\\.]*$', '', path)) # Recursively removes all file extension (e.g. .nii.gz)
}
#删除特定文件
remove_images <- function(path.image, pattern, process.num = 10) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.image), function(ID){
    path_image = glue('{path.image}/{ID}/{pattern}')
    system(glue('rm {path_image}'))
  })
}
#拷贝文件
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

# Reduced SEINAX
# Reference:  https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/blob/master/bb_structural_pipeline/bb_sienax
MNI152_2mm_BRAIN = drop_extension(PATH_MNI152_2mm_BRAIN)
MNI152_2mm_SKULL = drop_extension(PATH_MNI152_2mm_SKULL)
MNI152_2mm = drop_extension(PATH_MNI152_2mm)
MNI152_1mm = drop_extension(PATH_MNI152_1mm)
#执行线性配准，计算体积缩放因子并应用变换
sienax.dir <- function(path.input, process.num = 5) {
  # print(MNI152_2mm) # OK
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    create_when_absent('transform')
    if(file.exists('transform/vscale.txt')) return()
    
    # Linear registration of T1w to std space
    if(!file.exists('transform/T1_to_MNI_linear.mat')) {
	#将T1配准到标准空间
      system(str_c("pairreg ", MNI152_2mm_BRAIN, " T1/T1_brain ", MNI152_2mm_SKULL, " T1/T1_brain_skull transform/T1_to_MNI_linear.mat"))
    }
    
    # Calculates volume scaling factor计算缩放因子
    if(!file.exists('transform/T1_to_MNI_linear.avscale')) {
      #system(glue('avscale transform/T1_to_MNI_linear.mat {MNI152_2mm} > transform/T1_to_MNI_linear.avscale'))
      system(str_c("avscale transform/T1_to_MNI_linear.mat ", MNI152_2mm, " > transform/T1_to_MNI_linear.avscale"))
      if (!file.exists('transform/vscale.txt')) {
        xscale <- system('grep Scales transform/T1_to_MNI_linear.avscale | awk \'{print $4}\'', intern = T) %>% as.numeric()
        yscale <- system('grep Scales transform/T1_to_MNI_linear.avscale | awk \'{print $5}\'', intern = T) %>% as.numeric()
        zscale <- system('grep Scales transform/T1_to_MNI_linear.avscale | awk \'{print $6}\'', intern = T) %>% as.numeric()
        vscale <- xscale * yscale * zscale
        readr::write_lines(vscale, 'transform/vscale.txt')
      }
    }
    
    # Apply transformations应用变换
    # The results will be undergone QC and combined with FAST output
    if(!file.exists('transform/T1_to_MNI_linear.nii.gz')){
      #system(glue('flirt -in T1/T1_FOV  -ref {MNI152_1mm} -o transform/T1_to_MNI_linear -applyxfm -init transform/T1_to_MNI_linear.mat -interp spline'))
      system(str_c("flirt -in T1/T1_FOV  -ref ", MNI152_1mm, " -o transform/T1_to_MNI_linear -applyxfm -init transform/T1_to_MNI_linear.mat -interp spline"))
    }
    if(!file.exists('transform/T1_brain_skull_to_MNI_linear.nii.gz')) {
      #system(glue('flirt -in T1/T1_brain_skull  -ref {MNI152_1mm} -o transform/T1_brain_skull_to_MNI_linear -applyxfm -init transform/T1_to_MNI_linear.mat -interp trilinear')) 
      system(str_c("flirt -in T1/T1_brain_skull  -ref ", MNI152_1mm, " -o transform/T1_brain_skull_to_MNI_linear -applyxfm -init transform/T1_to_MNI_linear.mat -interp trilinear"))
    }
  })
}
#脑组织分割并制作掩膜，偏置场校正
fast.dir <- function(path.input, process.num = 3) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    create_when_absent('T1/T1_fast')
    if(file.exists('T1/T1_unbiased_brain.nii.gz')) return()
    # Run fast (180s)脑组织分割乘白质灰质脑脊液
    if(!file.exists('T1/T1_fast/T1_brain_bias.nii.gz')) {
      system('fast -b -o T1/T1_fast/T1_brain T1/T1_brain')
    }
 
    # Binarize PVE masks (10s)分割后的图像制作掩膜，0脑脊液1灰质2白质
    if(!file.exists('T1/T1_fast/T1_brain_WM_mask.nii.gz')) {
      system('fslmaths T1/T1_fast/T1_brain_pve_0.nii.gz -thr 0.5 -bin T1/T1_fast/T1_brain_CSF_mask.nii.gz')
      system('fslmaths T1/T1_fast/T1_brain_pve_1.nii.gz -thr 0.5 -bin T1/T1_fast/T1_brain_GM_mask.nii.gz')
      system('fslmaths T1/T1_fast/T1_brain_pve_2.nii.gz -thr 0.5 -bin T1/T1_fast/T1_brain_WM_mask.nii.gz')
    }

    # Apply bias field correction to T1 (10s)偏置场校正
    if(!file.exists('T1/T1_unbiased_brain.nii.gz')) {
      system('fslmaths T1/T1_FOV.nii.gz -div T1/T1_fast/T1_brain_bias.nii.gz T1/T1_unbiased.nii.gz')
      system('fslmaths T1/T1_brain.nii.gz -div T1/T1_fast/T1_brain_bias.nii.gz T1/T1_unbiased_brain.nii.gz')
    }
  })
}
#脑区分割并计算体积等指标，并生成result.txt文件
first.dir <- function(path.input, process.num = 3) {
  future::plan(future::multisession, workers = process.num)
  furrr::future_walk(list.files(path.input), function(ID){
    path_id <- glue('{path.input}/{ID}')
    system(glue('cd {path_id}'))
    setwd(path_id)
    create_when_absent('T1/T1_first')
    if(file.exists('T1/T1_first/result.txt')) return()
    if(!file.exists('T1/T1_first/T1_first_all_fast_firstseg.nii.gz')) {
	#脑区分割并计算体积等指标
      system('run_first_all -i T1/T1_unbiased_brain.nii.gz -b -o T1/T1_first/T1_first')
    }
    if(!file.exists('T1/T1_first/result.txt')) {
      #整理结果乘txt，system('fslstats T1/T1_first/T1_first_all_fast_firstseg -H 58 0.5 58.5 | sed /'s/\.000000//g/' | awk /'BEGIN { ORS = " " } { print }/'| awk /'{print $10 " " $49 " " $11 " " $50 " " $12 " " $51 " " $13 " " $52 " " $17 " " $53 " " $18 " " $54 " " $26 " " $58 " " $16 }' > T1/T1_first/result.txt')
      system("fslstats T1/T1_first/T1_first_all_fast_firstseg -H 58 0.5 58.5 | sed 's/\\.000000//g' | awk 'BEGIN { ORS = \" \" } { print }' | awk '{print $10 \" \" $49 \" \" $11 \" \" $50 \" \" $12 \" \" $51 \" \" $13 \" \" $52 \" \" $17 \" \" $53 \" \" $18 \" \" $54 \" \" $26 \" \" $58 \" \" $16 }' > T1/T1_first/result.txt")
    }
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
  image.list.B <- str_c(path.input, '/', list.files(path.input), '/', imageLine)
  image.list <- purrr::map_chr(image.list.B, function(x) stringr::str_c(x, ' ', imageBase))
  image.list <- paste(image.list, collapse = ' ')
  system(glue('slicesdir -o {image.list}'))
}

# Run reduced SIENAX ----------------------
sienax.dir(path.nifti, process.num = 10)

# QC of transformations of linear transformation
slices_dir_pair_fixed_line(path.input = path.nifti, 
                           imageBase = PATH_MNI152_1mm, 
                           imageLine = 'transform/T1_to_MNI_linear.nii.gz', 
                           path.report = glue('{path.QC}/sienax_T1_to_MNI_linear'))


# Run fast (5 min for each image)---------------------------------
fast.dir(path.input = path.nifti, process.num = 5)
slices_dir_single(path.input = path.nifti, image = 'T1/T1_fast/T1_brain_CSF_mask.nii.gz', path.report = str_c(path.QC, '/FAST_CSF')) # Checks CSF
slices_dir_single(path.input = path.nifti, image = 'T1/T1_fast/T1_brain_GM_mask.nii.gz', path.report = str_c(path.QC, '/FAST_GM')) # Checks GM
slices_dir_single(path.input = path.nifti, image = 'T1/T1_fast/T1_brain_WM_mask.nii.gz', path.report = str_c(path.QC, '/FAST_WM')) # Checks WM

# Run first (5min for each image) ---------------------
first.dir(path.input = path.nifti, process.num = 5)
slices_dir_pair(path.input = path.nifti, imageBase = 'T1/T1_unbiased_brain.nii.gz', imageLine = 'T1/T1_first/T1_first_all_fast_firstseg.nii.gz', path.report = str_c(path.QC, '/first_seg_subcortical'))

# Copy Related IDPs ----------------
copy_IDP(group_label = 'group_test', path.IDP.input = 'transform/vscale.txt', IDP_label = 'vscaling', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Volumetric scaling
copy_IDP(group_label = 'group_test', path.IDP.input = 'T1/T1_first/result.txt', IDP_label = 'FIRST', path.input = path.nifti, path.IDP.output = path.IDP, process.num = 10) # Subcortical segmentation
