# Created by Haojie Chen, Beijing Normal University, 20240104
library(magrittr)
library(stringr)
library(purrr)
library(dplyr)
library(glue)
#设置路径
path.nifti <<- '/home/xjx/sharefolder/IMAGE' # Absolute path of nifti
path.QC <<- '/home/xjx/sharefolder/QC' # Absolute path of quality control

# Related functions -------------
#当路径不存在时创建文件夹
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}
#删除特定的文件
remove_images <- function(path.image, pattern, process.num = 10) {
  #并行处理专用前置语言
  future::plan(future::multisession, workers = process.num)
  #并行括号内内容
  furrr::future_walk(list.files(path.image), function(ID){
    path_image = glue('{path.image}/{ID}/{pattern}')
    system(glue('rm {path_image}'))
  })
}
#根据表格进行FOV图像视野调整
fov_with_table <- function(fov_table, process.num = 5) {
  # fov_table contains: path.id（subjectID）, z（FOV参数） and adjust（0不调整，1调整）
  # path.id must be absolute directories
  future::plan(future::multisession, workers = process.num)
  furrr::future_pwalk(
    list(fov_table$path.id, fov_table$z, fov_table$adjust),
    function(path_id, z, adjust) {
      if(dir.exists(path_id)) {
        system(glue('cd {path_id}'))
        setwd(path_id)
        if(!file.exists('T1/data.nii')) return()
        # FOV adjustment
        if(adjust == 1) {
          system(glue('robustfov -i T1/data.nii -r T1/T1_FOV.nii.gz -b {z} > T1/fov_init'))
        }
      }
    }
  )
}
#根据表格进行BET剥头皮
bet_with_table <- function(bet_table, process.num = 10){
  # bet_table contains: path.id and f（拨头皮参数）
  # path.id must be absolute directories
  future::plan(future::multisession, workers = process.num)
  furrr::future_pwalk(
    list(bet_table$path.id, bet_table$f, bet_table$adjust),
    function(path_id, thresh_f, adjust) {
      if (adjust == 0) return()
      if(dir.exists(path_id)) {
        system(glue('cd {path_id}'))
        setwd(path_id)
        if(!file.exists('T1/data.nii')) return()
        if(!file.exists('T1/T1_FOV.nii.gz')) return()
        # Brain extraction with given fractional intensity threshold (f)
        system(
          glue('bet T1/T1_FOV.nii.gz T1/T1_brain.nii.gz -f {thresh_f} -B -o -m -s')
        )
      }
    }
  )
}
#生成html的质量控制的报告
slices_dir_single <- function(path.input, image, path.report) {
  path.report <- create_when_absent(path.report)
  system(glue('cd {path.report}'))
  setwd(path.report)
  image.list <- str_c(path.input, '/', list.files(path.input), '/', image)
  image.list <- paste(image.list, collapse = ' ')
  #slicedir时多视图切边，便于影像质量查看
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

# Removes images by pattern ONLY when REDO---------------------------------------

remove_images(path.nifti, pattern = '/T1/t1_*')
remove_images(path.nifti, pattern = '/T1/fov*')


# QC: Existence table --------------------------
path.id <- str_c(path.nifti, '/', list.files(path.nifti))
exist.id <- map_int(path.id, ~if_else(file.exists( str_c(.,'/T1/data.nii') ), 1, 0))
table_exist <- tibble(
  path.id = path.id
) %>% mutate(
  exist.t1 = map_int(path.id, ~if_else(file.exists( str_c(.,'/T1/data.nii') ), 1, 0)),
  exist.flair = map_int(path.id, ~if_else(file.exists( str_c(.,'/FLAIR/data.nii') ), 1, 0)),
  exist.dti = map_int(path.id, ~if_else(file.exists( str_c(.,'/DTI/data.nii') ), 1, 0)),
  exist.bold = map_int(path.id, ~if_else(file.exists( str_c(.,'/BOLD/data.nii') ), 1, 0))
)
readr::write_csv(table_exist, str_c(path.QC, '/exist_all_modal.csv'))


# FOV adjustment --------------------------------------------------

# Generate initial FOV table
fov_table_t1 <- readr::read_csv(str_c(path.QC, '/exist_all_modal.csv')) %>% 
  filter(exist.t1 == 1) %>% 
  select(path.id) %>% 
  mutate(z = 190) %>% # By default z is 170. The Larger, the lower, the more regions
  mutate(adjust = 1)
readr::write_csv(fov_table_t1, glue('{path.QC}/fov_table_z.csv'))

# MANUALLY and recursively do the FOV adjustment
round = 1 # PLUS one each re-run
fov_table_t1 <- readr::read_csv(glue('{path.QC}/fov_table_z.csv'))
num.adjust <- nrow(fov_table_t1 %>% filter(adjust == 1))
if(num.adjust > 1) {
  # Re-FOV adjust
  fov_with_table(fov_table = fov_table_t1, process.num = 10)
  fov_table_t1 <- mutate(fov_table_t1, adjust = 0)
  readr::write_csv(fov_table_t1, glue('{path.QC}/fov_table_z.csv'))
  slices_dir_single(path.nifti, 'T1/T1_FOV.nii.gz', glue('{path.QC}/T1_FOV_{round}'))
  
  # OURSIDE R: Change z accordingly based on FOV results
}



# BET: Brain extraction with initial fractional intensity threshold ---------------

bet_table_t1 <- readr::read_csv(str_c(path.QC, '/exist_all_modal.csv')) %>% 
  filter(exist.t1 == 1) %>% 
  select(path.id) %>% 
  mutate(f = 0.18) %>%  # Initial thresh
  mutate(adjust = 1)

readr::write_csv(bet_table_t1, glue('{path.QC}/bet_table_f.csv'))

# MANUALLY and recursively do the FOV adjustment
round = 1 # PLUS one each re-run
bet_table_t1 <- readr::read_csv(glue('{path.QC}/bet_table_f.csv'))
num.adjust <- nrow(bet_table_t1 %>% filter(adjust == 1))
if(num.adjust > 0) {
  # Re-FOV adjust
  bet_with_table(bet_table = bet_table_t1, process.num = 5)
  bet_table_t1 <- mutate(bet_table_t1, adjust = 0)
  readr::write_csv(bet_table_t1, glue('{path.QC}/bet_table_f.csv'))
  slices_dir_pair(path.input = path.nifti, imageBase = 'T1/T1_brain_overlay.nii.gz', imageLine = 'T1/T1_brain_skull.nii.gz', path.report = glue('{path.QC}/T1_BET_{round}'))
  
  # OURSIDE R: Change z accordingly based on FOV results
}

# QC: Image quality of raw nifty with slices dir ------------------
slices_dir_single(
  path.input = path.nifti, 
  image = 'T1/data.nii',
  path.report = '/home/xjx/sharefolder/QC/t1_raw'
)
slices_dir_single(
  path.input = path.nifti, 
  image = 'FLAIR/data.nii',
  path.report = '/home/xjx/sharefolder/QC/flair_raw'
)
slices_dir_single(
  path.input = path.nifti, 
  image = 'DTI/data.nii',
  path.report = '/home/xjx/sharefolder/QC/DTI_raw'
)
slices_dir_single(
  path.input = path.nifti, 
  image = 'BOLD/data.nii',
  path.report = '/home/xjx/sharefolder/QC/bold_raw'
)