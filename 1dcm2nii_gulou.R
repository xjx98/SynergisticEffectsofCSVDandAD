install.packages("installr") 
library("installr")  
installr()
# Created by Haojie Chen, Beijing Normal University, 20231220
#定义路径和并行处理数量
path.dcm = 'G:/programG_NC_CSVD_AD_DTI_BEIJING/image/CSVD/10'
path.nii = 'F:/programG_NC_CSVD_AD_DTI_BEIJING/imagenii/CSVD/10'
process.number = 10
path.transfer.table = 'F:/programG_NC_CSVD_AD_DTI_BEIJING/summary.csv'
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)
#定义辅助函数create_when_absent如果路径不存在则创建路径
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
  return(path)
}
#定义辅助函数dcm2nii_gulou用于DICOM转换为NIFTI
dcm2nii_gulou <- function(path.dcm, path.nii, process.number) {
  create_when_absent(path.nii)
  #使用多进程进行处理
  future::plan(future::multisession, workers = process.number)
  #遍历DICOM文件夹中的每个SUBJECT ID并创建对应的结果文件夹
  furrr::future_map_dfr(list.files(path.dcm), function(subject_id){
    path.nii.id = create_when_absent(glue('{path.nii}/{subject_id}'))
	#遍历每个subject ID 下的影像序列并创建对应的结果文件夹
    map_dfr(list.files(glue('{path.dcm}/{subject_id}')), function(image_modal){
      path.nii.id.modality = create_when_absent(glue('{path.nii.id}/{image_modal}'))
      path.dcm.id.modality = glue('{path.dcm}/{subject_id}/{image_modal}')
	  #使用dcm2niix将序列下的DICOM文件转为NIFTI
      system(glue('dcm2niix -o {path.nii.id.modality} -f data {path.dcm.id.modality}'))
	  #获取转换后的NIFTI文件的路径和大小
      path.nii.data = glue('{path.nii.id.modality}/data.nii')
      nii.size = ifelse(file.exists(path.nii.data), file.info(path.nii.data)$size / 1000^2, 0)
      #返回包含相关信息的tibble
	  return(tibble(
        id = subject_id,
        path_dcm = path.dcm.id.modality,
        path_nii = path.nii.data,
        modal = image_modal,
        num_volume = length(list.files(path.dcm.id.modality)),
        nii_size = nii.size
      ))
    })
  })
}
#运行函数并将结果保存在csv文件中
summary_transfer = dcm2nii_gulou(path.dcm, path.nii, process.number = process.number)
readr::write_csv(summary_transfer, path.transfer.table)
