library(terra)

fname <- "E:/00_GLASS/1LAI_AVHRR/GLASS01B02.V40.A2018361.2019358.hdf"



avhrr_path <- "E:/00_GLASS/1LAI_AVHRR/"
out_path <- "E:/00_GLASS/1LAI_China/"

file_list <- list.files(avhrr_path)

china_ext <- ext(72, 136, 17, 54)

for(fid in 599:length(file_list)){
  tmp_name <- paste0(avhrr_path, file_list[fid])
  tmp_ras <- rast(tmp_name)
  tmp_ras_cn <- crop(tmp_ras, china_ext)
  avhrr_name<- substr(file_list[fid], 1, 23)
  out_name <- paste0(out_path, avhrr_name, ".tif")
  writeRaster(tmp_ras_cn, out_name, overwrite=TRUE)
}



plot(tmp_ras_cn)
