library(raster)
library(parallel)
library(abind)
library(ncdf4)
library(phenex)
library(zoo)
library(terra)


# get the smooth result
getBiseVI <- function(vi){
  vi[vi <= 0] <- NA
  vi_fill <- na.approx(vi, na.rm = FALSE)
  vi_fill[is.na(vi_fill)] <- 0
  vi_sg <- modelNDVI(vi_fill, year.int=1995, multipleSeasons=FALSE,
                     correction="ravg", window.ravg = 5, 
                     method="SavGol", filter = 5, degree = 3)
  vi_smooth <- vi_sg[[1]]@correctedValues
  return(vi_smooth)
}


tiles <- c('h23v04','h23v05','h24v04','h24v05','h25v03',
           'h25v04','h25v05','h25v06','h26v03','h26v04',
           'h26v05','h26v06','h27v04','h27v05','h27v06',
           'h28v05','h28v06','h28v07','h29v05','h29v06')

days <- c('001','009','017','025','033','041','049','057','065',
          '073','081','089','097','105','113','121','129','137',
          '145','153','161','169','177','185','193','201','209',
          '217','225','233','241','249','257','265','273','281',
          '289','297','305','313','321','329','337','345','353','361')


cl.cores <- detectCores()
cl<-makeCluster(4)
clusterExport(cl,"na.approx",envir = environment())
clusterExport(cl,"modelNDVI",envir = environment())

vi_path <- "E:/00_modis/7MOD09A1_band/EVI_median/"
out_path = "E:/00_modis/7MOD09A1_band/EVI_median_smooth/MOD09A1.EVI_Median."

for(tid in 1:5){
  tmp_pat <- paste0(tiles[tid], ".tif")
  vi_file<-list.files(vi_path, pattern = tmp_pat)
  viArray <- array(0, dim = c(50, 50, 48))
  for(vi_idx in 1:48){
    if(vi_idx == 1){
      tmp_vi <- raster(paste0(vi_path, vi_file[46]))
      tmp_vi[is.na(tmp_vi)] <- 0
      vi_in <- tmp_vi[1:50, 1:50]
    }
    if (vi_idx == 48){
      tmp_vi <- raster(paste0(vi_path, vi_file[1]))
      tmp_vi[is.na(tmp_vi)] <- 0
      vi_in <- tmp_vi[1:50, 1:50]
    }
    if ((vi_idx > 1) & (vi_idx < 48)){
      tmp_vi <- raster(paste0(vi_path, vi_file[(vi_idx-1)]))
      tmp_vi[is.na(tmp_vi)] <- 0
      vi_in <- tmp_vi[1:50, 1:50]
    }
    #tmp_vi_array <- as.array(tmp_vi)
    tmp_vi_array <- as.array(vi_in)
    viArray[, , vi_idx] <- tmp_vi_array
  }
  viArray_sm <- parApply(cl, viArray, c(1,2), getBiseVI)
  
  tmp_ras <- rast(paste0(vi_path, vi_file[1]))
  tmp_ext <- ext(tmp_ras)
  tmp_crs <- crs(tmp_ras)
  
  for(i in 1:46){
    tmp_vism <- viArray_sm[(i+1), ,]
    tmp_out <- rast(tmp_vism, crs=tmp_crs, extent = tmp_ext)
    out_name <- paste0(out_path, tiles[tid], '.', days[i], '.sg.tif')
    writeRaster(tmp_out, out_name, overwrite=TRUE)
  }
}

stopCluster(cl)
