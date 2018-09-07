# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths(), libpathKiri))
library(tidyverse)
library(rasterVis)
library(ncdf4)

# The variables/constants
yearNumber     <- 2016
monthNumber    <- 04
dayNumber      <- 01
timeNumber     <- 00
variableNumber <- 265

dateNumber     <- yearNumber*10000 + monthNumber*100 + dayNumber
indir          <- "/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS"

# write a function that imports the raster:
import_ncdf <- function(afile, time, var){
  tmpr <- raster(afile, band = time, varname = var)
  return(tmpr)
}

# import the rasters and binds and saves them
savingData <- function(files){
    band <- 1
    tmprast1 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "aot500")
    tmpcroprast1 <- crop(tmprast1, extent(3.4, 7.1, 50.8, 53.5))
    tmprast2 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "ang_exp")
    tmpcroprast2 <- crop(tmprast2, extent(3.4, 7.1, 50.8, 53.5))
    tmprast3 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "o3col")
    tmpcroprast3 <- crop(tmprast3, extent(3.4, 7.1, 50.8, 53.5))
    tmpdf <- data.frame(lt = (band-1)*3, xcoor = coordinates(tmpcroprast1)[,1], ycoor = coordinates(tmpcroprast1)[,2], 
                        AOT_500 = values(tmpcroprast1), Ang_exp = values(tmpcroprast2), Ozone = values(tmpcroprast3))
  for (band in 2:9){
    tmprast1 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "aot500")
    tmpcroprast1 <- crop(tmprast1, extent(3.4, 7.1, 50.8, 53.5))
    tmprast2 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "ang_exp")
    tmpcroprast2 <- crop(tmprast2, extent(3.4, 7.1, 50.8, 53.5))
    tmprast3 <- import_ncdf(paste0("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/",files), band, "o3col")
    tmpcroprast3 <- crop(tmprast3, extent(3.4, 7.1, 50.8, 53.5))
    tmpdf2 <- data.frame(lt = (band-1)*3, xcoor = coordinates(tmpcroprast1)[,1], ycoor = coordinates(tmpcroprast1)[,2], 
                         AOT_500 = values(tmpcroprast1), Ang_exp = values(tmpcroprast2), Ozone = values(tmpcroprast3))
    tmpdf <- rbind(tmpdf, tmpdf2)
  }

  tmpfile <- tmpdf
  date <-gsub(".nc", "", gsub("CAMS_", "", basename(files)))
  saveRDS(tmpfile, file=paste0("/nobackup/users/bakker/Data/particlesvariables/",date,".rds"))
}

fnames    <- data.frame(files = list.files(path = indir, full.names = TRUE), stringsAsFactors = FALSE)
init_files  <- basename(fnames$files)[33:865]
init_files <- init_files[!init_files %in% c("CAMS_20171016_fc0-24h.nc","CAMS_20171017_fc0-24h.nc")]
for (i in 1:length(init_files)){
  savingData(init_files[i])
}