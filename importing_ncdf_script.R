# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths(), libpathKiri))
library(tidyverse)
library(rasterVis)
library(ncdf4)

indir <- "/net/pc150398/nobackup_1/users/meirink/CAMS_for_Kilian/"
filedirec <- "/nobackup/users/bakker/Data2/particlesvariables/"
bounds <- c(3.25, 7.25, 50.75, 53.5)
leadTimes <- seq(0,48,3)

#function that imports the raster:
import_ncdf <- function(afile, time, var){
  tmpr <- raster(afile, band = time, varname = var)
  return(tmpr)
}

LeadTimeImport <- function(band, file){
  tmprast1 <- import_ncdf(paste0(indir,file), band, "aot500")
  tmpcroprast1 <- crop(tmprast1, extent(bounds))
  tmprast2 <- import_ncdf(paste0(indir,file), band, "ang_exp")
  tmpcroprast2 <- crop(tmprast2, extent(bounds))
  tmprast3 <- import_ncdf(paste0(indir,file), band, "o3col")
  tmpcroprast3 <- crop(tmprast3, extent(bounds))
  tmpdf <- data.frame(lt = leadTimes[band], xcoor = coordinates(tmpcroprast1)[,1], ycoor = coordinates(tmpcroprast1)[,2], 
                      AOD_500 = values(tmpcroprast1), Ang_exp = values(tmpcroprast2), Ozone = values(tmpcroprast3))
  return(tmpdf)
}

savingData <- function(file){
   for (j in 1:length(leadTimes)){
     tmpdf <- LeadTimeImport(j,file)
     if (j == 1){
       tmpdf2 <- tmpdf
     } else {
       tmpdf2 <- rbind(tmpdf2,tmpdf)
     }
   }

  date <- gsub("0000.nc", "", gsub("CAMS_", "", basename(file)))
  tempDate <- paste0(substr(date,1,4), "-", substr(date,5,6), "-", substr(date,7,8))
  date2 <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
  saveRDS(tmpdf2, file=paste0(filedirec,date2,".rds"))
}

fnames    <- data.frame(files = list.files(path = indir, full.names = TRUE), stringsAsFactors = FALSE)
init_files  <- basename(fnames$files)[-length(fnames$files)]
for (i in 1:length(init_files)){
  print(i)
  Sys.sleep(0.001)
  savingData(init_files[i])
}