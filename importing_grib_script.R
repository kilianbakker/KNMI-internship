# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths(), libpathKiri))
library(tidyverse)
library(Rgrib2)
library(rasterVis)

# The variables/constants
yearNumber     <- 2016
monthNumber    <- 04
dayNumber      <- 01
timeNumber     <- 00
variableNumber <- 265

dateNumber     <- yearNumber*10000 + monthNumber*100 + dayNumber
indir          <- "/net/pc150388/nobackup_3/users/schmeits/HARM38/"

heights <- c(30000, 25000, 21000, 18437, 16838, 15557, 14486, 13564, 12753, 12026, 11365, 10756, 10182, 9638, 9122, 8631, 8163, 7716, 7289, 6881, 6491, 6117, 5759, 5417, 5088, 4774, 4473, 4185, 3909, 3646, 3395, 3155, 2927, 2710,
             2504, 2308, 2125, 1951, 1788, 1636, 1494, 1363, 1241, 1127, 1023, 926, 836, 753, 677, 606, 541, 481, 427, 375, 329, 287, 247, 
             211, 177, 146, 117, 89, 63, 38, 12, 0)
diffHeights <- abs(diff(heights))
heights2 <- c(11784.06, 9163.96, 7185.44, 5574.44, 4206.43, 3012.18, 1948.99, 1457.30, 988.50, 761.97, 540.34, 110.88, 0)
diffHeights2 <- abs(diff(heights2))

Pressures <- c(1000,950,925,900,850,800,700,600,500,400,300,200)
variableNumbersT <- c(280:291,66,270)
variableNumbersRH <- c(292:303,271)
variableNumbersCL <- c(67:131,265,269,268,267,272,279,304)
variableNumbersRAD <- c(277,275,276,274,273)

## functions to extract the date and leadtime
get_date <- function(afile){
  gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(afile)))
}
get_lt <- function(afile){
  gsub("00_GB", "", gsub("HA38_N25_SOLFC_[0-9]{12}_0", "", basename(afile)))
}

# write a function that imports the raster:
import_grib <- function(afile, param, adate, alt){
  tmpi      <- Gopen(afile)
  
  tmpg      <- Gdec(afile, param)
  tmpr      <- flip(raster(t(tmpg), 
                            xmn=attributes(tmpg)$domain$SW[1], 
                            xmx=attributes(tmpg)$domain$NE[1],
                            ymn=attributes(tmpg)$domain$SW[2], 
                            ymx=attributes(tmpg)$domain$NE[2]), 
                     direction="y")
  names(tmpr) <- paste0(tmpi$shortName[which(tmpi$position == param)], "_", adate, "_", alt)
  return(tmpr)
}

# import the rasters and binds and saves them
savingData <- function(dates){
  dataraster <- do.call(rbind, apply(files, 1, function(d) {
    for (j in 1:12){
      tmpcroprast1 <- crop(import_grib(afile = d, param = variableNumbersT[j], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
      tmpcroprast2 <- crop(import_grib(afile = d, param = variableNumbersRH[j], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
      tempt <- 1 - 373.15/values(tmpcroprast1)
      SVP <- 1013.25*exp(13.3185*tempt - 1.976*tempt^2 - 0.6445*tempt^3 - 0.1299*tempt^4)
      WaterVapor <- 10^6*values(tmpcroprast2) * SVP / Pressures[j]
      WaterVaporCorr[,j] <- WaterVapor * 10^(-6) * 18.02 * (Pressures[j]/10) / (8.3144 * values(tmpcroprast1))
      
      Temperatures[,j] <- values(tmpcroprast1)
      Humidities[,j] <- values(tmpcroprast2)
    }
    
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersT[13], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    T0meter <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersT[14], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    T2meter <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRH[13], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    RH2meter <- values(tmpcroprast)
    
    LowWaterVapor <- WaterVaporCorr[,1] * diffHeights2[12] + WaterVaporCorr[,2] * diffHeights2[11] + WaterVaporCorr[,3] * diffHeights2[10] + 
      WaterVaporCorr[,4] * diffHeights2[9] + WaterVaporCorr[,5] * diffHeights2[8] + WaterVaporCorr[,6] * diffHeights2[7]
    MediumWaterVapor <- WaterVaporCorr[,7] * diffHeights2[6] + WaterVaporCorr[,8] * diffHeights2[5] + WaterVaporCorr[,9] * diffHeights2[4]
    HighWaterVapor <- WaterVaporCorr[,10] * diffHeights2[3] + WaterVaporCorr[,11] * diffHeights2[2] + WaterVaporCorr[,12] * diffHeights2[1] 
    
    TotalWaterVapor <- LowWaterVapor + MediumWaterVapor + HighWaterVapor
    
    for (k in 1:65){
      tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[k], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
      CloudWaters[,k] <- values(tmpcroprast)*diffHeights[k]
    }
    
    HighCloudWater <- CloudWaters[,1]
    for (m in 2:22){
      HighCloudWater <- HighCloudWater + CloudWaters[,m]
    }
    
    MediumCloudWater <- CloudWaters[,23]
    for (m in 24:37){
      MediumCloudWater <- MediumCloudWater + CloudWaters[,m]
    }
    
    LowCloudWater <- CloudWaters[,38]
    for (m in 39:65){
      LowCloudWater <- LowCloudWater + CloudWaters[,m]
    }
    
    TotalCloudWater <- HighCloudWater + MediumCloudWater + LowCloudWater
    
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[66], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    TotalCloudCover <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[67], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    LowCloudCover <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[68], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    MediumCloudCover <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[69], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    HighCloudCover <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[70], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    RainAcc <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[71], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    CloudBases <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersCL[72], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    CloudTops <- values(tmpcroprast)
    
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRAD[1], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    GlobalRadiation <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRAD[2], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    DirectRadiationSURF <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRAD[3], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    DirectRadiationTOA <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRAD[4], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    NCSSURF <- values(tmpcroprast)
    tmpcroprast <- crop(import_grib(afile = d, param = variableNumbersRAD[5], adate = get_date(d), alt = get_lt(d)), extent(3.4, 7.1, 50.8, 53.5))
    NCSTOA <- values(tmpcroprast)
    
    xcoordinates <- coordinates(tmpcroprast)[,1]
    ycoordinates <- coordinates(tmpcroprast)[,2]
    
    tmpdf <- data.frame(lt = get_lt(d), xcoor = xcoordinates, ycoor = ycoordinates, T_0m = T0meter, T_2m = T2meter, 
                        T_1000 = Temperatures[,1], T_950 = Temperatures[,2], T_925 = Temperatures[,3], T_900 = Temperatures[,4], 
                        T_850 = Temperatures[,5], T_800 = Temperatures[,6], T_700 = Temperatures[,7], T_600 = Temperatures[,8], 
                        T_500 = Temperatures[,9], T_400 = Temperatures[,10], T_300 = Temperatures[,11], T_200 = Temperatures[,12],
                        lt = get_lt(d), xcoor = xcoordinates, ycoor = ycoordinates, RH_2m = RH2meter, 
                        RH_1000 = Humidities[,1], RH_950 = Humidities[,2], RH_925 = Humidities[,3], RH_900 = Humidities[,4], 
                        RH_850 = Humidities[,5], RH_800 = Humidities[,6], RH_700 = Humidities[,7], RH_600 = Humidities[,8], 
                        RH_500 = Humidities[,9], RH_400 = Humidities[,10], RH_300 = Humidities[,11], RH_200 = Humidities[,12],
                        lt = get_lt(d), xcoor = xcoordinates, ycoor = ycoordinates, 
                        CC_Total = TotalCloudCover, CC_Low = LowCloudCover, CC_Medium = MediumCloudCover, CC_High = HighCloudCover,
                        CW_Total = TotalCloudWater, CW_Low = LowCloudWater, CW_Medium = MediumCloudWater, CW_High = HighCloudWater, 
                        PW_Total = TotalWaterVapor, PW_Low = LowWaterVapor, PW_Medium = MediumWaterVapor, PW_High = HighWaterVapor,
                        Rain = RainAcc, Cloud_base = CloudBases, Cloud_top = CloudTops, 
                        lt = get_lt(d), xcoor = xcoordinates, ycoor = ycoordinates, Global = GlobalRadiation,
                        Direct_SURF = DirectRadiationSURF, Direct_TOA = DirectRadiationTOA, NCS_SURF = NCSSURF, NCS_TOA = NCSTOA)
    return(tmpdf)
  }))
  
  tmpfile1 <- dataraster[1:17]
  tmpfile2 <- dataraster[18:33]
  names(tmpfile2)[names(tmpfile2)=="lt.1"] <- "lt"
  names(tmpfile2)[names(tmpfile2)=="xcoor.1"] <- "xcoor"
  names(tmpfile2)[names(tmpfile2)=="ycoor.1"] <- "ycoor"
  tmpfile3 <- dataraster[34:51]
  names(tmpfile3)[names(tmpfile3)=="lt.2"] <- "lt"
  names(tmpfile3)[names(tmpfile3)=="xcoor.2"] <- "xcoor"
  names(tmpfile3)[names(tmpfile3)=="ycoor.2"] <- "ycoor"
  tmpfile4 <- dataraster[52:59]
  names(tmpfile4)[names(tmpfile4)=="lt.3"] <- "lt"
  names(tmpfile4)[names(tmpfile4)=="xcoor.3"] <- "xcoor"
  names(tmpfile4)[names(tmpfile4)=="ycoor.3"] <- "ycoor"
  
  saveRDS(tmpfile1, file=paste0("/nobackup/users/bakker/Data/temperaturevariables/",init_dates[dates],".rds"))
  saveRDS(tmpfile2, file=paste0("/nobackup/users/bakker/Data/humidityvariables/",init_dates[dates],".rds"))
  saveRDS(tmpfile3, file=paste0("/nobackup/users/bakker/Data/cloudvariables/",init_dates[dates],".rds"))
  saveRDS(tmpfile4, file=paste0("/nobackup/users/bakker/Data/radiationvariables/",init_dates[dates],".rds"))
}

yearNumber <- 2016
fnames    <- data.frame(files = list.files(path = paste0(indir, yearNumber), full.names = TRUE), stringsAsFactors = FALSE)
init_dates  <- unique(gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(fnames$files))))
init_dates <- init_dates[!init_dates %in% "20171211"]
init_dates <- init_dates[!init_dates %in% "20180117"]
WaterVaporCorr <- matrix(0,11817,12)
Temperatures <- matrix(0,11817,12)
Humidities <- matrix(0,11817,12)
CloudWaters <- matrix(0,11817,65)
for (i in 1:128){
  files <- as.matrix(c(as.matrix(fnames[grepl(init_dates[i], fnames$files), ])[5:21],as.matrix(fnames[grepl(init_dates[i], fnames$files), ])[29:45]))
  savingData(i)
}