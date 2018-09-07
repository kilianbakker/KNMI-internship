# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
library(tidyverse)
library(rasterVis)
library(ggplot2)
library(grid)
library(gridExtra)
library(fitdistrplus)
library(gamlss)
library(gamlss.tr)
library(ncdf4)
library(Rcpp)
sourceCpp("/usr/people/bakker/kilianbakker/R/crps_ensemble.cpp")
source("/usr/people/bakker/kilianbakker/R/functions.R")

# The variables/constants
stationycoor   <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                    52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor   <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                    6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber  <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                    344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName    <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                    "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                    "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                    "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
stationData    <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor)
stationPlace   <- 23
leadTime       <- 12
NumberofDays   <- 30

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

#formulas to get the other variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

xcoordinate <- coordinatesLargeGrid(stationData,stationPlace, c(0), c(0))[[1]]
ycoordinate <- coordinatesLargeGrid(stationData,stationPlace, c(0), c(0))[[2]]
xcoordinate2 <- coordinatesSmallGrid(stationData,stationPlace)[[1]]
ycoordinate2 <- coordinatesSmallGrid(stationData,stationPlace)[[2]]

plotcolors <- c("blue", "red", "green", "orange", "steelblue1", "purple", "deeppink", "darkblue", "darkgreen", "black")

##TESTING AREA##

#################################
##FORECAST COMPARISON##
#################################

forecastDataW <- read.table(paste0("/net/bhw444/nobackup/users/veenvds/CABAUW/",yearNumber,"/",sprintf("%02d",monthNumber),"/",sprintf("%02d",dayNumber)
                                ,"/Cabauw_SW_N",dateNumber,"0000_+48"), col.names=c("xcoor","ycoor","direct","globaal","dekkingsgraad"), fill=TRUE)
forecastDataW <- filter(forecastDataW, xcoor == 134, ycoor == 130, direct != -99.0)
diurnalCycle1 <- c(0,0,0,diff(forecastDataW[[4]][3:25]))/3600

forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",dateNumber,".rds"))
forecastData <- filter(forecastData, xcoor == xcoordinate, ycoor = ycoordinate, as.numeric(lt)<=25 & as.numeric(lt) >= 2)
diurnalCycle2 <- c(0,0,0,diff(forecastData[[7]][2:24]))/3600

time = c(0,1,forecastData[[2]][1:23])
plotData = data.frame(time,diurnalCycle1,diurnalCycle2)

#ggplot(data = plotData) +
#  geom_line(mapping = aes(x = time, y = diurnalCycle1, group = 1, color = "blue")) +
#  geom_line(mapping = aes(x = time, y = diurnalCycle2, group = 1, color = "red")) +
#  scale_colour_manual(labels = c("forecastsW", "forecasts"), values = c("blue", "red")) +
#  xlab("time") + ylab("radiation")

ggplot(data = plotData) +
  geom_line(mapping = aes(x = time, y = diurnalCycle1 - diurnalCycle2))

#################################
##CHECKING FOR 128 RESOLUTION##
#################################

yearNumber <- 2016
fnames    <- data.frame(files = list.files(path = paste0(indir, yearNumber), full.names = TRUE), stringsAsFactors = FALSE)
init_dates  <- unique(gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(fnames$files))))

rastfiles <- as.matrix(fnames[grepl(init_dates[30], fnames$files), ])
newraster1 <- import_grib(rastfiles[12], param = 275, adate = "20180110", alt = "01")
newraster2 <- import_grib(rastfiles[12], param = 276, adate = "20180110", alt = "01")
newraster3 <- import_grib(rastfiles[13], param = 275, adate = "20180110", alt = "01")
newraster4 <- import_grib(rastfiles[13], param = 276, adate = "20180110", alt = "01")

#gplot(newraster) + 
#  geom_tile(aes(fill = value)) +
#  geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = NA, color = "black", size = 0.1, linetype = "dashed") +
#  scale_fill_distiller("", palette = "Spectral", direction=1, na.value = "white") +
#  coord_fixed(xlim = c(3.4, 7.1), ylim = c(50.8, 53.5))

tmpdf <- data.frame(xcoor = coordinates(newraster3)[,1], ycoor = coordinates(newraster3)[,2], 
                    SWbottomcum = values(newraster3), SWtopcum = values(newraster4),
                    SWbottom = (values(newraster3) - values(newraster1))/3600, SWtop = (values(newraster4) - values(newraster2))/3600)

comparison1      <- init_dates[10]
comparison2      <- extract(newraster3, SpatialPoints(cbind(4.927, 51.971)))

for (i in 2: as.integer(length(init_dates)/10)){
  rastfiles  <- as.matrix(fnames[grepl(init_dates[10*i], fnames$files), ])
  newraster3 <- import_grib(rastfiles[13], param = 275, adate = "20180110", alt = "01")
  comparison1 <- c(comparison1,init_dates[10*i])
  comparison2 <- c(comparison2,extract(newraster3, SpatialPoints(cbind(4.927, 51.971))))
}

yearNumber <- 2017
fnames    <- data.frame(files = list.files(path = paste0(indir, yearNumber), full.names = TRUE), stringsAsFactors = FALSE)
init_dates  <- unique(gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(fnames$files))))

for (i in 1: as.integer(length(init_dates)/10)){
  rastfiles  <- as.matrix(fnames[grepl(init_dates[10*i], fnames$files), ])
  newraster3 <- import_grib(rastfiles[13], param = 275, adate = "20180110", alt = "01")
  comparison1 <- c(comparison1,init_dates[10*i])
  comparison2 <- c(comparison2,extract(newraster3, SpatialPoints(cbind(4.927, 51.971))))
}

yearNumber <- 2018
fnames    <- data.frame(files = list.files(path = paste0(indir, yearNumber), full.names = TRUE), stringsAsFactors = FALSE)
init_dates  <- unique(gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(fnames$files))))

for (i in 1: as.integer(length(init_dates)/10)){
  rastfiles  <- as.matrix(fnames[grepl(init_dates[10*i], fnames$files), ])
  newraster3 <- import_grib(rastfiles[13], param = 275, adate = "20180110", alt = "01")
  comparison1 <- c(comparison1,init_dates[10*i])
  comparison2 <- c(comparison2,extract(newraster3, SpatialPoints(cbind(4.927, 51.971))))
}

comparison <- data.frame(date = comparison1, radiation = comparison2, 
                         radiation128 = comparison2/128, test1 = as.integer(comparison2/128) == comparison2/128,
                         radiation64 = comparison2/64, test2 = as.integer(comparison2/64) == comparison2/64)

#################################
##CHECKING FOR ZERO DIRECT RADIATION##
#################################

rastfiles <- as.matrix(fnames[grepl(init_dates[30], fnames$files), ])

newraster1 <- import_grib(rastfiles[1], param = 275, adate = "20180110", alt = "01")
newraster2 <- import_grib(rastfiles[1], param = 276, adate = "20180110", alt = "01")
testdiurnalCycle <- data.frame(date = init_dates[30], lt = i-1, SWbottom <- extract(newraster1, SpatialPoints(cbind(9.68012500, 55.70508))),
                               SWtop <- extract(newraster2, SpatialPoints(cbind(9.68012500, 55.70508))))

for (i in 2:49){
  newraster1 <- import_grib(rastfiles[i], param = 275, adate = "20180110", alt = "01")
  newraster2 <- import_grib(rastfiles[i], param = 276, adate = "20180110", alt = "01")
  testdiurnalCycle <- rbind(testdiurnalCycle, c(init_dates[30], i-1, extract(newraster1, SpatialPoints(cbind(9.68012500, 55.70508))),
                                                extract(newraster2, SpatialPoints(cbind(9.68012500, 55.70508)))))
}

testdiurnalCycle <- data.frame(testdiurnalCycle, as.numeric(testdiurnalCycle[[3]])/128, as.numeric(testdiurnalCycle[[4]])/128)
testdiurnalCycle <- data.frame(testdiurnalCycle, c(0,diff(as.numeric(testdiurnalCycle[[3]]))), c(0,diff(as.numeric(testdiurnalCycle[[4]]))))

ggplot(data = testdiurnalCycle) +
  geom_line(mapping = aes(x = as.numeric(lt), y = c.0..diff.as.numeric.testdiurnalCycle..3.....)) +
  geom_line(mapping = aes(x = as.numeric(lt), y = c.0..diff.as.numeric.testdiurnalCycle..4.....))

#################################
##TESTING FOR THE NET CLEAR SKY INDEX##
#################################

forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",dateNumber,".rds"))
forecastData <- filter(forecastData, xcoor == xcoordinate, ycoor = ycoordinate, as.numeric(lt) <= 25)
SWclearskyTop <- c(0,diff(forecastData[[9]])/3600)
SWclearskyBottom <- c(0,diff(forecastData[[8]])/3600)
SWTop <- c(0,diff(forecastData[[6]])/3600)
SWBottom <- c(0,diff(forecastData[[5]])/3600)
SWglobal <- c(0,diff(forecastData[[7]])/3600)

SWclearskyBottom <- (1/0.77)*SWclearskyBottom

time <- c(0:24)
plotData <- data.frame(time,SWclearskyTop, SWclearskyBottom, SWTop, SWBottom, SWglobal)

for (i in 1:length(plotData[[2]])){
  if (is.na(plotData[[2]][[i]]) == FALSE){
    if (plotData[[2]][[i]] < 1){
      plotData[[1]][[i]] <- NA
      plotData[[2]][[i]] <- NA
      plotData[[3]][[i]] <- NA
      plotData[[4]][[i]] <- NA
      plotData[[5]][[i]] <- NA
      plotData[[6]][[i]] <- NA
    }
  }
}

new_plotData <- data.frame(time = plotData$time,
                           data = c(plotData$SWclearskyTop, plotData$SWclearskyBottom, plotData$SWTop, plotData$SWBottom, plotData$SWglobal),
                           radType = rep(c('SWclearskyTop', 'SWclearskyBottom', 'SWTop', 'SWBottom', 'SWglobal'), each = nrow(plotData)))

ggplot(data = na.omit(new_plotData)) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  scale_colour_manual( values = c("blue", "red", "green", "yellow", "black")) +
  xlab("time") + ylab("radiation") + ggtitle("20180225")

ggplot(data = na.omit(subset(new_plotData, radType %in% c('SWclearskyBottom','SWglobal')))) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  xlab("time") + ylab("radiation") + ggtitle("20180225")

#################################
##TESTING FOR THE ZENITH ANGLE##
#################################

ZenithAngles1 <- c(1:24)
ZenithAngles2 <- c(1:24)

for (i in 1:24){
  ZenithAngles1[i] <- zenithAngle(i,20160621, 51.967, 4.917,0,1)
  ZenithAngles2[i] <- zenithAngle(i,20160621, 51.967, 4.917,0,2)
}

ZenithAngles3 <- acos(c(-0.219592, -0.152665, -0.0539943, 0.0696992, 0.209988, 0.357315, 0.501643, 0.633138, 0.742843, 0.823283, 0.868978, 0.876815, 0.846259, 
                    0.779392, 0.680769, 0.557109, 0.416837, 0.269507, 0.125159, -0.0063754, -0.116134, -0.19664, -0.242408, -0.250506))/pi*180

time = c(1:24)
plotData <- data.frame(time = plotData$time,
                           data = c(plotData$ZenithAngles1, plotData$ZenithAngles2, plotData$ZenithAngles3),
                           radType = rep(c('ZenithAngles1', 'ZenithAngles2', 'ZenithAngles3'), each = length(time)))

ggplot(data = na.omit(plotData)) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  scale_colour_manual(values = c("blue", "red", "green")) +
  xlab("Time") + ylab("Zenith angle")

#ggplot(data = plotData) +
#  geom_line(mapping = aes(x = time, y = ZenithAngles1 - ZenithAngles2))

#################################
##TESTING FOR THE CLEAR SKY OR LARGE ERROR DAYS
#################################
Times <- c(5:19)

for (k in 1:length(init_dates)){
  observations <- filter(observationData, Station == stationData[[2]][[stationPlace]], Date == init_dates[k])[[4]][Times]
  clearskyRadiations <- filter(clearskyData, Station == stationData[[2]][[stationPlace]], Date == init_dates[k])[[4]][Times]
  
  forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",init_dates[k],".rds"))
  forecasts <- diff(filter(forecastData, xcoor == xcoordinate, ycoor == ycoordinate)[[4]][c(1,Times-3)])/3600

plotData <- data.frame(time = Times, data = c(observations,forecasts), Type = rep(c("observations", "forecasts"), each = length(Times)))

plot <- ggplot(data = plotData) +
  geom_line(mapping = aes(x = plotData$time, y = data, color = Type)) +
  scale_colour_manual(values = c("blue", "red")) +
  xlab("Time") + ylab("Radiation") + ggtitle(paste0(init_dates[k]))
print(plot)
}

#################################
##TESTING FOR THE DIFFERENT SOURCES OF OBSERVATION DATA
#################################

dateNumbers <- 20170600 + c(1:30)
mean_observations <- array(0, c(17))
mean_observations_min_halfhours <- array(0, c(17))
mean_observations_min_hours <- array(0, c(17))
mean_observations_min_meanhour <- array(0, c(17))
mean_forecasts <- array(0, c(17))
for (k in 1:length(dateNumbers)){
  dateNumber <- dateNumbers[i]
observationData2 <- read.csv(paste0("/nobackup/users/bakker/Data/observations_min/",dateNumber,".csv"))
observations <- filter(observationData, Date == dateNumber, Station == stationData[[2]][[stationPlace]])[[4]][4:20]
observations_min_meanhour <- c(4:20)
for (i in 4:20){
  observations_min_meanhour[i-3] <- mean(pmax(observationData2[[7]][((i-1)*60+1):(i*60)], array(0,c(60))), na.rm = TRUE)
}
observations_min_halfhours <- pmax(observationData2[[7]][seq(31,1411,60)], array(0,c(24)))[4:20]
observations_min_hours <- pmax(observationData2[[7]][seq(60,1440,60)], array(0,c(24)))[4:20]

forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",dateNumber,".rds"))
forecasts <- diff(c(0,filter(forecastData, xcoor == xcoordinate, ycoor == ycoordinate)[[4]][1:17]))/3600

mean_observations <- mean_observations + observations
mean_observations_min_halfhours <- mean_observations_min_halfhours + observations_min_halfhours
mean_observations_min_hours <- mean_observations_min_hours + observations_min_hours
mean_observations_min_meanhour <- mean_observations_min_meanhour + observations_min_meanhour
mean_forecasts <- mean_forecasts + forecasts
}

mean_observations <- mean_observations/length(dateNumbers)
mean_observations_min_halfhours <- mean_observations_min_halfhours/length(dateNumbers)
mean_observations_min_hours <- mean_observations_min_hours/length(dateNumbers)
mean_observations_min_meanhour <- mean_observations_min_meanhour/length(dateNumbers)
mean_forecasts <- mean_forecasts/length(dateNumbers)

Times <- c(4:20)
plotData <- data.frame(time = Times,
                           data = c(mean_observations, mean_observations_min_halfhours, mean_observations_min_hours, 
                                    mean_observations_min_meanhour, mean_forecasts),
                           obsType = rep(c('mean_observations', 'mean_observations_min_halfhours', 'mean_observations_min_hours', 
                                           'mean_observations_min_meanhour', 'mean_forecasts'), each = length(Times)))
ggplot(data = plotData) +
  geom_line(mapping = aes(x = time, y = data, color = obsType)) +
  scale_colour_manual(values = c("blue", "green", "orange", "red", "black")) +
  xlab("Time") + ylab("Radiation") + ggtitle("20170601-20170630")

#################################
##TESTING FOR THE REGRESSION OF ALL THE VARIABLES
#################################
dateNumber <- 20171021
observations <- filter(observationData, Date == dateNumber, Station == stationData[[2]][[stationPlace]])[[4]][4:20]
clearsky <- filter(clearskyData, Date == dateNumber, Station == stationData[[2]][[stationPlace]])[[4]][4:20]

forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/", dateNumber, ".rds"))
forecastData <- filter(forecastData, xcoor == xcoordinate, ycoor == ycoordinate)
diurnalCycle <- diff(c(0,forecastData[[4]][1:17]))/3600

variableData <- readRDS(file = paste0("/nobackup/users/bakker/Data/cloudcovervariables/", dateNumber, ".rds"))
variableData <- filter(variableData, xcoor == xcoordinate, ycoor == ycoordinate)
variable <- variableData[[6]][1:17]

variableData <- readRDS(file = paste0("/nobackup/users/bakker/Data/cloudwatervariables/", dateNumber, ".rds"))
variableData <- filter(variableData, xcoor == xcoordinate, ycoor == ycoordinate)
variable2 <- variableData[[5]][1:17]

time <- c(4:20)

plotData <- data.frame(observations, diurnalCycle, clearsky, variable, variable2, time)
for (i in 1:(length(plotData[[1]]))){
  if (plotData[[1]][[i]] < 10){
    for (j in 1:(length(plotData)-3)){
      plotData[[j]][[i]] <- 0
    }
  }
}

plotData <- data.frame(time = plotData$time,
                           data = c(plotData$observations, plotData$diurnalCycle, plotData$clearsky, plotData$variable, plotData$variable2),
                           radType = rep(c('observations', 'diurnalCycle', 'clearsky','CC', 'CW'), each = nrow(plotData)))

p1 <- ggplot(data = subset(plotData, radType %in% c('observations','diurnalCycle','clearsky'))) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  scale_colour_manual( values = c("blue", "red", "green")) +
  xlab("time") + ylab("radiation")

p2 <- ggplot(data = subset(plotData, radType %in% c('CC','CW'))) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  scale_colour_manual( values = c("black","blue")) +
  xlab("time") + ylab("variable")

grid.arrange(p1, p2, ncol = 1)

#################################
##TESTING FOR THE REGRESSION OF ALL THE VARIABLES ON CLEAR SKY DAYS
#################################
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskies <- which(clearskyHours == 1, arr.ind = T)
cutoff <- max(which(clearskies[,2] <= length(init_dates)))
clearskies <- clearskies[1:cutoff,]

tempvariableNamesCorr <- c(variableNames[Includedvariables], "T_Low", "T_Medium", "T_High", "RH_Low", "RH_Medium", "RH_High")
variableNamesCorr <- c("Global", "PW_Total", "AOD_500", "Ang_exp", "Ozone")
stationPlace <- 23

AllvariableData <- readRDS(file  = "/usr/people/bakker/kilianbakker/Data/Allstations_variables")

observations <- c(1:length(clearskies[,1]))
forecasts <- c(1:length(clearskies[,1]))
CSR <- c(1:length(clearskies[,1]))
PW <- c(1:length(clearskies[,1]))
AOD <- c(1:length(clearskies[,1]))
ANG <- c(1:length(clearskies[,1]))
OZ <- c(1:length(clearskies[,1]))
for (i in 1:length(clearskies[,1])){
  observations[i] <- filter(observationData, Date == init_dates[clearskies[i,2]], Station == stationData[[2]][[stationPlace]])[[4]][clearskies[i,1]]
  CSR[i] <- ClearskyRadiation(clearskies[i,1] + 0.5, init_dates[clearskies[i,2]],stationData[[4]][stationPlace],stationData[[3]][stationPlace],0, 1, 1)

  variableData <- AllvariableData[length(init_dates)*(stationPlace - 1) + clearskies[i,2],  which(tempvariableNamesCorr %in% variableNamesCorr), clearskies[i,1]]
  forecasts[i] <- variableData[1]
  PW[i] <- variableData[2]
  AOD[i] <- variableData[3]
  ANG[i] <- variableData[4]
  OZ[i] <- variableData[5]
}
  
plotData <- data.frame(observations, forecasts = forecasts*CSR, CSR, PW, AOD, ANG, OZ)

ggplot(data = plotData) +
  geom_point(mapping = aes(x = observations, y = OZ))

#################################
##TESTING FOR THE DISTRIBUTION OF OBSERVATIONS
#################################
init_dates <-as.numeric(init_dates)
observations <- filter(observationData, Station == stationData[[2]][[stationPlace]], Time == leadTime, Date %in% init_dates)[[4]]
CSR <- filter(clearskyData, Station == stationData[[2]][[stationPlace]], Time == leadTime, Date %in% init_dates)[[4]]
Zenith <- filter(zenithangleData, Station == stationData[[2]][[stationPlace]], Time == leadTime, Date %in% init_dates)[[4]]
obs_TOA <- 1367*cos(Zenith*pi/180)

CSI_obs <- observations/CSR
Ratio_obs <- observations/obs_TOA

predictandData <- data.frame(init_dates, observations, CSI_obs, Ratio_obs)
colnames(predictandData) <- c("Dates", "observations", "CSI_obs", "Ratio_obs")
months1 <- c(201612,201701,201702,201712,201801,201802)
months2 <- c(201604,201605,201703,201704,201705,201803)
months3 <- c(201606,201607,201608,201706,201707,201708)
months4 <- c(201609,201610,201611,201709,201710,201711)

temppredictandData1 <- filter(predictandData, floor(Dates/100) %in% months1)
p1 <- ggplot(mapping = aes(x = CSI_obs), data = temppredictandData1) +
  geom_histogram(aes(y=..density..), binwidth = 0.05) +
  geom_density(kernel = "gaussian") + ggtitle("Winter")

temppredictandData2 <- filter(predictandData, floor(Dates/100) %in% months2)
p2 <- ggplot(mapping = aes(x = CSI_obs), data = temppredictandData2) +
  geom_histogram(aes(y=..density..), binwidth = 0.05) +
  geom_density(kernel = "gaussian") + ggtitle("Spring")

temppredictandData3 <- filter(predictandData, floor(Dates/100) %in% months3)
p3 <- ggplot(mapping = aes(x = CSI_obs), data = temppredictandData3) +
  geom_histogram(aes(y=..density..), binwidth = 0.05) +
  geom_density(kernel = "gaussian") + ggtitle("Summer")

temppredictandData4 <- filter(predictandData, floor(Dates/100) %in% months4)
p4 <- ggplot(mapping = aes(x = CSI_obs), data = temppredictandData4) +
  geom_histogram(aes(y=..density..), binwidth = 0.05) +
  geom_density(kernel = "gaussian") + ggtitle("Fall")

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

gen.trun(par=c(min(CSI_obs),max(CSI_obs)),"NO", type="both")
fitcontrol <- gamlss.control(c.crit = 0.01, n.cyc = 200, mu.step = 1, sigma.step = 0.1, nu.step = 1, tau.step = 1)
fitcontrol2 <- glim.control(cc = 0.01, cyc = 50, glm.trace = F, bf.cyc = 200, bf.tol = 0.01, bf.trace = F)

histos <- par(mfrow = c(2,2))
hist1 <- histDist(temppredictandData1[[3]], family = GB2, control = fitcontrol, i.control = fitcontrol2, nbins=50, density = TRUE, main = "GB2", xlab = "CSI_obs")
hist2 <- histDist(temppredictandData1[[3]], family = BCPE, control = fitcontrol, i.control = fitcontrol2, nbins=50, density = TRUE, main = "BCPE", xlab = "CSI_obs")
hist3 <- histDist(temppredictandData1[[3]], family = GG, control = fitcontrol, i.control = fitcontrol2, nbins=50, density = TRUE, main = "GG", xlab = "CSI_obs")
hist4 <- histDist(temppredictandData1[[3]], family = NOtr, control = fitcontrol, i.control = fitcontrol2, xlim = c(min(CSI_obs),max(CSI_obs)),nbins=50, density = TRUE, main = "NOtr", xlab = "CSI_obs")
par(histos)

m1<-fitDist(temppredictandData1[[3]], type="realplus", extra=c("NOtr","NO"))
m1$fits
m1$failed

plotdist(CSI_obs, distr = "norm", para=list(mean=mean(CSI_obs),sd=sd(CSI_obs)), histo = TRUE, demp = TRUE)

fit <- fitdist(CSI_obs, "lnorm")
summary(fit)
plot(fit)
boot <- bootdist(fit, niter = 1000)
summary(boot)
plot(boot)
descdist(CSI_obs, discrete=FALSE, boot=500)

plot(density(CSI_obs))

h<-hist(CSI_obs, breaks=50, col="black", xlab="CSI_obs")
xfit<-seq(min(CSI_obs),max(CSI_obs),length=40) 
yfit<-dNO(xfit,mu=mean(CSI_obs),sigma=sd(CSI_obs)) 
yfit2<-dNOtr(xfit,mu=mean(CSI_obs),sigma=sd(CSI_obs)) 
yfit <- yfit*diff(h$mids[1:2])*length(CSI_obs)
yfit2 <- yfit2*diff(h$mids[1:2])*length(CSI_obs)
lines(xfit, yfit, col="blue", lwd=2)
lines(xfit, yfit2, col="blue", lwd=2)

plot(function(y) dNOtr(y, mu=0.5 ,sigma=1), 0,1)
plot(function(y) pNOtr(y, mu=0.5 ,sigma=1), 0, 1)
plot(function(y) qNOtr(y, mu=0.5 ,sigma=1), 0.1, 0.9)

m1<-fitDist(CSI_obs, type="real0to1", extra=c("NOtr"))
m1$fits
m1$failed
tempData <- runif(1000, min = 0.01, max = 0.99)
Predictionset <- fitDistPred(CSI_obs, type="real0to1", extra=c("NOtr"), newdata=tempData)
Predictionset$fits
Predictionset$failed

t1 <- chooseDist(fitmethod4, type="realline", parallel="snow", ncpus=4)
# ordering according to BIC
getOrder(t1,3)

#################################
##TESTING FOR GAMLSS
#################################
X1 <- rnorm(1000,1,2)
X2 <- rnorm(1000,2,1)
X3 <- runif(1000,1,2)
X4 <- runif(1000,1,5)
X5 <- rnorm(1000,2.5,2.5)
mu <- 1.5 + X1 + 3*X3 + 4*X5
sigma <- 2*X4

X1 <- runif(1000,0.5,0.9)
X2 <- runif(1000,0.1,0.8)
X3 <- runif(1000,0.1,0.9)
X4 <- runif(1000,0.2,0.3)
X5 <- runif(1000,0.1,0.7)
mu <- 0.5 + 0.1*X1 + 0.3*X3 - 0.4*X5
sigma <- 0.7 + 2*X4 - 1*X1

predictand <- c(1:1000)
for (i in 1:1000){
  predictand[i] <- rnorm(1,mu[i],sigma[i])
  #predictand[i] <- runif(1,mu[i] - sigma[i], mu[i] + sigma[i])
}

fit0 <- gamlss(formula = predictand ~ 1, data = data.frame(X1,X2,X3,X4,X5,predictand), family = NO(mu.link = "identity", sigma.link = "identity"))
fitscope <- gamlss.scope(data.frame(X1,X2,X3,X4,X5,predictand), response = 6, smoother = "cs", form = T)
fit <- stepGAIC(fit0, scope= list(upper = as.formula(paste0("~", paste(c("X1","X2","X3","X4","X5"), collapse = "+")))), 
                       sigma.scope = list(upper = as.formula(paste0("~", paste(c("X1","X2","X3","X4","X5"), collapse = "+")))), 
                       direction = "both", trace = T)
fit <- stepGAICAll.B(fit0, what="mu", scope = list(upper = as.formula(paste0("~", paste(c("X1","X2","X3","X4","X5"), collapse = "+")))), direction = "both", trace = T)
fit <- stepGAICAll.B(fit, what="sigma", scope = list(upper = as.formula(paste0("~", paste(c("X1","X2","X3","X4","X5"), collapse = "+")))), direction = "both", trace = T)
summary(fit)

#################################
##TESTING FOR THE PREDICTION SETS
#################################

RegMethodTexts <- c("RAW", "GLM",  "Lasso", "RF", "GBM", "NN")
RegMethods <- c(0,1,3,4,5,6)
clearskydateNumbers <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128, 20170409,20170526, 20171015, 20180207,20180225)
tmprds <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
largeerrordateNumbers <- init_dates[as.numeric(which(apply(tmprds, 2, mean, na.rm = T) > 0.75))]
dateNumbers <- largeerrordateNumbers

clearskyRadiations <- array(0, c(15, length(dateNumbers)))
observations <- array(0, c(15,length(dateNumbers)))
forecasts <- array(0, c(15,length(dateNumbers), length(RegMethods)+1))
Times <- c(5:19)
for (j in 1:length(dateNumbers)){
  clearskyRadiations[,j] <- filter(clearskyData, Station == stationData[[2]][[stationPlace]], Date == dateNumbers[j])[[4]][Times]
  observations[,j] <- filter(observationData, Station == stationData[[2]][[stationPlace]], Date == dateNumbers[j])[[4]][Times]
  
  for (RegMethod in 1:length(RegMethods)){
    forecastData <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_2_1_", RegMethods[RegMethod], "_0.rds"))
    forecasts[,j,RegMethod] <- forecastData[which(init_dates == dateNumbers[j]),Times-4]
  }
  forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",dateNumbers[j],".rds"))
  tempforecasts <- filter(forecastData, xcoor == xcoordinate, ycoor == ycoordinate)[[4]]
  forecasts[,j,length(RegMethods) + 1] <- diff(tempforecasts[1:16])/3600
}

ColorRowMean <- 1
meanForecasts <- array(0, c(length(Times), length(RegMethods)+1))
for (RegMethod in 1:(length(RegMethods)+1)){
  meanForecasts[,RegMethod] <- apply(as.matrix(forecasts[,,RegMethod]),ColorRowMean,sum,na.rm = T)/length(dateNumbers)
}

RadiationData <- data.frame(apply(clearskyRadiations,ColorRowMean,mean), apply(observations, ColorRowMean, mean), meanForecasts)
plotData <- data.frame(time = Times, data = c(as.matrix(RadiationData[,2:9])),
                       radType = rep(c("clearskyRadiations", "observations", RegMethodTexts, "Raw2")[2:9], each = length(Times)))

ggplot(data = plotData) +
  geom_line(mapping = aes(x = time, y = data, color = radType)) +
  scale_colour_manual(values = plotcolors[1:length(colnames(RadiationData[,2:9]))]) +
  xlab("Time") + ylab("Radiation") + ggtitle(dateNumbers)

#################################
##TESTING FOR THE (SKILL) SCORES
#################################
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
LargeErrorHours <-readRDS("/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
LargeErrorHours <- LargeErrorHours[,which(as.numeric(colnames(LargeErrorHours)) %in% init_dates)]

typeCS <- 1
typeZenith <- 1
discretizeWidth <- 0.1
stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
ProbMethodTexts <- c("Deterministic (RMSE SS)", "Probabilistic (CRPSS)")
leadTimesDays <- list(c(5:19), c(29:43))
leadTimeIndices <- list(c(leadTimesDays[[1]]-4),c(leadTimesDays[[2]]-13))
leadTimes <- c(leadTimesDays[[1]],leadTimesDays[[2]])
leadTimesTexts <- c(leadTimesDays[[1]],"-",leadTimesDays[[2]])
settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GLM",  "Lasso", "RF", "GBM", "NN", "RAW", "NO", "NOtr[min,max)", "NOtr[0,Inf]", "Lasso", "QRF", "GBM", "QRNN")
settings <- settings[7:13]
MethodsTexts <- MethodsTexts[7:13]

stationPlace <- 23
Scores <- array(NA, c(length(leadTimes), length(settings),4))
Skillscores <- array(NA, c(length(leadTimes), length(settings),4))
for (setting in 1:length(settings)){
  ContMethod <- settings[[setting]][1]
  ProbMethod <- settings[[setting]][2]
  RegMethod <- settings[[setting]][3]
  DistMethod <- settings[[setting]][4]

  importingforecasts <- 1
for (time in 1:length(leadTimes)){
  leadTime <- leadTimes[time]
  if (leadTime <= 24){
    temp_init_dates <- init_dates
  } else {
    temp_init_dates <- c(1:length(init_dates))
    for (d in 1:length(init_dates)){
      tempDate <- paste0(substr(init_dates[d],1,4), "-", substr(init_dates[d],5,6), "-", substr(init_dates[d],7,8))
      temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
    }
  }

  observations <- filter(observationData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]]
  CSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]]
  
  #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128, 20170409,20170526, 20171015, 20180207,20180225)
  #temp_init_dates2 <- init_dates[which(clearskyHours[leadTime,] == 1)]
  #temp_init_dates2 <- init_dates[which(LargeErrorHours[leadTime,] == 1)]
  #temp_init_dates2 <- init_dates[floor(init_dates/100) %in% seasons[[season]]]
  temp_init_dates2 <- init_dates
  
  indices <- which(init_dates %in% temp_init_dates2)
  if (length(indices) > 0){

  tempobservations <- observations[indices]
  tempCSR <- CSR[indices]

  if (importingforecasts == 1){
    forecasts1 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
    forecasts2 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_year_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
    forecasts3 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_allstations_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
    forecasts4 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_allstations_year_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
    importingforecasts <- 0  
  }
  
  if (ProbMethod == 1){
    clim_scores_det <- calculating_sampleclimatologies(ContMethod, ProbMethod, tempobservations, tempCSR, temp_init_dates2, quants, discretizeWidth)
    Scores[time,setting,1] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts1[indices,stationPlace,time], DiscretizeWidth)
    Skillscores[time,setting,1] <- 1 - Scores[time, setting,1]/clim_scores_det
    Scores[time,setting,2] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts2[indices,stationPlace,time], DiscretizeWidth)
    Skillscores[time,setting,2] <- 1 - Scores[time, setting,2]/clim_scores_det
    Scores[time,setting,3] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts3[indices,stationPlace,time], DiscretizeWidth)
    Skillscores[time,setting,3] <- 1 - Scores[time, setting,3]/clim_scores_det
    Scores[time,setting,4] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts4[indices,stationPlace,time], DiscretizeWidth)
    Skillscores[time,setting,4] <- 1 - Scores[time, setting,4]/clim_scores_det
  } else if (ProbMethod == 2){
    clim_scores_prob <- calculating_sampleclimatologies(ContMethod, ProbMethod, tempobservations, tempCSR, temp_init_dates2, quants, discretizeWidth)
    Scores[time,setting,1] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts1[indices,stationPlace,time,], DiscretizeWidth)
    Skillscores[time,setting,1] <- 1 - Scores[time, setting,1]/clim_scores_prob
    Scores[time,setting,2] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts2[indices,stationPlace,time,], DiscretizeWidth)
    Skillscores[time,setting,2] <- 1 - Scores[time, setting,2]/clim_scores_prob
    Scores[time,setting,3] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts3[indices,stationPlace,time,], DiscretizeWidth)
    Skillscores[time,setting,3] <- 1 - Scores[time, setting,3]/clim_scores_prob
    Scores[time,setting,4] <- calculating_scores(ContMethod, ProbMethod, tempobservations, forecasts4[indices,stationPlace,time,], DiscretizeWidth)
    Skillscores[time,setting,4] <- 1 - Scores[time, setting,4]/clim_scores_prob
  }
  }
}
}

tempProbMethod <- 2
MethodsData <- list(list(c(),c()),list(c(),c()),list(c(),c()),list(c(),c()))
for (i in 1:4){
  tempMethodsTexts <- c()
  for (setting in 1:length(settings)){
    if (settings[[setting]][2] == tempProbMethod){
      MethodsData[[i]][[1]] <- c(MethodsData[[i]][[1]], Scores[leadTimeIndices[[1]],setting,i], NA, Scores[leadTimeIndices[[2]],setting,i])
      MethodsData[[i]][[2]] <- c(MethodsData[[i]][[2]], Skillscores[leadTimeIndices[[1]],setting,i], NA, Skillscores[leadTimeIndices[[2]],setting,i])
      tempMethodsTexts <- c(tempMethodsTexts, MethodsTexts[setting])
    }
  }
}

plotData <- data.frame(leadTime = leadTimesTexts, data = MethodsData[[2]],
                           Method = rep(tempMethodsTexts, each = length(leadTimesTexts)))

tempPlot <- ggplot(data = plotData) +
  geom_line(mapping = aes(x = rep(seq(leadTimesTexts),length(tempMethodsTexts)), y = data, color = Method)) +
  scale_colour_manual(values = plotcolors[1:(length(tempMethodsTexts))]) +
  scale_x_continuous(breaks=seq(leadTimesTexts), labels=leadTimesTexts) + 
  ylim(0,0.75) + xlab("leadTime") + ylab("Skill scores") + ggtitle(paste0(ProbMethodTexts[tempProbMethod]))
print(tempPlot)

