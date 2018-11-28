# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
library(maps)
library(maptools)
library(ncdf4)
library(tidyverse)
library(rasterVis)
library(ggplot2)
library(caret)
library(leaps)
library(grid)
library(gridExtra)
library(glmnet)
library(quantreg)
library(rqPen)
library(gamlss)
library(gamlss.tr)
library(randomForest)
library(quantregForest)
library(grf)
library(gbm)
library(neuralnet)
library(qrnn)
library(e1071)
library(qrsvm)
library(rpart)
library(gptk)
library(class)
library(pryr)
library(devtools)
library(Rcpp)
sourceCpp("/usr/people/bakker/kilianbakker/R/crps_ensemble.cpp")
source("/usr/people/bakker/kilianbakker/R/functions.R")

# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                    52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                    6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                    344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                    "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                    "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                    "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
varSize <- varSizes[vars]
AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))

hyperparamText <- "musteps"
samps <- c(1/6,2/6,3/6,4/6,5/6)
mtries <- c(1/6,2/6,3/6,4/6,5/6)
nodesizes <- c(5,10,20,50)
ntrees <- c(100,500,1000,2000)
shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
depths <- c(1,5,10,20,40)
hiddens1 <- c(1,2,3)
hiddens2 <- c(1,2,3)
iters <- c(2,10,100)
musteps <- c(1,3,5,10)
sigmasteps <- c(0,1,3,5)
hyperparameters <- musteps
for (hyper in 1:length(hyperparameters)){

stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
leadTimes <- c(5:19, 29:43)
sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                 c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))

NumberofCV <- 3
setting <- 14
stationPlaces <- c(1:30)
time <- 8
season <- 1
CVmonth <- 1
repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
variableIndices <- c(1:length(variableNames))
QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))

for (setting in c(10,12,17)){
  ContMethod <- settings[[setting]][1]
  ProbMethod <- settings[[setting]][2] 
  RegMethod <- settings[[setting]][3]
  DistMethod <- settings[[setting]][4]

if (ProbMethod == 1){
  Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
} else if (ProbMethod == 2){
  Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
}

for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
  leadTime <- leadTimes[time]
  
  temp_init_dates <- init_dates
  if (leadTime > 24){
    for (d in 1:length(temp_init_dates)){
      tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
      temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
    }
  }
  
  temp_init_dates2 <- init_dates
  #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
  DateIndices <- which(repeated_init_dates %in% temp_init_dates2)

  tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
  obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
  tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
  CSR <- tempCSData[order(tempCSData$Date),][[4]]
  tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
  Zenith <- tempZAData[order(tempZAData$Date),][[4]]
  obs_TOA <- 1367*cos(Zenith*pi/180)

  tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
  variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
  for (date in 1:length(init_dates)){
    variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
  }
  
  for (var in 1:5){
    variableData[,var] <- variableData[,var]/CSR
  }

  if (ContMethod == 1){
    CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
    Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
    obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
    CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
  } else if (ContMethod == 2){
    CSI_obs = obs_SURF/CSR
    Ratio_obs = obs_SURF/obs_TOA
    CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
  }
  RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
  All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
  colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")

  for (season in 1:length(seasons)){
  DataMonths <- seasons[[season]]
  tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
  tempIndices <- sample(tempIndices)
  if (ProbMethod == 1){
    tempPredictionset <- array(NA,c(length(tempIndices)))
  } else if (ProbMethod == 2){
    tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
  }
  if (length(tempIndices) >= 50){
  for (CVmonth in 1:NumberofCV){
    testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
    TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
    TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
    #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
    #if (length(zeroVarColumns) > 0){
    #  TrainingData <- TrainingData[,-zeroVarColumns]
    #  TestingData <- TestingData[,-zeroVarColumns]
    #}
    
    Predictors <- TrainingData[,-length(colnames(TrainingData))]
    Predictand <-  TrainingData[,length(colnames(TrainingData))]
    
    #normal linear regression
    #variableNumber <- 1
    #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
    #FitPlot(plotData, 300, -0.5, 0.1)
    
    form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
    formula1 <- paste0(predictandNames[predictand], "~.")
    formula2 <-  paste0(predictandNames[predictand], "~1")
    formula3 <- paste0(predictandNames[predictand], form)

    #stepwise linear regression
   if (RegMethod == 1){
         fitcontrol <- gamlss.control(c.crit = 0.01, n.cyc = 50, mu.step = 1, sigma.step = 1, nu.step = 1, 
                                     tau.step = 1, gd.tol = Inf, iter = 0, trace = F, autostep = TRUE, save = T)
         fitcontrol2 <- glim.control(cc = 0.01, cyc = 50, glm.trace = F, bf.cyc = 50, bf.tol = 0.01, bf.trace = F)
        if (DistMethod %in% c(-1,0,1)){
          fit0 <- gamlss(formula = CSI_obs~1, sigma.formula = CSI_obs~1, data = TrainingData, family = NO(mu.link = "identity", sigma.link = "identity"), 
                         control = fitcontrol,i.control = fitcontrol2, trace = F)
          if (DistMethod == -1){
            fit1 <- fit0
          } else if (DistMethod == 0){
            fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = paste0("~Global")), direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          } else if (DistMethod == 1){
            fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = musteps[hyper], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
            fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmasteps[2], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          }
          Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
          Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.001, c(length(testIndices))))
          for (i in 1:length(testIndices)){
            tempPredictionset[testIndices[i],] <- unname(qNO(quants, mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
          }
        } else if (DistMethod %in% c(2,3)){
          if (DistMethod == 2){
            gen.trun(par=c(min(RadiationData[predictand][[1]][tempIndices], na.rm = T),
                           max(RadiationData[predictand][[1]][tempIndices], na.rm = T)),"NO", type="both")
          } else if (DistMethod == 3){
            gen.trun(par=c(0),"NO", type="left")
          }
          fit0 <- gamlss(formula = CSI_obs~1, sigma.formula = CSI_obs~1, data = TrainingData, family = NOtr(mu.link = "identity", sigma.link = "identity"), control = fitcontrol,i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = musteps[hyper], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmasteps[2], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
          Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.001, c(length(testIndices))))
          for (i in 1:length(testIndices)){
            tempPredictionset[testIndices[i],] <- unname(qNOtr(quants, mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
          }
        }
      } 
      tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
    if (RegMethod > 0){
      saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
      rm(fit1, TrainingData, TestingData)
    }
    } 
  }
  for (date in 1:length(init_dates)){
    tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
    if (ProbMethod == 1){
      Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
    } else if (ProbMethod == 2){
      Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
    }
  }
}
rm(variableData, All_Data, RadiationData)
}
saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText, hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
rm(Predictionset)
gc()
}
}
}







# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "sigmasteps"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- sigmasteps
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 14
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(10,12)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 1){
                  fitcontrol <- gamlss.control(c.crit = 0.01, n.cyc = 50, mu.step = 1, sigma.step = 1, nu.step = 1, 
                                               tau.step = 1, gd.tol = Inf, iter = 0, trace = F, autostep = TRUE, save = T)
                  fitcontrol2 <- glim.control(cc = 0.01, cyc = 50, glm.trace = F, bf.cyc = 50, bf.tol = 0.01, bf.trace = F)
                  if (DistMethod %in% c(-1,0,1)){
                    fit0 <- gamlss(formula = CSI_obs~1, sigma.formula = CSI_obs~1, data = TrainingData, family = NO(mu.link = "identity", sigma.link = "identity"), 
                                   control = fitcontrol,i.control = fitcontrol2, trace = F)
                    if (DistMethod == -1){
                      fit1 <- fit0
                    } else if (DistMethod == 0){
                      fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = paste0("~Global")), direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
                    } else if (DistMethod == 1){
                      fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = musteps[2], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
                      fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmasteps[hyper], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
                    }
                    Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
                    Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.001, c(length(testIndices))))
                    for (i in 1:length(testIndices)){
                      tempPredictionset[testIndices[i],] <- unname(qNO(quants, mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
                    }
                  } else if (DistMethod %in% c(2,3)){
                    if (DistMethod == 2){
                      gen.trun(par=c(min(RadiationData[predictand][[1]][tempIndices], na.rm = T),
                                     max(RadiationData[predictand][[1]][tempIndices], na.rm = T)),"NO", type="both")
                    } else if (DistMethod == 3){
                      gen.trun(par=c(0),"NO", type="left")
                    }
                    fit0 <- gamlss(formula = CSI_obs~1, sigma.formula = CSI_obs~1, data = TrainingData, family = NOtr(mu.link = "identity", sigma.link = "identity"), control = fitcontrol,i.control = fitcontrol2, trace = F)
                    fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = musteps[2], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
                    fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmasteps[hyper], direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
                    Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
                    Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.001, c(length(testIndices))))
                    for (i in 1:length(testIndices)){
                      tempPredictionset[testIndices[i],] <- unname(qNOtr(quants, mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
                    }
                  }
                }
              tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
library(maps)
library(maptools)
library(ncdf4)
library(tidyverse)
library(rasterVis)
library(ggplot2)
library(caret)
library(leaps)
library(grid)
library(gridExtra)
library(glmnet)
library(quantreg)
library(rqPen)
library(gamlss)
library(gamlss.tr)
library(randomForest)
library(quantregForest)
library(grf)
library(gbm)
library(neuralnet)
library(qrnn)
library(e1071)
library(qrsvm)
library(rpart)
library(gptk)
library(class)
library(pryr)
library(devtools)
library(Rcpp)
sourceCpp("/usr/people/bakker/kilianbakker/R/crps_ensemble.cpp")
source("/usr/people/bakker/kilianbakker/R/functions.R")

# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "ntrees"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- ntrees
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(14,15,18)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 4){
                  fit1 <- quantregForest(Predictors, Predictand, sampsize=ceiling(nrow(Predictors)*samps[3]), mtry = ceiling(ncol(Predictors)*mtries[2]), nodesize = nodesizes[1], ntree = ntrees[hyper])
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, what = quants)
                  #QRFimp[time,season,CVmonth,] <- importance(fit1)
                } else if (RegMethod == 5){
                  for (q in 1:length(quants)){
                    fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = ntrees[hyper], shrinkage=shrinkages[3],
                                interaction.depth=depths[1], bag.fraction = samps[3], train.fraction = 1, n.minobsinnode = nodesizes[1], verbose=FALSE)
                    tempPredictionset[testIndices,q] <- predict(fit1, newdata = TestingData, n.trees = ntrees[hyper])
                    #GBMimp[time,season,CVmonth,q,] <- varImp(fit1,gbmTrees)[[1]]
                  }
                } else if (RegMethod == 9){
                  fit1 <- quantile_forest(Predictors, Predictand, quantiles = quants, regression.splitting = FALSE, sample.fraction = samps[3], 
                                          mtry = ceiling(ncol(Predictors)*mtries[2]), num.trees = ntrees[hyper], min.node.size = nodesizes[1], honesty = F, honesty.fraction = NULL)
                  #split_frequencies(fit1, max.depth = 5)
                  #GRFimp[time,season,CVmonth,] <- variable_importance(fit1)[,1]
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, quantiles = quants)
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "mtries"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- mtries
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(14,18)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 4){
                  fit1 <- quantregForest(Predictors, Predictand, sampsize=ceiling(nrow(Predictors)*samps[3]), mtry = ceiling(ncol(Predictors)*mtries[hyper]), nodesize = nodesizes[1], ntree = ntrees[2])
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, what = quants)
                  #QRFimp[time,season,CVmonth,] <- importance(fit1)
                } else if (RegMethod == 5){
                  for (q in 1:length(quants)){
                    fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = ntrees[2], shrinkage=shrinkages[3],
                                interaction.depth=depths[1], bag.fraction = samps[3], train.fraction = 1, n.minobsinnode = nodesizes[1], verbose=FALSE)
                    tempPredictionset[testIndices,q] <- predict(fit1, newdata = TestingData, n.trees = ntrees[2])
                    #GBMimp[time,season,CVmonth,q,] <- varImp(fit1,gbmTrees)[[1]]
                  }
                } else if (RegMethod == 9){
                  fit1 <- quantile_forest(Predictors, Predictand, quantiles = quants, regression.splitting = FALSE, sample.fraction = samps[3], 
                                          mtry = ceiling(ncol(Predictors)*mtries[hyper]), num.trees = ntrees[2], min.node.size = nodesizes[1], honesty = F, honesty.fraction = NULL)
                  #split_frequencies(fit1, max.depth = 5)
                  #GRFimp[time,season,CVmonth,] <- variable_importance(fit1)[,1]
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, quantiles = quants)
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "nodesizes"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- nodesizes
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(14,15,18)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 4){
                  fit1 <- quantregForest(Predictors, Predictand, sampsize=ceiling(nrow(Predictors)*samps[3]), mtry = ceiling(ncol(Predictors)*mtries[2]), nodesize = nodesizes[hyper], ntree = ntrees[2])
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, what = quants)
                  #QRFimp[time,season,CVmonth,] <- importance(fit1)
                } else if (RegMethod == 5){
                  for (q in 1:length(quants)){
                    fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = ntrees[2], shrinkage=shrinkages[3],
                                interaction.depth=depths[1], bag.fraction = samps[3], train.fraction = 1, n.minobsinnode = nodesizes[hyper], verbose=FALSE)
                    tempPredictionset[testIndices,q] <- predict(fit1, newdata = TestingData, n.trees = ntrees[2])
                    #GBMimp[time,season,CVmonth,q,] <- varImp(fit1,gbmTrees)[[1]]
                  }
                } else if (RegMethod == 9){
                  fit1 <- quantile_forest(Predictors, Predictand, quantiles = quants, regression.splitting = FALSE, sample.fraction = samps[3], 
                                          mtry = ceiling(ncol(Predictors)*mtries[2]), num.trees = ntrees[2], min.node.size = nodesizes[hyper], honesty = F, honesty.fraction = NULL)
                  #split_frequencies(fit1, max.depth = 5)
                  #GRFimp[time,season,CVmonth,] <- variable_importance(fit1)[,1]
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, quantiles = quants)
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}





# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "samps"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- samps
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(14,15,18)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 4){
                  fit1 <- quantregForest(Predictors, Predictand, sampsize=ceiling(nrow(Predictors)*samps[hyper]), mtry = ceiling(ncol(Predictors)*mtries[2]), nodesize = nodesizes[1], ntree = ntrees[2])
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, what = quants)
                  #QRFimp[time,season,CVmonth,] <- importance(fit1)
                } else if (RegMethod == 9){
                  fit1 <- quantile_forest(Predictors, Predictand, quantiles = quants, regression.splitting = FALSE, sample.fraction = samps[hyper], 
                                          mtry = ceiling(ncol(Predictors)*mtries[2]), num.trees = ntrees[2], min.node.size = nodesizes[1], honesty = F, honesty.fraction = NULL)
                  #split_frequencies(fit1, max.depth = 5)
                  #GRFimp[time,season,CVmonth,] <- variable_importance(fit1)[,1]
                  tempPredictionset[testIndices,] <- predict(fit1, newdata = TestingData, quantiles = quants)
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
library(maps)
library(maptools)
library(ncdf4)
library(tidyverse)
library(rasterVis)
library(ggplot2)
library(caret)
library(leaps)
library(grid)
library(gridExtra)
library(glmnet)
library(quantreg)
library(rqPen)
library(gamlss)
library(gamlss.tr)
library(randomForest)
library(quantregForest)
library(grf)
library(gbm)
library(neuralnet)
library(qrnn)
library(e1071)
library(qrsvm)
library(rpart)
library(gptk)
library(class)
library(pryr)
library(devtools)
library(Rcpp)
sourceCpp("/usr/people/bakker/kilianbakker/R/crps_ensemble.cpp")
source("/usr/people/bakker/kilianbakker/R/functions.R")

# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "shrinkages"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- shrinkages
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(15)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 5){
                  for (q in 1:length(quants)){
                    fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = ntrees[2], shrinkage=shrinkages[hyper],
                                interaction.depth=depths[1], bag.fraction = samps[3], train.fraction = 1, n.minobsinnode = nodesizes[1], verbose=FALSE)
                    tempPredictionset[testIndices,q] <- predict(fit1, newdata = TestingData, n.trees = ntrees[2])
                    #GBMimp[time,season,CVmonth,q,] <- varImp(fit1,gbmTrees)[[1]]
                  }
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}



# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "depths"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- depths
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(15)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 5){
                  for (q in 1:length(quants)){
                    fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = ntrees[2], shrinkage=shrinkages[3],
                                interaction.depth=depths[hyper], bag.fraction = samps[3], train.fraction = 1, n.minobsinnode = nodesizes[1], verbose=FALSE)
                    tempPredictionset[testIndices,q] <- predict(fit1, newdata = TestingData, n.trees = ntrees[2])
                    #GBMimp[time,season,CVmonth,q,] <- varImp(fit1,gbmTrees)[[1]]
                  }
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "iters"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- iters
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(16)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 6){
                  fit1 <- mcqrnn.fit(as.matrix(Predictors), as.matrix(Predictand), n.hidden = hiddens1[1], n.hidden2 = hiddens2[1], tau = quants, iter.max=iters[hyper], n.trials=2,
                                     trace = F, Th = sigmoid, Th.prime = sigmoid.prime)
                  tempPredictionset[testIndices,] <- as.matrix(mcqrnn.predict(as.matrix(TestingData), fit1, tau = quants))
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "hiddens1"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- hiddens1
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(16)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
             if (RegMethod == 6){
                  fit1 <- mcqrnn.fit(as.matrix(Predictors), as.matrix(Predictand), n.hidden = hiddens1[hyper], n.hidden2 = hiddens2[1], tau = quants, iter.max=iters[2], n.trials=2,
                                     trace = F, Th = sigmoid, Th.prime = sigmoid.prime)
                  tempPredictionset[testIndices,] <- as.matrix(mcqrnn.predict(as.matrix(TestingData), fit1, tau = quants))
                }
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}




# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "hiddens2"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  hyperparameters <- hiddens2
  for (hyper in 1:length(hyperparameters)){
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(16)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 6){
                  fit1 <- mcqrnn.fit(as.matrix(Predictors), as.matrix(Predictand), n.hidden = hiddens1[1], n.hidden2 = hiddens2[hyper], tau = quants, iter.max=iters[2], n.trials=2,
                                     trace = F, Th = sigmoid, Th.prime = sigmoid.prime)
                  tempPredictionset[testIndices,] <- as.matrix(mcqrnn.predict(as.matrix(TestingData), fit1, tau = quants))
                } 
                tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}





# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data2/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)
stationNumber <- c(235, 240, 249, 251, 257, 260, 267, 269, 270, 273, 275, 277, 278, 279, 280, 283, 286, 290, 310, 319, 330, 
                   344, 348, 350, 356, 370, 375, 377, 380, 391)
stationName <- c("De Kooy", "Schiphol", "Berkhout", "Hoorn (Terschelling)", "Wijk aan zee", 
                 "De Bilt", "Stavoren", "Lelystad", "Leeuwarden", "Marknesse", "Deelen", "Lauwersoog", "Heino", 
                 "Hoogeveen", "Eelde", "Hupsel", "Nieuw Beerta", "Twenthe", "Vlissingen", "Westdorpe",
                 "Hoek van Holland", "Rotterdam", "Cabauw", "Gilze-Rijen", "Herwijnen", "Eindhoven", "Vonkel", "Ell", "Maastricht", "Arcen")
CoastDistance <- calculating_distances(stationxcoor, stationycoor, "Coast")
WaterDistance <- calculating_distances(stationxcoor, stationycoor, "Water")
InlandDistance <- calculating_distances(stationxcoor, stationycoor, "Inland")
stationData <- data.frame(Station = stationName, Number = stationNumber, xcoor = stationxcoor, ycoor = stationycoor, 
                          DistToCoast = CoastDistance, DistToWater = WaterDistance, DistToInland = InlandDistance)

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)

varSizes <- c("variables_1x1gridbox_1lt","variables_1x1gridbox_3lt","variables_1x1gridbox_3lt_RadNoAvg","variables_3x3gridbox_1lt", "variables_3x3gridbox_3lt", 
              "variables_3x3gridbox_3lt_RadNoAvg","variables_5x5gridbox_1lt", "variables_5x5gridbox_3lt", "variables_5x5gridbox_3lt_RadNoAvg",
              "variables_7x7gridbox_1lt", "variables_7x7gridbox_3lt", "variables_7x7gridbox_3lt_RadNoAvg","variables_9x9gridbox_1lt", "variables_9x9gridbox_3lt",
              "variables_9x9gridbox_3lt_RadNoAvg")[c(14)]
for (vars in 1:length(varSizes)){
  varSize <- varSizes[vars]
  AllvariableData <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/",varSize,".rds"))
  
  hyperparamText <- "penalties"
  samps <- c(1/6,2/6,3/6,4/6,5/6)
  mtries <- c(1/6,2/6,3/6,4/6,5/6)
  nodesizes <- c(5,10,20,50)
  ntrees <- c(100,500,1000,2000)
  shrinkages <- c(0.001,0.01,0.1,0.2,0.5,1)
  depths <- c(1,5,10,20,40)
  hiddens1 <- c(1,2,3)
  hiddens2 <- c(1,2,3)
  iters <- c(2,10,100)
  musteps <- c(1,3,5,10)
  sigmasteps <- c(0,1,3,5)
  penalties <- c(0,0.01,0.1,1)
  hyperparameters <- penalties
  for (hyper in 1:length(hyperparameters)){
    print("a")
    Sys.sleep(0.01)
    
    stepsize <- 1/50
    quants <- seq(stepsize,1-stepsize, stepsize)
    predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
    predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
    heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
    leadTimes <- c(5:19, 29:43)
    sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
    seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                    c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
    clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
    clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
    settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),
                     c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0),c(2,2,8,0),c(2,2,9,0),c(2,2,10,0))
    
    NumberofCV <- 3
    setting <- 4
    stationPlaces <- c(1:30)
    time <- 8
    season <- 1
    CVmonth <- 1
    repeated_init_dates <- rep(init_dates,each = length(stationPlaces))
    repeated_places <- rep(stationData[[2]][stationPlaces], length(init_dates))
    variableIndices <- c(1:length(variableNames))
    QRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    GBMimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(quants),length(variableIndices)))
    GRFimp <- array(NA, c(length(leadTimes),length(seasons),NumberofCV,length(variableIndices)))
    
    for (setting in c(16)){
      ContMethod <- settings[[setting]][1]
      ProbMethod <- settings[[setting]][2] 
      RegMethod <- settings[[setting]][3]
      DistMethod <- settings[[setting]][4]
      
      if (ProbMethod == 1){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces),length(leadTimes)))
      } else if (ProbMethod == 2){
        Predictionset <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants)))
      }
      
      for (time in c(4,8,12,19,23,27)){#1:length(leadTimes)){
        leadTime <- leadTimes[time]
        
        temp_init_dates <- init_dates
        if (leadTime > 24){
          for (d in 1:length(temp_init_dates)){
            tempDate <- paste0(substr(temp_init_dates[d],1,4), "-", substr(temp_init_dates[d],5,6), "-", substr(temp_init_dates[d],7,8))
            temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
          }
        }
        
        temp_init_dates2 <- init_dates
        #temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
        DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
        
        tempObsData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        obs_SURF <- tempObsData[order(tempObsData$Date),][[4]]
        tempCSData <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        CSR <- tempCSData[order(tempCSData$Date),][[4]]
        tempZAData <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
        Zenith <- tempZAData[order(tempZAData$Date),][[4]]
        obs_TOA <- 1367*cos(Zenith*pi/180)
        
        tempvariableData <- AllvariableData[, stationPlaces, time, variableIndices]
        variableData <- array(NA, c(length(init_dates)*length(stationPlaces),length(variableIndices)))
        for (date in 1:length(init_dates)){
          variableData[(length(stationPlaces)*(date - 1) + 1):(length(stationPlaces)*date),] <- tempvariableData[date,,]
        }
        
        for (var in 1:5){
          variableData[,var] <- variableData[,var]/CSR
        }
        
        if (ContMethod == 1){
          CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
          Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
          obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
          CSI_for = round(variableData[,which(variableNames[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
        } else if (ContMethod == 2){
          CSI_obs = obs_SURF/CSR
          Ratio_obs = obs_SURF/obs_TOA
          CSI_for = variableData[,which(variableNames[variableIndices] == "Global")]
        }
        RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = repeated_places, Date = repeated_init_dates)
        All_Data <- data.frame(variableData, RadiationData[predictand][[1]], repeated_places, repeated_init_dates)
        colnames(All_Data) <- c(variableNames[variableIndices], predictandNames[predictand], "Station", "Date")
        
        for (season in 1:length(seasons)){
          DataMonths <- seasons[[season]]
          tempIndices <- intersect(intersect(which(CSR > 20), which(floor(repeated_init_dates/100) %in% DataMonths)), DateIndices)
          tempIndices <- sample(tempIndices)
          if (ProbMethod == 1){
            tempPredictionset <- array(NA,c(length(tempIndices)))
          } else if (ProbMethod == 2){
            tempPredictionset <- array(NA,c(length(tempIndices),length(quants)))
          }
          if (length(tempIndices) >= 50){
            for (CVmonth in 1:NumberofCV){
              testIndices <- c(floor(length(tempIndices)*(CVmonth-1)/NumberofCV)+1):(floor(length(tempIndices)*CVmonth/NumberofCV))
              TrainingData <- All_Data[tempIndices,][-testIndices, 1:(length(colnames(All_Data))-2)]
              TestingData <- All_Data[tempIndices,][testIndices, 1:(length(colnames(All_Data))-3)]
              #zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
              #if (length(zeroVarColumns) > 0){
              #  TrainingData <- TrainingData[,-zeroVarColumns]
              #  TestingData <- TestingData[,-zeroVarColumns]
              #}
              
              Predictors <- TrainingData[,-length(colnames(TrainingData))]
              Predictand <-  TrainingData[,length(colnames(TrainingData))]
              
              #normal linear regression
              #variableNumber <- 1
              #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
              #FitPlot(plotData, 300, -0.5, 0.1)
              
              form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
              formula1 <- paste0(predictandNames[predictand], "~.")
              formula2 <-  paste0(predictandNames[predictand], "~1")
              formula3 <- paste0(predictandNames[predictand], form)
              
              #stepwise linear regression
              if (RegMethod == 6){
                fit1 <- mcqrnn.fit(as.matrix(Predictors), as.matrix(Predictand), n.hidden = hiddens1[1], n.hidden2 = hiddens2[1], tau = quants, iter.max=iters[2], n.trials=2,
                                   trace = F, Th = sigmoid, Th.prime = sigmoid.prime, penalty = penalties[hyper])
                tempPredictionset[testIndices,] <- as.matrix(mcqrnn.predict(as.matrix(TestingData), fit1, tau = quants))
              } 
              tempPredictionset[testIndices,] <- t(apply(tempPredictionset[testIndices,],1,sort))
              if (RegMethod > 0){
                saveRDS(fit1, file = paste0("/nobackup/users/bakker/Data2/fits/fit_", CVmonth, "_", season, "_", time, "_", hyperparamText, hyperparameters[hyper], "_", RegMethod,"_", DistMethod, ".rds"))
                rm(fit1, TrainingData, TestingData)
              }
            } 
          }
          for (date in 1:length(init_dates)){
            tempIndices2 <- which(tempIndices %in% c((length(stationPlaces)*(date - 1)+1):(length(stationPlaces)*date)))
            if (ProbMethod == 1){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time] <- tempPredictionset[tempIndices2]*CSR[tempIndices[tempIndices2]]
            } else if (ProbMethod == 2){
              Predictionset[date,tempIndices[tempIndices2] - length(stationPlaces)*(date - 1),time,] <- tempPredictionset[tempIndices2,]*array(CSR[tempIndices[tempIndices2]], c(length(tempIndices2), length(quants)))
            }
          }
        }
        rm(variableData, All_Data, RadiationData)
      }
      saveRDS(Predictionset, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",hyperparamText,hyperparameters[hyper], "_", RegMethod, "_", DistMethod,".rds"))
      rm(Predictionset)
      gc()
    }
  }
}