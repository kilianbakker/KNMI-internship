# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
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
source("/usr/people/bakker/kilianbakker/R/functions.R")

# The variables/constants
fnames         <- data.frame(files = list.files(path = "/nobackup/users/bakker/Data/radiationvariables", full.names = TRUE), stringsAsFactors = FALSE)
init_dates     <- as.numeric(unique(gsub(".rds", "", basename(fnames$files))))

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

variableNames <- c("T_0m", "T_2m", "T_1000", "T_950", "T_925", "T_900", "T_850", "T_800", "T_700", "T_600", "T_500", "T_400", "T_300", "T_200",
                   "RH_2m", "RH_1000", "RH_950", "RH_925", "RH_900", "RH_850", "RH_800", "RH_700", "RH_600", "RH_500", "RH_400", "RH_300", "RH_200",
                   "Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High",
                   "CW_Total", "CW_Low", "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "Cloud_base", "Cloud_top",
                   "AOD_500", "Ang_exp", "Ozone")

leadTimeSizes <- c(-1,0,1)
xDomainSizes <- c(-1,0,1)
yDomainSizes <- c(-1,0,1)
typeCS <- 1
typeZenith <- 1
threshold <- 0.8
stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
Includedvariables <- c(28:45,48:50)
DiscretizeWidth <- 0.1
predictand <- 2 #1 = obs_SURF, 2 = obs_SURF/CSR, 3 = obs_SURF/obs_TOA
predictandNames <- c("Radiation_obs", "CSI_obs", "Ratio_obs")
heights2 <- c(0, 2, 110.88, 540.34, 761.97, 988.50, 1457.30, 1948.99, 3012.18, 4206.43, 5574.44, 7185.44, 9163.96, 11784.06)
leadTimes <- c(5:19, 29:43)
variableNamesCorr <- c(variableNames[Includedvariables], "T_Low", "T_Medium", "T_High", "RH_Low", "RH_Medium", "RH_High")
averagingMethod <- 4
sigmaPredictors <- c(0,0,0,1,1,1,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2,2,1,1,1,0,0,0)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0))

NumberofCV <- 6
setting <- 6
stationPlaces <- c(23)
time <- 8
season <- 1
CVmonth <- 1
repeated_init_dates <- rep(init_dates, length(stationPlaces))

temp_init_dates2 <- repeated_init_dates
#temp_init_dates2 <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170526,20171015,20180207,20180225)
DateIndices <- which(repeated_init_dates %in% temp_init_dates2)
variableIndices <- c(1:length(variableNamesCorr))

clearskyData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
zenithangleData <- readRDS(file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
observationData[[4]] <- observationData[[4]]*(10000/3600)
AllvariableData <- readRDS(file  = "/usr/people/bakker/kilianbakker/Data/Allstations_3x3gridbox_3lt_variables.rds")

for (setting in c(11,14)){#1:length(settings)){
  ContMethod <- settings[[setting]][1]
  ProbMethod <- settings[[setting]][2] 
  RegMethod <- settings[[setting]][3]
  DistMethod <- settings[[setting]][4]

if (ProbMethod == 1){
  Predictionset <- array(NA, c(length(repeated_init_dates[DateIndices]), length(leadTimes)))
} else if (ProbMethod == 2){
  Predictionset <- array(NA, c(length(repeated_init_dates[DateIndices]), length(quants), length(leadTimes)))
}

for (time in 1:length(leadTimes)){
  leadTime <- leadTimes[time]
  
  if (leadTime <= 24){
    temp_init_dates <- init_dates
  } else {
    temp_init_dates <- c(init_dates[-1],20180317)
  }

  tempobservationData <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)
  obs_SURF <- tempobservationData[[4]][DateIndices]
  CSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]][DateIndices]
  Zenith <- filter(zenithangleData, Station %in% stationData[[2]][stationPlaces], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]][DateIndices]
  obs_TOA <- 1367*cos(Zenith*pi/180)
  cosZenith <- cos(Zenith*pi/180)
  
  variableData <- AllvariableData[DateIndices, variableIndices, time]
  
  if (ContMethod == 1){
    CSI_obs <- round(obs_SURF/CSR, digits = log(1/DiscretizeWidth)/log(10))
    Ratio_obs <- round(obs_SURF/obs_TOA, digits = log(1/DiscretizeWidth)/log(10))
    obs_SURF <- round(obs_SURF, digits = log(1/DiscretizeWidth)/log(10))
    CSI_for = round(variableData[,which(variableNamesCorr[variableIndices] == "Global")], digits = log(1/DiscretizeWidth)/log(10))
  } else if (ContMethod == 2){
    CSI_obs = obs_SURF/CSR
    Ratio_obs = obs_SURF/obs_TOA
    CSI_for = variableData[,which(variableNamesCorr[variableIndices] == "Global")]
  }
  RadiationData <- data.frame(obs_SURF, CSI_obs, Ratio_obs, CSI_for, Station = tempobservationData[[1]][DateIndices], Date = repeated_init_dates[DateIndices])
  All_Data <- data.frame(variableData, cosZenith, RadiationData[predictand][[1]], tempobservationData[[1]][DateIndices], repeated_init_dates[DateIndices])
  colnames(All_Data) <- c(variableNamesCorr[variableIndices], "cosZenith", predictandNames[predictand], "Station", "Date")

  for (season in 1:1){
  DataMonths <- c(seasons[[1]],seasons[[2]],seasons[[3]],seasons[[4]])
  init_dates_temp <- Reduce(intersect, list(repeated_init_dates[DateIndices][which(CSR > 20)], repeated_init_dates[floor(repeated_init_dates/100) %in% DataMonths]))
  #init_dates_temp <- Reduce(intersect, list(repeated_init_dates[which(CSR > 20)], repeated_init_dates[DateIndices]))
  init_dates_temp <- rep(init_dates_temp, length(stationPlaces))
  if (ProbMethod == 1){
    tempPredictionset <- array(NA,c(length(init_dates_temp)))
  } else if (ProbMethod == 2){
    tempPredictionset <- array(NA,c(length(init_dates_temp),length(quants)))
  }
  if (length(init_dates_temp) >= 50){
  for (CVmonth in 1:NumberofCV){
    tempIndices <- c(floor(length(init_dates_temp)*(CVmonth-1)/NumberofCV)+1):(floor(length(init_dates_temp)*CVmonth/NumberofCV))
    TrainingData <- filter(All_Data, Date %in% init_dates_temp)[-tempIndices, 1:(length(colnames(All_Data))-2)]
    TestingData <- filter(All_Data, Date %in% init_dates_temp)[tempIndices, 1:(length(colnames(All_Data))-3)]
    zeroVarColumns <- c(which(as.numeric(apply(TrainingData,2,var)) <= 0.0001))
    if (length(zeroVarColumns) > 0){
      TrainingData <- TrainingData[,-zeroVarColumns]
      TestingData <- TestingData[,-zeroVarColumns]
    }
    
    #normal linear regression
    #variableNumber <- 1
    #plotData <- TrainingData[,c(variableNumber, length(TrainingData[1,])-1)]
    #FitPlot(plotData, 300, -0.5, 0.1)
    
    form <- paste0("~", paste(colnames(TestingData), collapse = "+"))
    formula1 <- paste0(predictandNames[predictand], "~.")
    formula2 <-  paste0(predictandNames[predictand], "~1")
    formula3 <- paste0(predictandNames[predictand], form)
    
    #stepwise linear regression
    if (ProbMethod == 1){
      if (RegMethod == 0){
        tempPredictionset[tempIndices] <- filter(RadiationData, Date %in% init_dates_temp)[[4]][tempIndices]
      } else if (RegMethod == 1){
      fit0 <- glm(formula2, family = gaussian(link = "identity"), data = TrainingData)
      fit1 <- stepAIC(fit0, scope = list(upper = form), steps = 3, direction = "both", trace = FALSE)
      tempPredictionset[tempIndices] <- unname(predict(fit1,TestingData), force = FALSE)
      } else if (RegMethod == 2){
      fitmethod1_1 <- regsubsets(formula1, TrainingData, method = "forward")
      which.max(summary(fit1)$adjr2)
      plotData <- data.frame(predictors = 1:8, adj_R2 = results$adjr2, Cp = results$cp, BIC = results$bic)
      ggplot(data = plotData, mapping = aes(x = predictors, y = adj_R2)) + geom_line() + geom_point()
      TestingData <- model.matrix(formula1, data = TestingData)
      tempPredictionset[tempIndices] <- unname(TestingData[,names(coef(fitmethod1_1, id = 5))], force = FALSE)
      } else if (RegMethod == 3){
      #alpha = 1 (lasso), alpha = 0 (Ridge regression)
      fit1 <- cv.glmnet(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))]), 
                        alpha = 1, family = "gaussian", type.measure = "mse" , nfolds = 10)
      tempPredictionset[tempIndices] <- c(predict(fit1, as.matrix(TestingData), s = fit1$lambda.min))
      } else if (RegMethod == 4){
      #tries <- tuneRF(TrainingData[,-length(colnames(TrainingData))], TrainingData[,length(colnames(TrainingData))],
      #                   ntreeTry=50, stepFactor=2, improve=0.05, plot = F)
      fit1 <- randomForest(TrainingData[,-length(colnames(TrainingData))], TrainingData[,length(colnames(TrainingData))], ntree=250, mtry = 5, nodesize = 5)
      tempPredictionset[tempIndices] <- predict(fit1, newdata = TestingData)
      } else if (RegMethod == 5){
        gbmTrees <- 100
        fit1 <- gbm(CSI_obs~., data=TrainingData, distribution="gaussian", n.trees=gbmTrees, shrinkage=0.05,
                                  interaction.depth=1, bag.fraction = 0.75, train.fraction = 1, n.minobsinnode = 5, verbose = F)
        #best.iter <- gbm.perf(fit1,method="cv", plot.it = F)
        tempPredictionset[tempIndices] <- predict(fit1, newdata = TestingData, n.trees = gbmTrees)
      } else if (RegMethod == 6){
        #maxs <- apply(TrainingData, 2, max)
        #mins <- apply(TrainingData, 2, min)
        #TrainingData <- as.data.frame(scale(TrainingData, center = mins, scale = maxs - mins))
        #TestingData <- as.data.frame(scale(TestingData, center = mins[-length(mins)], scale = maxs[-length(mins)] - mins[-length(mins)]))
        
        #fit1 <- neuralnet(formula3, data = TrainingData, hidden = c(1), threshold = 0.001, stepmax = 1e+06, rep = 1, algorithm = "rprop+")
        #tempPredictionset[tempIndices] <- c(compute(fit1,TestingData)$net.result*
        #                     (maxs[length(maxs)]-mins[length(mins)])+mins[length(mins)])
        
        fit1 <- qrnn.fit(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))])
                 , n.hidden = 1, tau = c(0.5), iter.max=50, n.trials=1, bag=T, trace = F, Th = linear, Th.prime = linear.prime)
        tempPredictionset[tempIndices] <- qrnn.predict(as.matrix(TestingData), fit1)
      } else if (RegMethod == 7){
        fit1 <- svm(formula1, data = TrainingData, cost = 20, gamma = 0.0001, eps = 0.1)
        tempPredictionset[tempIndices] <- predict(fit1, newdata = TestingData)
      } else if (RegMethod == 8){
        fit1 <- rpart(formula1, data = TrainingData)
        tempPredictionset[tempIndices] <- predict(fit1, newdata = TestingData)
      } else if (RegMethod <- 9){
        #preProcess(TrainingData[-length(colnames(TrainingData))], method = c("center", "scale"))
        fitGrid <-  expand.grid(n.trees = c(50,100,150,200), interaction.depth = c(1,2,3), shrinkage = c(0.025,0.05,0.1,0.2), n.minobsinnode = c(10,20))
        fitcontrol <- trainControl(method = "repeatedcv", number = 10, p = 0.75)
        fit1 <- train(TrainingData[-length(colnames(TrainingData))], TrainingData[length(colnames(TrainingData))][[1]], method = paste0(RegMethodTexts[CaretRegMethod]), 
                      trControl = fitcontrol, tuneGrid = NULL, preProcess = NULL)
        tempPredictionset[tempIndices] <- predict(fit1, newdata = TestingData)
      }
    } else if (ProbMethod == 2){
      if (RegMethod == 0){
        tempPredictionset[tempIndices,] <- array(filter(RadiationData, Date %in% init_dates_temp)[[4]][tempIndices],c(length(tempIndices),length(quants)))
      } else if (RegMethod == 1){
        fitcontrol <- gamlss.control(c.crit = 0.02, n.cyc = 200, mu.step = 1, sigma.step = 1, nu.step = 1, tau.step = 1)
        fitcontrol2 <- glim.control(cc = 0.02, cyc = 200, glm.trace = F, bf.cyc = 200, bf.tol = 0.02, bf.trace = F)
        if (DistMethod == 1){
          fit0 <- gamlss(formula = CSI_obs~1, data = TrainingData, family = NO(mu.link = "identity", sigma.link = "identity"), control = fitcontrol,i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = 5, parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmaPredictors[time], parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
          Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.0001, c(length(tempIndices))))
          for (i in 1:length(tempIndices)){
            tempPredictionset[tempIndices[i],] <- unname(qNO(quants, 
                                       mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
          }
        } else if (DistMethod == 2){
          gen.trun(par=c(min(RadiationData[predictand][[1]][which(repeated_init_dates %in% init_dates_temp)], na.rm = T),
                         max(RadiationData[predictand][[1]][which(repeated_init_dates %in% init_dates_temp)], na.rm = T)),"NO", type="both")
          fit0 <- gamlss(formula = CSI_obs~1, data = TrainingData, family = NOtr(mu.link = "identity", sigma.link = "identity"), control = fitcontrol,i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = 5, parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmaPredictors[time], parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
          Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE),array(0.0001, c(length(tempIndices))))
          for (i in 1:length(tempIndices)){
            tempPredictionset[tempIndices[i],] <- unname(qNOtr(quants, 
                                      mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
          }
        } else if (DistMethod == 3){
          gen.trun(par=c(0),"NO", type="left")
          fit0 <- gamlss(formula = CSI_obs~1, data = TrainingData, family = NOtr(mu.link = "identity", sigma.link = "identity"), control = fitcontrol,i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit0, what="mu", scope = list(upper = form), steps = 5, parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          fit1 <- stepGAIC(fit1, what="sigma", scope = list(upper = form), steps = sigmaPredictors[time], parallel = "multicore", ncpus = 4, direction = "both", control = fitcontrol, i.control = fitcontrol2, trace = F)
          Predictionset_mu <- unname(predict(fit1, newdata = TestingData, what = c("mu")), force = FALSE)
          Predictionset_sigma <- pmax(unname(predict(fit1, newdata = TestingData, what = c("sigma")), force = FALSE), array(0.0001, c(length(tempIndices))))
          for (i in 1:length(tempIndices)){
            tempPredictionset[tempIndices[i],] <- unname(qNOtr(quants, 
                                       mu = Predictionset_mu[i], sigma = Predictionset_sigma[i]))
          }
        }
      } else if (RegMethod == 3){
        fit1 <- rq.lasso.fit.mult(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))]), 
                          tau_seq=quants,lambda = 0.2)
        for (q in 1:length(quants)){
          tempPredictionset[tempIndices,q] <- c(predict(fit1, as.matrix(TestingData))[[q]])
        }
      } else if (RegMethod == 4){
        fit1 <- quantregForest(TrainingData[,-length(colnames(TrainingData))], TrainingData[,length(colnames(TrainingData))], nodesize=5,sampsize=25)
        tempPredictionset[tempIndices,] <- predict(fit1, newdata = TestingData, what = quants)
      } else if (RegMethod == 5){
        gbmTrees      <- 200
        for (q in 1:length(quants)){
          fit1 <- gbm(CSI_obs~., data=TrainingData, distribution=list(name = "quantile",alpha = quants[q]), n.trees = gbmTrees, shrinkage=0.05,
                      interaction.depth=1, bag.fraction = 0.75, train.fraction = 1, n.minobsinnode = 5, verbose=FALSE)
          tempPredictionset[tempIndices,q] <- predict(fit1, newdata = TestingData, n.trees = gbmTrees)
        }
      } else if (RegMethod == 6){
        #for (q in 1:length(quants)){
        #  fit1 <- qrnn.fit(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))]),
        #          n.hidden = 1, tau = quants[q], iter.max=50, n.trials=1, bag=T, trace = F, Th = linear, Th.prime = linear.prime)
        #  tempPredictionset[tempIndices,q] <- qrnn.predict(as.matrix(TestingData), fit1)
        #}

        fit1 <- mcqrnn.fit(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))]),
                          n.hidden = 1, n.hidden2 = 1, tau = quants, iter.max=50, n.trials=1,  trace = F, Th = sigmoid, Th.prime = sigmoid.prime)
        tempPredictionset[tempIndices,] <- as.matrix(mcqrnn.predict(as.matrix(TestingData), fit1))
      } else if (RegMethod == 7){
        fit1 <- qtSVM(as.matrix(TrainingData[,-length(colnames(TrainingData))]), as.matrix(TrainingData[,length(colnames(TrainingData))]), weights = quants[-length(quants)], max_gamma = 625, min_lambda  = 3.2e-07, clipping = -1,do.select = TRUE)
        tempPredictionset[tempIndices,] <- cbind(predict(fit1, TestingData), predict(fit1, TestingData)[,length(quants)-1])
      }
      tempPredictionset[tempIndices,] <- t(apply(t(tempPredictionset[tempIndices,]),2,sort))
    }
    if (RegMethod > 0){
      saveRDS(fit1, file = paste0("/usr/people/bakker/kilianbakker/Data/fits/fit_", CVmonth, "_", season, "_", time, "_", ContMethod, "_", ProbMethod, "_", RegMethod,"_", DistMethod, ".rds"))
      rm(fit1, TrainingData, TestingData)
    }
    } 
  }
  DateIndices2 <- which(repeated_init_dates[DateIndices] %in% init_dates_temp)
  if (ProbMethod == 1){
    Predictionset[DateIndices2,time] <- tempPredictionset*CSR[DateIndices2]
  } else if (ProbMethod == 2){
    Predictionset[DateIndices2, ,time] <- tempPredictionset*array(CSR[DateIndices2], c(length(CSR[DateIndices2]), length(quants)))
  }
}
rm(variableData, All_Data, RadiationData)
}
saveRDS(Predictionset, file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets/PS_year_",ContMethod, "_", ProbMethod, "_", RegMethod, "_", DistMethod,".rds"))
rm(Predictionset)
gc()
}