# set library paths:
libpathKiri = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
libpath = "/usr/people/bakker/R/x86_64-redhat-linux-gnu-library/3.4"
#libpath = "/tmp/RtmpSDMBdh/downloaded_packages"
.libPaths(c(.libPaths()))#, libpathKiri))
library(tidyverse)
library(rasterVis)
library(ggplot2)
library(ggpubr)
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

tmprds <- readRDS(file = "/nobackup/users/bakker/Data/temperaturevariables/20160401.rds")
xcoordinate <- coordinatesLargeGrid(stationData,stationPlace, c(0), c(0), tmprds)[[1]]
ycoordinate <- coordinatesLargeGrid(stationData,stationPlace, c(0), c(0), tmprds)[[2]]
xcoordinate2 <- coordinatesSmallGrid(stationData,stationPlace)[[1]]
ycoordinate2 <- coordinatesSmallGrid(stationData,stationPlace)[[2]]

plotcolors <- c("blue", "red", "green", "orange", "purple", "black", "steelblue1")
plotlines <- c("solid", "dashed")
plotsettings <- list(rep(plotcolors, each = length(plotlines)), rep(plotlines, length(plotcolors)))

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
##TESTING FOR THE VARIABLES
#################################
variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

AllvariableData <- readRDS(file  = "/usr/people/bakker/kilianbakker/Data/variables_3x3gridbox_3lt.rds")

variableNumber <- 21
Times <- c(1:15)
stationPlace <- 23
dateNumbers <- c(20160401,20160402,20160403,20160404,20160505)
dateIndices <- which(init_dates %in% dateNumbers)
Data <- c()
for (i in 1:length(dateIndices)){
  Data <- c(Data, AllvariableData[dateIndices[i],stationPlace, Times,variableNumber])
}
plotData <- data.frame(time = Times, data = Data, type = rep(as.character(dateNumbers), each = length(Times)))
tempPlot <- ggplot(data = plotData) + 
  geom_line(mapping = aes(x = time, y = data, color = type)) +
  scale_colour_manual(values = plotcolors[1:length(dateNumbers)]) +
  xlab("Time") + ylab(variableNames[variableNumber]) + ggtitle(paste0(variableNames[variableNumber]," comparison"))

ggsave(paste0(variableNames[variableNumber],"_comparison.pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/variable_plots/")

#################################
##TESTING FOR REGRESSION BETWEEN VARIABLES
#################################
variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

AllvariableData <- readRDS(file  = "/usr/people/bakker/kilianbakker/Data/variables_3x3gridbox_3lt.rds")

clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskies <- which(clearskyHours == 1, arr.ind = T)

variableNumber <- 21
Times <- c(1:15)
leadTimes <- c(5:19,29:43)
stationPlace <- 23
tempclearskies <- clearskies[(clearskies[,1] %in% leadTimes[Times]) & (clearskies[,2] <= length(init_dates)),]
observations <- c()
forecasts <- c()
variableData <- c()
for (cs in 1:length(tempclearskies[,1])){
  date <- unname(tempclearskies[cs,2])
  time <- Times[which(leadTimes == unname(tempclearskies[cs,1]))]
#for (date in 1:length(init_dates)){
#  for (time in Times){
    observations <- c(observations, filter(observationData, Date %in% init_dates[date], Station == stationData[[2]][stationPlace], Time %in% leadTimes[time])[[4]])
    forecasts <- c(forecasts, AllvariableData[date,stationPlace,time,1])
    variableData <- c(variableData, AllvariableData[date,stationPlace,time,variableNumber])
#  }
#}
}

plotData <- data.frame(RadiationError = observations - forecasts, variableData)
tempPlot <- FitPlot(plotData, 0.25, -200, 40) + xlab(variableNames[variableNumber]) + ylab("Radiation Error") + ggtitle(paste0(variableNames[variableNumber]," scatterplot"))

ggsave(paste0(variableNames[variableNumber],"_scatterplot_CSdays.pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/regression_plots/particles/")

#################################
##TESTING FOR IMPORTANCES
#################################

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

Method <- "GBM"
Importances1 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/",Method,"imp.rds"))
Importances2 <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/",Method,"imp_v2.rds"))
tempImportances <- array(0, c(length(variableNames),2))
for (var in 1:length(variableNames)){
  tempImportances[var,1] <- sum(Importances1[,,,,var], na.rm = T)
  tempImportances[var,2] <- sum(Importances2[,,,,var], na.rm = T)
}
  
version <- 2
plotData <- data.frame(Names =  variableNames[order(tempImportances[,version])], Importance = tempImportances[,version][order(tempImportances[,version])])
plotData$Names <- factor(plotData$Names,levels = variableNames[order(tempImportances[,version])])
tempPlot <- ggplot(data = plotData) +
  geom_bar(aes(x = Names, y = Importance), stat = "identity") +
  coord_flip() + xlab("predictors") + ylab("Importance") + ggtitle(paste0("Importance measure from ", Method))
print(tempPlot)
ggsave(paste0("importances_",Method, "_v",version,".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/")

#################################
##TESTING FOR THE PREDICTION SETS
#################################
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
LargeErrorHours <-readRDS("/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
LargeErrorHours <- LargeErrorHours[,which(as.numeric(colnames(LargeErrorHours)) %in% init_dates)]

stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
seasonTexts <- c("winter", "spring", "summer", "fall", "year")
fitMethodTexts <- c("seasons-fit", "year-fit")
leadTimesDays <- list(c(5:19), c(29:43))
TimeIndices <- list(c(1:15),c(16:30))
leadTimes <- c(leadTimesDays[[1]],leadTimesDays[[2]])
leadTimesTexts <- c(leadTimesDays[[1]],"-",leadTimesDays[[2]])
settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GLM",  "Lasso", "RF", "GBM", "NN", "RAW", "NO_0", "NO_1", "NO_+", "NOtr[min,max)", "NOtr[0,Inf]", "Lasso", "QRF", "GBM", "QRNN")
plotQuantiles <- c(0.1,0.5,0.9)
stationPlaces <- c(1:30)

tempsetting <- c(15)
tempMethodsTexts <- MethodsTexts[tempsetting]

ContMethod <- settings[[tempsetting]][1]
ProbMethod <- settings[[tempsetting]][2]
RegMethod <- settings[[tempsetting]][3]
DistMethod <- settings[[tempsetting]][4]

tempforecasts <- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants),length(fitMethodTexts)))
tempforecasts[,,,,1] <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets2/PS_allstations_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
tempforecasts[,,,,2] <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets2/PS_allstations_year_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
Rawforecasts <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets2/PS_allstations_2_2_0_0.rds"))

tempPlaces <- 23
tempTimes <- c(7:15)
tempDates <- 20170624
DateIndices <- which(init_dates %in% tempDates)
fit_method <- 1
Labels <- c("Obs", "RAW", tempMethodsTexts)
observations <- array(NA, c(length(DateIndices), length(tempTimes)))
rawforecast <- Rawforecasts[DateIndices, tempPlaces, tempTimes, which(quants == 0.5)]
for (date in 1:length(DateIndices)){
  observations[date,] <- filter(observationData, Date %in% tempDates[date], Time %in% leadTimes[tempTimes], Station %in% stationData[[2]][tempPlaces])[[4]]
}
CImiddle <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == 0.5),fit_method]
CImin <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == plotQuantiles[1]),fit_method]
CImax <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == plotQuantiles[length(plotQuantiles)]),fit_method]

if (length(DateIndices) > 1){
  observations <- apply(observations,2,mean, na.rm = T)
  rawforecast <- apply(rawforecast,2,mean, na.rm = T)
  CImiddle <- apply(CImiddle,2,mean, na.rm = T)
  CImin <- apply(CImin,2,mean, na.rm = T)
  CImax <- apply(CImax,2,mean, na.rm = T)
}

plotData <- data.frame(time = leadTimes[tempTimes], data = c(observations, rawforecast, CImiddle), labels = rep(Labels, each = length(tempTimes)))

tempPlot <- ggplot(data = plotData) +
  geom_line(mapping = aes(x = time, y = data, color = labels)) +
  scale_colour_manual(values = plotcolors[1:length(Labels)]) +
  geom_ribbon(data = data.frame(time = leadTimes[tempTimes], CImin, CImax), mapping = aes(x = time, ymin = CImin, ymax = CImax), fill="gray", alpha="0.5") +
  xlab("time") + ylab("Radiation (W/m2)") + ggtitle(paste0(tempDates, fitMethodTexts[fit_method]))

#tempPlot <- ggarrange(tempPlot[[1]], tempPlot[[2]], tempPlot[[3]], ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
print(tempPlot)
ggsave(paste0(tempMethodsTexts, "_", tempDates, "_comparison_", fit_method, ".pdf"), plot = tempPlot, device = "pdf", 
       path = "/usr/people/bakker/kilianbakker/plots/comparison_plots/", width = 20, height = 14, units = "cm")

#################################
##TESTING FOR THE (SKILL) SCORES
#################################
clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
LargeErrorHours <-readRDS("/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
LargeErrorHours <- LargeErrorHours[,which(as.numeric(colnames(LargeErrorHours)) %in% init_dates)]
init_dates_CS <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170525,20170526,20171015,20180207,20180225)
init_dates_LE <- c(20160929,20170414,20170516,20170624,20170930,20171017,20171024,20171031,20171120,20171224,20171231,20180123)

thresholds <- c(0.8)
NumberofBins <- 10
stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
seasonTexts <- c("winter", "spring", "summer", "fall", "year", "CSdays", "CShours", "LEdays", "LEhours")
fitMethodTexts <- c("seasons-fit", "year-fit")
seasonPlotLimits <- list(c(0.4,0.65),c(0.3,0.65),c(0.2,0.55),c(0.4,0.75),c(0.3,0.75),c(0.25,0.75),c(0.25,0.75),c(0,0.75),c(0,0.75))
leadTimesDays <- list(c(5:19), c(29:43))
TimeIndices <- list(c(1:15),c(16:30))
leadTimes <- c(leadTimesDays[[1]],leadTimesDays[[2]])
leadTimesTexts <- c(leadTimesDays[[1]],"-",leadTimesDays[[2]])
settings <- list(c(2,1,0,0),c(2,1,1,0),c(2,1,3,0),c(2,1,4,0),c(2,1,5,0),c(2,1,6,0),c(2,2,0,0),c(2,2,1,-1),c(2,2,1,0),c(2,2,1,1),c(2,2,1,2),c(2,2,1,3),c(2,2,3,0),c(2,2,4,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GLM",  "Lasso", "RF", "GBM", "NN", "RAW", "NO_0", "NO_1", "NO_+", "NOtr[min,max)", "NOtr[0,Inf]", "Lasso", "QRF", "GBM", "QRNN")

seasonOptions <- c(1:5)
stationPlaces <- c(1:30)
Times <- c(1:length(leadTimes))
Scores <- array(NA,c(length(leadTimes),length(stationPlaces),length(settings),length(seasonOptions),length(thresholds),length(fitMethodTexts),NumberofBins))
Skillscores <- array(NA,c(length(leadTimes),length(stationPlaces),length(settings),length(seasonOptions),length(thresholds),length(fitMethodTexts),NumberofBins))

tempsettings <- c(7,10:12,14:16)
tempMethodsTexts <- MethodsTexts[tempsettings]
tempPlaces <- c(23)
tempTimes <- c(8)
ScoringMethod <- c(1,"ReliabilityPlot")

for (setting in tempsettings){
  ContMethod <- settings[[setting]][1]
  ProbMethod <- settings[[setting]][2]
  RegMethod  <- settings[[setting]][3]
  DistMethod <- settings[[setting]][4]
  
  tempforecasts<- array(NA, c(length(init_dates), length(stationPlaces), length(leadTimes), length(quants),length(fitMethodTexts)))
  tempforecasts[,,,,1] <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets2/PS_allstations_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
  tempforecasts[,,,,2] <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/predictionsets2/PS_allstations_year_",ContMethod, "_" ,ProbMethod, "_", RegMethod, "_", DistMethod, ".rds"))
  
  for (time in tempTimes){
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
  
    for (place in tempPlaces){
      stationPlace <- stationPlaces[place]
      observations <- filter(observationData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]]
      CSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates)[[4]]
  
      for (season in seasonOptions){
        if (season <= length(seasons)){
          temp_init_dates2 <- init_dates[floor(init_dates/100) %in% seasons[[season]]]
        } else if (season == (length(seasons)+1)){
          temp_init_dates2 <- init_dates
        } else if (season == (length(seasons)+2)){
          temp_init_dates2 <- init_dates_CS #clear sky days
        } else if (season == (length(seasons)+3)){
          temp_init_dates2 <- init_dates[which(clearskyHours[leadTime,] == 1)] #clear sky hours
        } else if (season == (length(seasons)+4)){
          temp_init_dates2 <- init_dates_LE #large error days
        } else if (season == (length(seasons)+5)){
          temp_init_dates2 <- init_dates[which(LargeErrorHours[leadTime,] == 1)] #large error hours
        }
  
        Dateindices <- which(init_dates %in% temp_init_dates2)
        if (length(Dateindices) > 1){
          tempobservations <- observations[Dateindices]
          tempCSR <- CSR[Dateindices]
          clim_scores <- calculating_sampleclimatologies(ContMethod, ProbMethod, tempobservations, tempCSR, temp_init_dates2, quants, discretizeWidth)
  
          for (k in 1:length(fitMethodTexts)){
            for (thres in 1:length(thresholds)){
              threshold <- thresholds[thres]
              Scores[time,place,setting,season,thres,k,] <- calculating_scores(ScoringMethod, tempobservations, tempforecasts[Dateindices,place,time,,k], tempCSR)
              Skillscores[time,place,setting,season,thres,k,] <- 1 - Scores[time,place,setting,season,thres,k,]/clim_scores
              #Skillscores[time,place,setting,season,thres,k] <- 1 - Scores[time,place,setting,season,thres,k]/Scores[time,place,which(MethodsTexts == "RAW")[2],season,k]
            }
          }
        }
      }
    }
  }
}

tempseason <- 5
fit_method <- 1
RegMethod <- 4
plotSort <- 3
saving <- 0
tempthreshold <- 1

if (plotSort == 1){
  SSData <- c()
  for (setting in tempsettings){
    SSData <- c(SSData, Skillscores[TimeIndices[[1]],tempPlaces,setting,tempseason,tempthreshold,fit_method], NA, Skillscores[TimeIndices[[2]],tempPlaces,setting,tempseason,tempthreshold,fit_method])
  }
  plotData <- data.frame(data = SSData, Method = rep(tempMethodsTexts, each = length(leadTimesTexts)))
  plotName <- paste0(ScoringMethod[2], "_", stationData[[1]][tempPlaces], "_ltAll_", seasonTexts[tempseason], "_", fitMethodTexts[fit_method])
  tempPlot <- LinePlot(plotData, plotName, plotsettings, leadTimesTexts, c(0,0.7))
  print(tempPlot)
  if (saving == 1){
    ggsave(paste0(plotName, ".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/scores_plots/")
  }
} else if (plotSort == 2){
  tempMethodsTexts2 <- array(tempMethodsTexts[RegMethod], c(length(stationPlaces)))
  for (place in 1:length(stationPlaces)){
    tempMethodsTexts2[place] <- tempMethodsTexts[which(max(Skillscores[tempTimes,place,tempsettings,tempseason,fit_method],na.rm = T) == Skillscores[tempTimes,place,tempsettings,tempseason,fit_method])]
  }
  plotData <- data.frame(latitude = stationData[[4]], longitude = stationData[[3]], Method = tempMethodsTexts2, Labels = round(apply(Skillscores[tempTimes,,,tempseason,fit_method],1,max, na.rm = T), digits = 2))
  plotName <- paste0(ScoringMethod[2], "_", "NL_lt", leadTimes[tempTimes], "_", seasonTexts[[tempseason]], "_", fitMethodTexts[fit_method])
  tempPlot <- MapPlot(plotData, plotName, plotsettings)
  print(tempPlot)
  if (saving == 1){
    ggsave(paste0(plotName, ".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/scores_plots/")
  }
} else if (plotSort == 3){
  ReliabilityData <- c()
  for (setting in tempsettings){
    ReliabilityData <- c(ReliabilityData, Scores[tempTimes, tempPlaces, setting, tempseason, tempthreshold, fit_method,])
  }
  plotData <- data.frame(xaxis = seq(1/(2*NumberofBins),1-1/(2*NumberofBins),1/NumberofBins), data = ReliabilityData, Method = rep(tempMethodsTexts, each = NumberofBins))
  plotName <- paste0(ScoringMethod[2], "_", stationData[[1]][tempPlaces], "_lt", leadTimes[tempTimes], "_", seasonTexts[[tempseason]], "_", fitMethodTexts[fit_method])
  tempPlot <- Reliability_plot(plotData, plotName, plotsettings)
  print(tempPlot)
  if (saving == 1){
    ggsave(paste0(plotName, "_t", thresholds[tempthreshold]*100, ".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/scores_plots/")
  }
}
