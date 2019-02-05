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
library(verification)
library(xtable)
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

plotcolors <- c("black", "blue", "red", "green", "orange", "purple", "steelblue1", "darkgreen", "yellowgreen")
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

variableNames2 <- c("Global", "Direct (SURF)", "Direct (TOA)", "NCS (SURF)", "NCS (TOA)", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZen", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

variableNames <- c("Global", "Direct_SURF", "Direct_TOA", "NCS_SURF", "NCS_TOA", "CC_Total", "CC_Low", "CC_Medium", "CC_High", "CW_Total", "CW_Low", 
                   "CW_Medium", "CW_High", "PW_Total", "PW_Low", "PW_Medium", "PW_High", "Rain", "AOD_500", "Ang_exp", "Ozone", "T_Low", "T_Medium", "T_High", 
                   "RH_Low", "RH_Medium", "RH_High", "CosZenith", "Lat", "Lon", "DistToCoast", "DistToWater", "DistToInland", "DoY")

source("/usr/people/bakker/kilianbakker/R/functions.R")

###########################
##Calculating importances##
leadTimes <- c(5:19, 29:43)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
NumberofCV <- 3
RegMethod <- 6
DistMethod <- 0
MuSig <- c(0,1)
gbmTrees <- 100
tempImp <- CalculateImp(RegMethod,DistMethod,leadTimes,seasons,NumberofCV,variableNames, MuSig, quants, gbmTrees)
saveRDS(tempImp, file = paste0("/nobackup/users/bakker/Data2/predictionsets2/Imp_final_",RegMethod,"_",DistMethod,".rds"))
###########################

settings <- list(c(2,2,0,0),c(2,2,1,1),c(2,2,1,2),c(2,2,8,0),c(2,2,4,0),c(2,2,9,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GA", "NOTR", "QR", "QRF", "GRF", "GBM", "QRNN")

for (i in 2:3){
tempsetting <- i
tempMethodsTexts <- MethodsTexts[tempsetting]

RegMethod <- settings[[tempsetting]][3]
DistMethod <- settings[[tempsetting]][4]

Importances <- readRDS(file = paste0("/nobackup/users/bakker/Data2/predictionsets2/Imp_final_",RegMethod,"_",DistMethod,".rds"))
tempImportances <- array(0, c(length(variableNames2)))
for (var in 1:length(variableNames2)){
  tempImportances[var] <- mean(Importances[,,,var], na.rm = T)
}
  
plotData <- data.frame(Names =  variableNames2[order(tempImportances)], Importance = tempImportances[order(tempImportances)])
plotData$Names <- factor(plotData$Names,levels = variableNames2[order(tempImportances)])
tempPlot <- ggplot(data = plotData) +
  geom_bar(aes(x = Names, y = Importance), stat = "identity") +
  coord_flip() + xlab("predictors") + ylab("Importance") + theme_bw()
print(tempPlot)
ggsave(paste0("importances_",tempMethodsTexts,".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/final_plots/",  width = 12.5, height = 10, units = "cm")
}

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
seasonTexts <- c( "year", "winter", "spring", "summer", "fall")
fitMethodTexts <- c("seasons-fit", "year-fit")
leadTimesDays <- list(c(5:19), c(29:43))
TimeIndices <- list(c(1:15),c(16:30))
leadTimes <- c(leadTimesDays[[1]],leadTimesDays[[2]])
leadTimesTexts <- c(leadTimesDays[[1]],"-",leadTimesDays[[2]])
settings <- list(c(2,2,0,0),c(2,2,1,1),c(2,2,1,2),c(2,2,8,0),c(2,2,4,0),c(2,2,9,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GA", "NOTR", "QR", "QRF", "GRF", "GBM", "QRNN")
plotQuantiles <- c(0.02,0.5,0.98)
stationPlaces <- c(1:30)

tempsetting <- c(2)
tempMethodsTexts <- MethodsTexts[tempsetting]

RegMethod <- settings[[tempsetting]][3]
DistMethod <- settings[[tempsetting]][4]

tempforecasts <- readRDS(file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_final_", RegMethod, "_", DistMethod, ".rds"))
Rawforecasts <- readRDS(file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_final_0_0.rds"))

tempPlaces <- 23
tempTimes <- c(1:15)
tempDates <- 20171017
DateIndices <- which(init_dates %in% tempDates)
Labels <- c("Obs", "Raw", tempMethodsTexts)
Tempplotcolors <- plotcolors[2:4]
names(Tempplotcolors) <- Labels
observations <- array(NA, c(length(DateIndices), length(tempTimes)))
rawforecast <- Rawforecasts[DateIndices, tempPlaces, tempTimes, which(quants == 0.5)]
for (date in 1:length(DateIndices)){
  observations[date,] <- filter(observationData, Date %in% tempDates[date], Time %in% leadTimes[tempTimes], Station %in% stationData[[2]][tempPlaces])[[4]]
}
CImiddle <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == 0.5)]
CImin <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == plotQuantiles[1])]
CImax <- tempforecasts[DateIndices,tempPlaces,tempTimes,which(quants == plotQuantiles[length(plotQuantiles)])]

if (length(DateIndices) > 1){
  observations <- apply(observations,2,mean, na.rm = T)
  rawforecast <- apply(rawforecast,2,mean, na.rm = T)
  CImiddle <- apply(CImiddle,2,mean, na.rm = T)
  CImin <- apply(CImin,2,mean, na.rm = T)
  CImax <- apply(CImax,2,mean, na.rm = T)
}

plotData <- data.frame(time = leadTimes[tempTimes], data = c(observations, rawforecast, CImiddle), Method = rep(Labels, each = length(tempTimes)))
tempPlot <- ggplot(data = plotData) + 
  geom_line(mapping = aes(x = time, y = data, color = Method)) + 
  scale_colour_manual(values =  Tempplotcolors) + 
  geom_ribbon(data = data.frame(time = leadTimes[tempTimes], CImin, CImax), mapping = aes(x = time, ymin = CImin, ymax = CImax), fill="gray", alpha="0.5") +
  xlab("time") + ylab("Radiation (W/m2)") + theme_bw() + ylim(0,600)

#tempPlot <- ggarrange(tempPlot[[1]], tempPlot[[2]], tempPlot[[3]], ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
print(tempPlot)
ggsave(paste0(tempMethodsTexts, "_", tempDates, ".pdf"), plot = tempPlot, device = "pdf", 
       path = "/usr/people/bakker/kilianbakker/plots/final_plots/", width = 12.5, height = 10, units = "cm")

#################################
##TESTING FOR THE (SKILL) SCORES
#################################
#clearskyHours <-readRDS("/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
#clearskyHours <- clearskyHours[,which(as.numeric(colnames(clearskyHours)) %in% init_dates)]
#LargeErrorHours <-readRDS("/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
#LargeErrorHours <- LargeErrorHours[,which(as.numeric(colnames(LargeErrorHours)) %in% init_dates)]
#init_dates_CS <- c(20160605,20160817,20160824,20160825,20160913,20160914,20161005,20161128,20170409,20170525,20170526,20171015,20180207,20180225)
#init_dates_LE <- c(20160929,20170414,20170516,20170624,20170930,20171017,20171024,20171031,20171120,20171224,20171231,20180123)

stepsize <- 1/50
quants <- seq(stepsize,1-stepsize, stepsize)
seasons <- list(c(201612,201701,201702,201712,201801,201802), c(201604,201605,201703,201704,201705,201803),
                c(201606,201607,201608,201706,201707,201708), c(201609,201610,201611,201709,201710,201711))
seasonTexts <- c( "year", "winter", "spring", "summer", "fall")
leadTimes <- c(5:19,29:43)
leadTimesTexts <- c(5:19,"-",29:43)
TimeIndices <- list(c(1:15),c(16:30))
hyperparamTimes <- c(4,8,12,19,23,27)
hyperparamLeadtimes <- c(8,12,16,32,36,40)
stationPlaces <- c(1:30)
settings <- list(c(2,2,0,0),c(2,2,1,1),c(2,2,1,2),c(2,2,8,0),c(2,2,4,0),c(2,2,9,0),c(2,2,5,0),c(2,2,6,0))
MethodsTexts <- c("RAW", "GA", "NOTR", "QR", "QRF", "GRF", "GBM", "QRNN")

varSizes <- c("1x1gridbox_1lt","1x1gridbox_3lt","3x3gridbox_1lt","3x3gridbox_3lt","5x5gridbox_1lt","5x5gridbox_3lt",
              "7x7gridbox_1lt","7x7gridbox_3lt","9x9gridbox_1lt", "9x9gridbox_3lt")[10]
hyperparamTexts <- c("samps","mtries","nodesizes","ntrees","shrinkages","depths","hiddens1","hiddens2","iters","musteps","sigmasteps","penalties")
hyperparameters <- list(c(2/6,3/6,4/6)[2], c(2/6,3/6,4/6)[1], c(5,10,50)[1], c(100,500,2000)[2], c(0.01,0.1,1)[2], c(1,5,20)[1], c(1,2,3)[1], 
                        c(0,1,2,3)[1], c(2,10,100)[2], c(1,3,5,10)[3], c(0,1,3,5)[2], c(0,0.01,0.1,1))
readingText <- paste0(hyperparameters[[1]],"_",hyperparameters[[2]],"_",hyperparameters[[3]],"_",hyperparameters[[4]],"_",hyperparameters[[5]],"_",
                      hyperparameters[[6]],"_",hyperparameters[[7]],"_",hyperparameters[[8]],"_",hyperparameters[[9]],"_",hyperparameters[[10]],"_",
                      hyperparameters[[11]],"_",hyperparameters[[12]])
#readingText <- paste0(hyperparameters[[1]],"_",hyperparameters[[2]],"_",hyperparameters[[3]],"_",rep(hyperparameters[[4]],each = 9),"_",
#                      rep(rep(hyperparameters[[5]], each = 3), 3),"_",rep(hyperparameters[[6]],9),"_",hyperparameters[[7]],"_",
#                      hyperparameters[[8]],"_",hyperparameters[[9]],"_",hyperparameters[[10]],"_",hyperparameters[[11]],"_",hyperparameters[[12]])
#readingText <- paste0(hyperparameters[[1]],"_",hyperparameters[[2]],"_",hyperparameters[[3]],"_",hyperparameters[[4]],"_",hyperparameters[[5]],"_",
#                      hyperparameters[[6]],"_",rep(hyperparameters[[7]],each = 9),"_",rep(rep(hyperparameters[[8]], each = 3),3),"_",rep(hyperparameters[[9]],9),"_",hyperparameters[[10]],"_",
#                      hyperparameters[[11]],"_",hyperparameters[[12]])
#LabelText <- paste0(rep(hyperparameters[[7]],each = 9),"_",
#                    rep(rep(hyperparameters[[8]], each = 3), 3),"_",rep(hyperparameters[[9]],9))
paramNumber <- 12
temphyperparameters <- hyperparameters[[paramNumber]]
hyperparamText <- hyperparamTexts[paramNumber]
fitmethods <- list(readingText, varSizes, c(paste0(varSizes),paste0(varSizes,"_year")), c("final"))[[1]]
regMethods <- c(8)
tempMethodsTexts <- MethodsTexts[regMethods]
seasonOptions <- c(1:length(seasonTexts))[1]
Places <- c(1:length(stationPlaces))
PlaceOption <- c("group","separate")[1]
Times <- c(1:length(leadTimes))
thresholds <- list(seq(0.1,0.9,0.1),c(0.2,0.5,0.9),c(0))[[3]]
extraOptions <- list(seq(0.01,1,0.01),seq(0.1,1,0.1),c(1))[[3]]

ScoringMethods <- list(c(0,"RMSE"),c(0,"MAE"),c(0,"CRPS"),c(1,"BS"),c(1,"PEV"),c(1,"ReliabilityPlot"))[3]
scoreSettings <- c(1,3)[2]
Percs <- c(1,2,3)[1]
scaledPredictand <- 1

Scores <- array(NA,c(length(Times),length(Places),length(seasonOptions),length(regMethods)+1,length(fitmethods),
                     length(thresholds),length(extraOptions),length(ScoringMethods),length(scoreSettings),length(Percs)))
CLratiosupd <- array(NA, c(length(thresholds),length(extraOptions)))

for (s1 in 1:length(ScoringMethods)){
  ScoringMethod <- ScoringMethods[[s1]]
for (s2 in 1:length(scoreSettings)){
  scoreSetting <- scoreSettings[s2]
for (s3 in 1:length(Percs)){
  Perc <- Percs[s3]
for (r in 1:length(regMethods)){
  tempregMethod <- regMethods[r]
  for (f in 1:length(fitmethods)){
  tempfitmethod <- fitmethods[f]
  if (fitmethods == readingText[1]){
    tempforecasts <- readRDS(file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_",
                settings[[tempregMethod]][3], "_", settings[[tempregMethod]][4], "_",tempfitmethod,".rds"))
  } else {
    tempforecasts <- readRDS(file = paste0("/nobackup/users/bakker/Data2/predictionsets2/PS_", tempfitmethod,"_",
                                         settings[[tempregMethod]][3], "_", settings[[tempregMethod]][4],".rds"))
  }

  for (time in 1:length(Times)){
    leadTime <- leadTimes[Times[time]]
    if (leadTime <= 24){
      temp_init_dates <- init_dates
    } else {
      temp_init_dates <- c(1:length(init_dates))
      for (d in 1:length(init_dates)){
        tempDate <- paste0(substr(init_dates[d],1,4), "-", substr(init_dates[d],5,6), "-", substr(init_dates[d],7,8))
        temp_init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
      }
    }
    
    for (tempseason in 1:length(seasonOptions)){
      season <- seasonOptions[tempseason]
      if (season == 1){
        temp_init_dates2 <- init_dates
      } else if (season %in% c(2,3,4,5)){
        temp_init_dates2 <- init_dates[floor(init_dates/100) %in% seasons[[season-1]]]
      } 
      DateIndices <- which(init_dates %in% temp_init_dates2)
      
    if (PlaceOption == "group"){
        observations <- filter(observationData, Station %in% stationData[[2]][stationPlaces[Places]], Time == (leadTime %% 24), Date %in% temp_init_dates[DateIndices])[[4]]
        CSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces[Places]], Time == (leadTime %% 24), Date %in% temp_init_dates[DateIndices])[[4]]
        forecasts <- array(NA,c(length(Places)*length(DateIndices),length(quants)))
        for (place in 1:length(Places)){
          forecasts[(length(DateIndices)*(place-1)+1):(length(DateIndices)*place),] <- tempforecasts[DateIndices,Places[place],Times[time],]
        }
        for (t in 1:length(thresholds)){
          tempscores <- calculating_scores(ScoringMethod, observations, forecasts, CSR, quants, thresholds[t], scaledPredictand, Perc)
          Scores[time,1,tempseason,r,f,t,,s1,s2,s3] <- tempscores[[scoreSetting]]
          if (ScoringMethod[2] == "PEV"){
            CLratiosupd[t,] <- tempscores[[2]]
          } else {
            Scores[time,1,tempseason,length(regMethods)+1,f,t,,s1,s2,s3] <- tempscores[[2]]
          }
        }
    } else if (PlaceOption == "separate"){
    for (place in 1:length(Places)){
      stationPlace <- stationPlaces[Places[place]]
      observations <- filter(observationData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates[DateIndices])[[4]]
      CSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), Date %in% temp_init_dates[DateIndices])[[4]]
      forecasts <- tempforecasts[Dateindices,Places[place],Times[time],]
      for (t in 1:length(thresholds)){
        tempscores <- calculating_scores(ScoringMethod, observations, forecasts, CSR, quants, thresholds[t], scaledPredictand, Perc)
        Scores[time,place,tempseason,r,f,t,,s1,s2,s3] <- tempscores[[scoreSetting]]
        Scores[time,place,tempseason,length(regMethods)+1,f,t,,s1,s2,s3] <- tempscores[[2]]
      }
    }
    }
    }
  }
  }
}
}
}
}

#hyperparameters
for (k in 1:length(regMethods)){
for (i in 1:length(seasonOptions)){
  ScoreData <- c()
for (j in 1:length(fitmethods)){
  ScoreData <- c(ScoreData,apply(Scores[hyperparamTimes,,i,k,j,,,,,],1,mean,na.rm = T))
}

plotData <- data.frame(xaxis = rep(seq(hyperparamLeadtimes),length(temphyperparameters)), yaxis = ScoreData, 
                       Value = rep(as.character(round(temphyperparameters, digits = 3)), each = length(hyperparamLeadtimes)))
tempplotcolors <- plotcolors[1:length(temphyperparameters)]
names(tempplotcolors) <- as.character(round(temphyperparameters, digits = 3))
tempPlot <- LinePlot(plotData,tempplotcolors, c(0.1,0.5), "CRPSS", "Lead Time")
print(tempPlot)
ggsave(paste0("CRPSS_",tempMethodsTexts[k],"_", hyperparamText,"_",seasonTexts[i],".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/comparison_plots/hyperparameters/", width = 12.5, height = 10, units = "cm")
}
}

#space/time
for (k in 1:length(regMethods)){
  for (i in 1:length(seasonOptions)){
    ScoreData <- c()
    for (j in 1:length(fitmethods)){
      ScoreData <- c(ScoreData,apply(Scores[TimeIndices[[1]],,i,k,j,,,,,],1,mean,na.rm = T),NA,apply(Scores[TimeIndices[[2]],,i,k,j,,,,,],1,mean,na.rm = T))
    }
    
    plotData <- data.frame(xaxis = rep(seq(leadTimesTexts),length(varSizes)), yaxis = ScoreData, 
                           Size = rep(varSizes, each = length(leadTimesTexts)))
    tempplotsettings <- list(plotsettings[[1]][c(1:length(varSizes))],plotsettings[[2]][c(1:length(varSizes))])

    tempPlot <- LinePlot2(plotData, tempplotsettings, c(0.2,0.6), "CRPSS", "Lead Time")
    print(tempPlot)
    ggsave(paste0("CRPSS_spacetimecomparison_",tempMethodsTexts[k],"_",seasonTexts[i],".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/comparison_plots/", width = 12.5, height = 10, units = "cm")
  }
}

#yearfit/seasonfit
for (k in 1:length(regMethods)){
  for (i in 1:length(seasonOptions)){
    ScoreData <- c()
    for (j in 1:length(fitmethods)){
      ScoreData <- c(ScoreData,apply(Scores[TimeIndices[[1]],,i,k,j,,,,,],1,mean,na.rm = T),NA,apply(Scores[TimeIndices[[2]],,i,k,j,,,,,],1,mean,na.rm = T))
    }
    
    plotData <- data.frame(xaxis = rep(seq(leadTimesTexts),2), yaxis = ScoreData, 
                           Fit = rep(c("per season","complete dataset"), each = length(leadTimesTexts)))
    
    tempplotcolors <- plotcolors[1:2]
    names(tempplotcolors) <- c("per season","complete dataset")
    tempPlot <- LinePlot(plotData,tempplotcolors, c(0.2,0.6), "CRPSS", "Lead Time")
    print(tempPlot)
    ggsave(paste0("CRPSS_fitcomparison_",tempMethodsTexts[k],"_",seasonTexts[i],".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/comparison_plots/", width = 12.5, height = 10, units = "cm")
  }
}

#final plots
tempMethodsTexts <- c(tempMethodsTexts,"CLIM") #Scores or skillscores?
tempMethodsTextsExt <- paste0(rep(tempMethodsTexts, each = 2), rep(c("_RMSE","_MAE"), length(tempMethodsTexts)))
tempplotcolors <- plotcolors[1:length(tempMethodsTexts)]
names(tempplotcolors) <- tempMethodsTexts
plotSort <- 1
timetemp <- 1
placetemp <- 1
seasontemp <- 1
threstemp <- 1
for (seasontemp in 1:length(seasonOptions)){
if (plotSort == 1){
  if (seasontemp == 1){
    temp_init_dates <- init_dates
  } else if (seasontemp %in% c(2,3,4,5)){
    temp_init_dates <- init_dates[floor(init_dates/100) %in% seasons[[seasontemp-1]]]
  } 
  DateIndices <- which(init_dates %in% temp_init_dates)
    
  tempobservations <- filter(observationData, Station %in% stationData[[2]][stationPlaces], Time %in% c(5:19), Date %in% init_dates[DateIndices])[[4]]
  tempCSR <- filter(clearskyData, Station %in% stationData[[2]][stationPlaces], Time %in% c(5:19), Date %in% init_dates[DateIndices])[[4]]
  tempObs <- (tempobservations/tempCSR)[which(tempCSR >= 20)]

  ScoreData <- c()
  for (r in 1:(length(regMethods)+1)){
    ScoreData <- c(ScoreData, Scores[TimeIndices[[1]],placetemp,seasontemp,r,,,,1,,1], NA, Scores[TimeIndices[[2]],placetemp,seasontemp,r,,,,1,,1],
                   Scores[TimeIndices[[1]],placetemp,seasontemp,r,,,,2,,1], NA, Scores[TimeIndices[[2]],placetemp,seasontemp,r,,,,2,,1])
  }
  plotData <- data.frame(xaxis = rep(seq(leadTimesTexts), length(tempMethodsTextsExt)), yaxis = ScoreData,
                         Method = rep(tempMethodsTextsExt, each = length(leadTimesTexts)))
  tempplotsettings <- list(plotsettings[[1]][c(1:length(tempMethodsTextsExt))],plotsettings[[2]][c(1:length(tempMethodsTextsExt))])
  
  plotName <- paste0("stationAvg_ltAll_", seasonTexts[seasontemp])
  tempPlot <- LinePlot2(plotData, tempplotsettings, c(0.1,0.3), "Score", "leadTime") +  
                                scale_y_continuous(sec.axis = sec_axis(~.*100/mean(tempObs), name = "Relative Score (%)"))
  print(tempPlot)
} else if (plotSort == 2){
  ScoreData <- c()
  for (r in 1:length(regMethods)){
    ScoreData <- c(ScoreData, Scores[TimeIndices[[1]],placetemp,seasontemp,r,,,,2,,1], NA, 
                Scores[TimeIndices[[2]],placetemp,seasontemp,r,,,,2,,1])
  }
  plotData <- data.frame(xaxis = rep(seq(leadTimesTexts), length(tempMethodsTexts)), yaxis = ScoreData, Method = rep(tempMethodsTexts, each = length(leadTimesTexts)))
  plotName <- paste0("stationAvg_ltAll_", seasonTexts[seasontemp])
  tempPlot <- LinePlot(plotData, tempplotcolors[2:8], c(0.1,0.45), "CRPSS", "leadTime")
  print(tempPlot)
} else if (plotSort == 3){
  ScoreData <- round(apply(Scores[timetemp,,seasontemp,,,,,,,,],1,max, na.rm = T), digits = 2)
  tempMethodsTexts2 <- array(NA, c(length(stationPlaces)))
  for (p in 1:length(stationPlaces)){
    tempMethodsTexts2[p] <- tempMethodsTexts[which(max(Scores[timetemp,p,seasontemp,,,,,,,],na.rm = T) == Scores[timetemp,p,seasontemp,,,,,,,,])]
  }
  plotData <- data.frame(longitude = stationData[[3]], latitude = stationData[[4]], Method = tempMethodsTexts2, Labels = ScoreData)
  plotName <- paste0(ScoringMethod[2], "_", "NL_lt", leadTimes[timetemp], "_", seasonTexts[[seasontemp]])
  tempPlot <- MapPlot(plotData, tempplotcolors)
  print(tempPlot)
} else if (plotSort == 4){
  ScoreData <- c()
  for (r in 1:length(regMethods)){
    ScoreData <- c(ScoreData, Scores[timetemp,placetemp,seasontemp,r,,threstemp,,,,])
  }
  Binaxis <- seq(1/(2*length(NumberofBins)),1-1/(2*length(NumberofBins)),1/length(NumberofBins))
  plotData <- data.frame(xaxis = rep(Binaxis, length(tempMethodsTexts)), yaxis = ScoreData, Method = rep(tempMethodsTexts, each = length(NumberofBins)))
  plotName <- paste0("stationAvg_lt", leadTimes[Times], "_", seasonTexts[[1]], "_", thresholds[threstemp])
  tempPlot <- Reliability_plot(plotData, tempplotcolors[2:8])
  print(tempPlot)
} else if (plotSort == 5){
  ScoreData <- c()
  for (r in 1:length(regMethods)){
    ScoreData <- c(ScoreData, Scores[timetemp,placetemp,seasontemp,r,,threstemp,,,,])
  }
  plotData <- data.frame(xaxis = rep(CLratiosupd[threstemp,], length(tempMethodsTexts)), yaxis = ScoreData, 
                         Method = rep(tempMethodsTexts, each = length(CLratiosupd[threstemp,])))
  plotName <- paste0("stationAvg_lt", leadTimes[Times], "_", seasonTexts[[1]], "_", thresholds[threstemp])
  tempPlot <- PEVplot(plotData, tempplotcolors,c(-0.1,0.7))
  print(tempPlot)
} else if (plotSort == 6){
  ScoreData <- c()
  for (r in 1:length(regMethods)){
    ScoreData <- c(ScoreData, Scores[timetemp,placetemp,seasontemp,r,,,,,,])
  }
  plotData <- data.frame(xaxis = rep(thresholds, length(tempMethodsTexts)), yaxis = ScoreData, Method = rep(tempMethodsTexts, each = length(thresholds)))
  plotName <- paste0("stationAvg_lt", leadTimes[Times], "_", seasonTexts[[1]])
  tempPlot <- LinePlot(plotData, tempplotcolors[2:8], c(-0.2,0.45), "BSS", "Threshold")
  print(tempPlot)
}
#ggsave(paste0("CRPSS_final_", plotName,".pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/final_plots/", width = 12.5, height = 10, units = "cm")
}