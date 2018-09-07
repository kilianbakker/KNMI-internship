ClearskyRadiation <- function(ToD, date, lat, lon, timezone, typeCS, typeZenith){
  ###########################
  #VARIABLES
  ###########################
  zenithAngle <- zenithAngle(ToD, date, lat, lon, timezone, typeZenith)
  z     <- zenithAngle*pi/180
  mu    <- cos(z)
  gamma <- 90 - zenithAngle
  I0    <- 1367
  p0    <- 1013
  p     <- 1013
  eps   <- 1
  TLmethod <- 1
  
  month <- as.numeric(gsub(floor(date/10000),"",floor(date/100)))
  if (TLmethod == 1){
    #Climatological monthly average turbidity values taken from https://hal.archives-ouvertes.fr/hal-00465791/document
    Turbidities <- c(2.7,2.5,3.6,3.2,3.4,3.9,4.1,4.6,4.0,3.6,2.4,2.7)
    TL <- Turbidities[month]
  } else if (TLmethod == 2){
    #Climatological monthly average turbidity values taken from http://www.meteonorm.com/images/uploads/downloads/ieashc36_report_TL_AOD_climatologies.pdf
    Turbidities <- c(2.85,3.54,4.43,4.12,4.33,4.34,4.16,3.95,4.03,3.79,3.12,3.43)
    TL <- Turbidities[month]
  } else if (TLmethod == 3){
    #Model for calculating TL from forecasts
    tmpdata1 <- readRDS(file = paste0("/nobackup/users/bakker/Data/particlesvariables/",date,".rds"))
    tmpdata2 <- readRDS(file = paste0("/nobackup/users/bakker/Data/cloudwatervariables/",date,".rds"))
    
    coordinates1   <- unique(tmpdata2[[2]])
    coordinates2   <- unique(tmpdata2[[3]])
    xcoordinate    <- coordinates1[which(min(abs(coordinates1 - lon)) == abs(coordinates1 - lon))]
    ycoordinate    <- coordinates2[which(min(abs(coordinates2 - lat)) == abs(coordinates2 - lat))]
    
    coordinates1   <- unique(tmpdata1[[2]])
    coordinates2   <- unique(tmpdata1[[3]])
    xcoordinate2    <- coordinates1[which(min(abs(coordinates1 - lon)) == abs(coordinates1 - lon))]
    ycoordinate2    <- coordinates2[which(min(abs(coordinates2 - lat)) == abs(coordinates2 - lat))]
    
    closestLT1 <- which(min(abs(ToD - unique(tmpdata1[[1]]))) == abs(ToD - unique(tmpdata1[[1]])))
    closestLT2 <- which(min(abs(ToD - c(4:20))) == abs(ToD - c(4:20)))
    
    tmpdata1 <- filter(tmpdata1, xcoor == xcoordinate2, ycoor == ycoordinate2)
    tmpdata2 <- filter(tmpdata2, xcoor == xcoordinate, ycoor == ycoordinate)
    
    w <- tmpdata2[[9]][[closestLT2]]/10
    beta <- tmpdata1[[4]][[closestLT1]]/(0.5^(-tmpdata1[[5]][[closestLT1]]))
    TL <- 1.8494 + 0.2425*w - 0.0203*w^2 + (15.427 + 0.3153*w - 0.0254*w^2)*beta
  }
  
  ###########################
  #MODEL
  ###########################
  #ESRA CLEAR SKY MODEL
  if (typeCS == 1){
  #solar altitude angle corrected for refraction
  gamma2 <- gamma + 0.06*(180/pi)*(0.16 + 1.123*(pi/180)*gamma + 0.066*(pi/180)^2*gamma^2)/(1 + 28.93*(pi/180)*gamma + 277.4*(pi/180)^2*gamma^2)
  gamma2 <- gamma2*pi/180
  #relative optical air mass
  m <- (p/p0)/(sin(gamma2) + 0.51*(gamma2 + 6.08)^(-1.64))
  #Rayleigh optical thickness
  if (m <= 20){
    delta <- 1/(6.63 + 1.75*m - 0.12*m^2 + 0.0065*m^3 - 0.00013*m^4)
  } else {
    delta <- 1/(10.4 + 0.718*m)
  }
  
  #diffuse transmission function at zenith
  Transmission <- -0.0158 + 0.0305 * TL + 0.00038 * TL^2
  #diffuse angular function
  A0 <- 0.265 - 0.062 * TL + 0.0031 * TL^2
  if (A0 < 0.002/Transmission){
    A0 <- 0.002/Transmission
  }
  A1 <- 2.04 + 0.0189 * TL - 0.0112 * TL^2
  A2 <- -1.3 + 0.039 * TL + 0.0085 * TL^2
  Fang <- A0 + A1*sin(gamma*pi/180) + A2*(sin(gamma*pi/180))^2
  
  #calculation of direct, diffuse and global (clear sky) radiation
  Direct <- I0 * eps * sin(gamma*pi/180) * exp(-0.8662*TL*m*delta)
  Diffuse <- I0 * eps * Transmission * Fang
  
  CSR <- Direct + Diffuse
  } else if (typeCS == 2){
    a1 <- (1092 + 1151)/2
    a2 <- (1.238 + 1.312)/2
    CSR <- a1 * mu^a2
  } else if (typeCS == 3){
    CSR <-  1.333 + 345.203*mu + 2043.370*mu^2 - 2199.560*mu^3 + 860.286*mu^4
  }
  CSR <- max(CSR,0)
  return(CSR)  
}

zenithAngle <- function(ToD, date, lat, lon, timezone, typeZenith){
  ###########################
  #VARIABLES
  ###########################
  date2 <- paste0(substr(date,1,4), "-", substr(date,5,6), "-", substr(date,7,8))
  DoY <- as.numeric(strftime(date2, format = "%j"))
  lat <- lat*pi/180
  tempTOD <- ToD - c(0:59)/60

  ###########################
  #MODEL
  ###########################
  ZAdeg <- array(0,c(length(tempTOD)))
  for (t in 1:length(tempTOD)){
  if (typeZenith == 1 ){
    #number of days since 0 UTC, 01-01-2000
    Ndays <- (as.numeric(substr(date,3,4)))*365 + floor((as.numeric(substr(date,3,4))-1)/4) + DoY
    #Mean longitude (in degrees)
    L <- ((280.46 + 0.9856474*Ndays) %% 360)
    #Mean anomaly (in radians)
    g <- ((357.528 + 0.9856003*Ndays) %% 360)*pi/180
    #Ecliptic longitude (in radians)
    lambda <- ((L + 1.915*sin(g) + 0.02*sin(2*g)) %% 360)*pi/180
    #Obliquity of ecliptic (in radians)
    eps <- (23.439 - 0.0000004*Ndays)*pi/180
    #Right ascension (in radians)
    num <- cos(eps)*sin(lambda)
    den <- cos(lambda)
    RA <- atan(num/den)
    if (den < 0){
      RA <- RA + pi
    }
    if (den >= 0 & num < 0){
      RA <- RA + 2*pi
    }
    #Declination angle (in radians)
    DA <- asin(sin(eps)*sin(lambda))
    #Greenwich mean sidereal time
    gmst <- (6.697375 + 0.0657098242*Ndays + tempTOD[t]) %% 24
    #Local mean sidereal time
    lmst <- ((gmst + lon/15) %% 24)*15*pi/180
    #Hour angle (in radians)
    HA <- ((lmst - RA + pi) %% (2*pi)) - pi
    #solar elevation angle (in radians)
    EA <- asin(sin(lat)*sin(DA) + cos(lat)*cos(DA)*cos(HA))
    #solar zenith angle (in degrees)
    ZAdeg[t] <- 90 - EA*180/pi
    #solar azimuth angle (clockwise from north) (in radians)
    #AA <- asin((-cos(DA)*sin(HA))/(cos(EA)))
  } else if (typeZenith == 2){
    #fractional year
    g <- 2*pi/365*(DoY - 1 + (tempTOD[t] - 12)/24)
    #equation of time (in minutes)
    ET <- 229.18*(0.000075 + 0.00187*cos(g) - 0.0321*sin(g) - 0.0146*cos(2*g) - 0.041*sin(2*g))
    #solar declination angle (in radians)
    DA <- 0.0069 - 0.4*cos(g) + 0.07*sin(g) - 0.00676*cos(2*g) + 0.0009*sin(2*g) - 0.0027*cos(3*g) + 0.00148*sin(3*g)
    #true solar time (in minutes)
    tst <- tempTOD[t]*60 + ET - 4*lon + 60*timezone
    #solar hour angle (in radians)
    HA <- (tst/4 - 180)*2*pi/360
    #solar zenith angle (in radians)
    ZA <- acos(sin(lat)*sin(DA) + cos(lat)*cos(DA)*cos(HA))
    #solar zenith angle (in degrees)
    ZAdeg[t] <- ZA*180/pi
    #solar azimuth angle (clockwise from north) (in radians)
    #AA <- acos((sin(lat)*cos(ZA) - sin(DA))/(cos(lat)*sin(ZA)))
  }
  }
  return(mean(ZAdeg))
}

#Reading in zenith angles and clear sky radiations
reading_ZA_CSR <- function(){
  stationPlaces <- c(1:30)
  init_dates <- as.numeric(sub("-","",sub("-","",seq(as.Date('2016-04-01'),as.Date('2018-08-29'),by = "days"))))
  CSR <- array(0, c(length(stationPlaces)*length(init_dates)*24,4))
  ZA <- array(0, c(length(stationPlaces)*length(init_dates)*24,4))
  
  for (place in 1:length(stationPlaces)){
    for (date in 1:length(init_dates)){
      for (time in 1:24){
        CSR[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,1] <- stationData[[2]][place]
        CSR[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,2] <- init_dates[date]
        CSR[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,3] <- time
        CSR[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,4] <- ClearskyRadiation(time,init_dates[date],
                                                                                             stationData[[4]][place],stationData[[3]][place],0, 1, 1)
        
        ZA[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,1] <- stationData[[2]][place]
        ZA[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,2] <- init_dates[date]
        ZA[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,3] <- time
        ZA[length(init_dates)*24*(place - 1) + 24*(date - 1) + time,4] <- zenithAngle(time,init_dates[date],
                                                                                             stationData[[4]][place],stationData[[3]][place],0,1)
      }
    }
  }
  CSR <- data.frame(CSR)
  ZA <- data.frame(ZA)
  colnames(CSR) <- c("Station", "Date", "Time", "CSRadiation")
  colnames(ZA) <- c("Station", "Date", "Time", "ZenithAngles")
  saveRDS(CSR, file = "/usr/people/bakker/kilianbakker/Data/clearsky_data.rds")
  saveRDS(ZA, file = "/usr/people/bakker/kilianbakker/Data/zenithangles_data.rds")
}

ME <- function(var1, var2){
  MEvalue <- mean(var1 - var2, na.rm = TRUE)
  MEperc1 <- MEvalue/mean(var1, na.rm = TRUE)*100
  MEperc2 <- MEvalue/mean(var2, na.rm = TRUE)*100
  return(c(MEvalue, MEperc1, MEperc2))
}

MAE <- function(var1, var2){
  MAEvalue <- mean(abs(var1 - var2), na.rm = TRUE)
  MAEperc1 <- MAEvalue/mean(var1, na.rm = TRUE)*100
  MAEperc2 <- MAEvalue/mean(var2, na.rm = TRUE)*100
  return(c(MAEvalue, MAEperc1, MAEperc2))
}

RMSE <- function(var1, var2){
  n <- length(var1)
  RMSEvalue <- (sum((var1 - var2)^2, na.rm = TRUE)/(n-1))^(1/2)
  RMSEperc1 <-  RMSEvalue/mean(var1, na.rm = TRUE)*100
  RMSEperc2 <-  RMSEvalue/mean(var2, na.rm = TRUE)*100
  return(c(RMSEvalue, RMSEperc1, RMSEperc2))
}

CRPS <- function(var1, var2){
  tempCRPSvalue <- crps_ensemble(as.matrix(var2), as.numeric(var1))
  CRPSvalue <- sum(tempCRPSvalue)/length(var1)
  CRPSperc1 <-  CRPSvalue/mean(var1, na.rm = TRUE)*100
  return(c(CRPSvalue, CRPSperc1))
}

HeidkeScores <- function(observations, forecasts, DiscretizeWidth){
  observations <- floor(observations/DiscretizeWidth)*DiscretizeWidth
  forecasts <- floor(forecasts/DiscretizeWidth)*DiscretizeWidth
  
  start <- round(min(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  end <- round(max(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  
  HS <- 0
  for (i in (start/DiscretizeWidth):(end/DiscretizeWidth)){
    HS <- HS + length(observations[observations == i*DiscretizeWidth & forecasts == i*DiscretizeWidth])
  }
  HS <- HS/length(observations)
  return(HS)
}

PeirceScores <- function(observations, forecasts, DiscretizeWidth){
  observations <- floor(observations/DiscretizeWidth)*DiscretizeWidth
  forecasts <- floor(forecasts/DiscretizeWidth)*DiscretizeWidth

  start <- round(min(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  end <- round(max(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  
  PS <- 0
  for (i in (start/DiscretizeWidth):(end/DiscretizeWidth)){
    PS <- PS + length(observations[observations == i*DiscretizeWidth & forecasts == i*DiscretizeWidth])
  }
  PS <- PS/length(observations)
  return(PS)
}

GandinMurphySkillScores <- function(observations, forecasts, DiscretizeWidth){
  observations_clim <- floor(mean(observations)/DiscretizeWidth)*DiscretizeWidth
  observations <- floor(observations/DiscretizeWidth)*DiscretizeWidth
  forecasts <- floor(forecasts/DiscretizeWidth)*DiscretizeWidth

  start <- round(min(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  end <- round(max(c(observations, forecasts)), digits = log(1/DiscretizeWidth)/log(10))
  
  D <- c((start/DiscretizeWidth):(end/DiscretizeWidth))
  for (i in 1:((end-start)/DiscretizeWidth+1)){
    tempNumber <- 0
    for (k in 1:i){
      tempNumber <- tempNumber + length(observations[observations_clim == (k-1)*DiscretizeWidth-start])
    }
    tempNumber <- tempNumber/length(observations)
    D[i] <- (1-tempNumber)/tempNumber
  }
  S <- array(0, c(((end-start)/DiscretizeWidth+1),((end-start)/DiscretizeWidth+1)))
  for (i in 1:((end-start)/DiscretizeWidth+1)){
    for (j in 1:((end-start)/DiscretizeWidth+1)){
      if (i == j){
        S[i,j] <- (sum(1/D[1:(j-1)]) + sum(D[j:((end-start)/DiscretizeWidth)]))/((end-start)/DiscretizeWidth)
      } else if (i < j){
        S[i,j] <- (sum(1/D[1:(i-1)]) + sum(D[j:((end-start)/DiscretizeWidth)]) - (j-i))/((end-start)/DiscretizeWidth)
        S[j,i] <- S[i,j] 
      }
    }
  }
  
  GM_SS_for <- 0
  for (i in (start/DiscretizeWidth):(end/DiscretizeWidth)){
    for (j in (start/DiscretizeWidth):(end/DiscretizeWidth)){
      GM_SS_for <- GM_SS_for + length(observations[observations == i*DiscretizeWidth & forecasts == j*DiscretizeWidth])*S[i+1-start/DiscretizeWidth,j+1-start/DiscretizeWidth]
    }
  }
  GM_SS <- GM_SS_for/length(observations)
  return(GM_SS)
}

RPS <- function(observations, forecasts, DiscretizeWidth, quants){
  #Input forecasts are in (discretized) CDF form 
  observations <- floor(observations/DiscretizeWidth)*DiscretizeWidth
  forecasts <- floor(forecasts/DiscretizeWidth)*DiscretizeWidth

  start <- round(min(observations), digits = log(1/DiscretizeWidth)/log(10))
  end <- round(max(observations), digits = log(1/DiscretizeWidth)/log(10))
  DiscretizationInterval <- seq(start,end,DiscretizeWidth)
  
  temp_obs <- array(0, c(length(observations),length(DiscretizationInterval)))
  for (f in 1:length(DiscretizationInterval)){
    for (i in 1:length(observations)){
      if (observations[i] > DiscretizationInterval[f]){
        temp_obs[i,f] <- 1
      }
    }
  }
  
  temp_RPS_for <- array(0,c(length(observations), length(DiscretizationInterval)))
  for (i in 1:length(observations)){
    for (f in 1:length(DiscretizationInterval)){
      temp_RPS_for[i,f] <- (forecasts[i,f] - temp_obs[i,f])^2
    }
  }
  RPS <- sum(temp_RPS_for)/length(observations)
  return(RPS)
}

POVscores <- function(observations, forecasts, cost, loss){
  PEV <- c(1:length(observations))
  for (i in 1:length(observations)){
    if (observations[i] == 1){
      PEVscore[i] <- Cost*forecasts[i,1]
    } else if (observations[i] == 0){
      PEVscore[i] <- Cost*forecasts[i,1] + Loss*forecasts[i,2]
    }
  }
  PEVscore <- sum(PEV)/length(observations)
  return(PEVscore)
}

QEVscores <- function(observations, forecasts, quants){
  check_loss_function <- function(quantile, observation, forecast){
    if (observation - forecast >= 0){
      value <- quantile*abs(observation - forecast)
    } else if (observation - forecast < 0){
      value <- (1-quantile)*abs(observation - forecast)
    }
  }
  QEV <- array(0, c(length(observations), length(quants)))
  for (i in 1:length(observations)){
    for (j in 1:length(quants)){
      QEV[i,j] <- check_loss_function(quants[j], observations[i], forecasts[i,j])
    }
  }
  QEVscore <- sum(QEV)/(length(observations)*length(quants))
  return(QEVscore)
}

FitPlot <- function(plotData, xpos1, ypos1, delta_ypos1){
  fit <- lm(plotData[[1]] ~ plotData[[2]], data = plotData)
  
  ggplot(data = plotData, mapping = aes(x = plotData[[2]], y = plotData[[1]])) +
    geom_point() +
    geom_smooth(method='lm') +
    annotate("text", x = xpos1, y = ypos1, label = paste0(names(plotData)[[1]], " = ", round(coefficients(fit)[[1]], digits = 2), " + ", 
                                                          round(coefficients(fit)[[2]], digits = 2), " * ", names(plotData)[[2]])) +
    annotate("text", x = xpos1, y = ypos1 - delta_ypos1, label = paste0("R^2 = ", round(summary(fit)[[8]], digits = 4))) +
    xlab(names(plotData)[[2]]) + ylab(names(plotData)[[1]])
}

RadiationPlot <- function(plotData, xpos1, ypos1, delta_ypos1){
  biases <- MAE(plotData[[1]], plotData[[2]])
  scores <- RMSE(plotData[[1]],plotData[[2]])
  fit <- lm(plotData[[1]] ~ plotData[[2]], data = plotData)
  
  ggplot(data = plotData, mapping = aes(x = plotData[[2]], y = plotData[[1]])) +
    geom_point() +
    geom_smooth(method='lm') +
    annotate("text", x = xpos1, y = ypos1, label = paste0("For = ", round(coefficients(fit)[[1]], digits = 2), " + ", 
                                                          round(coefficients(fit)[[2]], digits = 2), " * Obs")) +
    annotate("text", x = xpos1, y = ypos1 - delta_ypos1, label = paste0("R^2 = ", round(summary(fit)[[8]], digits = 4))) + 
    annotate("text", x = xpos1, y = ypos1 - 3*delta_ypos1, label = paste0("meanObs = ", round(mean(plotData[[1]],na.rm = TRUE), digits = 2))) +
    annotate("text", x = xpos1, y = ypos1 - 4*delta_ypos1, label = paste0("meanFor = ", round(mean(plotData[[2]],na.rm = TRUE), digits = 2))) +
    annotate("text", x = xpos1, y = ypos1 - 5*delta_ypos1, label = paste0("bias = ", round(biases[1], digits = 2), " (", round(biases[2],digits = 2), "%) ")) +
    annotate("text", x = xpos1, y = ypos1 - 6*delta_ypos1, 
             label = paste0("RMSE = ", round(scores[1], digits = 2), " (", round(scores[2],digits = 2), "%)")) + 
    xlab(names(plotData)[[2]]) + ylab(names(plotData)[[1]]) + ggtitle(paste0(names(plotData)[[2]], " VS observations"))
}

importing_forecasts <- function(init_dates_CS, leadTimes){
  tempforecasts <- c()
  for (k in 1:length(init_dates_CS)){
    yearNumber <- substring(as.character(init_dates_CS[k]),1,4)
    monthNumber <- substring(as.character(init_dates_CS[k]),5,6)
    dayNumber <- substring(as.character(init_dates_CS[k]),7,8)
    filename <- paste0(yearNumber,"/",monthNumber,"/",dayNumber,"/")
    url <- paste0("ftp://csr_knmi:C5r-Knm1!@bbc.knmi.nl/MEMBERS/knmi/bsrn/harmonie/",filename)
    filenames <- getURL(url, ftp.use.epsv = TRUE ,dirlistonly = TRUE, verbose = TRUE) 
    filenames <- paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
    forecastData <-read.table(filenames[1], col.names=c("xgridbox","ygridbox","Direct","Global", "CC"), fill=TRUE)
    tempforecasts2 <- diff(filter(forecastData, xgridbox == 134, ygridbox == 130)[[4]][c(leadTimes[k], leadTimes[k] + 1)])/3600
    tempforecasts <- c(tempforecasts,tempforecasts2)
  }
  return(tempforecasts)
}

coordinatesLargeGrid <- function(stationData, stationPlaces, xDomainSizes, yDomainSizes){
  xcoordinates <- c()
  ycoordinates <- c()
  tmprds         <- readRDS(file = "/nobackup/users/bakker/Data/temperaturevariables/20160401.rds")
  tempxcoordinates   <- unique(tmprds[[2]])
  tempycoordinates   <- unique(tmprds[[3]])
  for (i in 1:length(stationPlaces)){
    position      <- which(min(abs(tempxcoordinates - stationData[[3]][[stationPlaces[i]]])) == abs(tempxcoordinates - stationData[[3]][[stationPlaces[i]]]))[1]
    xcoordinates <- c(xcoordinates,tempxcoordinates[position + xDomainSizes])
    position2     <- which(min(abs(tempycoordinates - stationData[[4]][[stationPlaces[i]]])) == abs(tempycoordinates - stationData[[4]][[stationPlaces[i]]]))[1]
    ycoordinates    <- c(ycoordinates,tempycoordinates[position2 + yDomainSizes])
  }
  return(list(xcoordinates, ycoordinates))
}

coordinatesSmallGrid <- function(stationData, stationPlaces){
  xcoordinates <- c()
  ycoordinates <- c()
  tmpraster      <- crop(raster("/net/pc150398/nobackup_1/users/meirink/msgdev/CAMS/CAMS_20170808.nc", varname = "aot500"), extent(3.4, 7.1, 50.8, 53.5))
  tempxcoordinates   <- unique(coordinates(tmpraster)[,1])
  tempycoordinates   <- unique(coordinates(tmpraster)[,2])
    for (i in 1:length(stationPlaces)){
    position      <- which(min(abs(tempxcoordinates - stationData[[3]][[stationPlaces[i]]])) == abs(tempxcoordinates - stationData[[3]][[stationPlaces[i]]]))[1]
    xcoordinates   <- c(xcoordinates,tempxcoordinates[position])
    position2      <- which(min(abs(tempycoordinates - stationData[[4]][[stationPlaces[i]]])) == abs(tempycoordinates - stationData[[4]][[stationPlaces[i]]]))[1]
    ycoordinates   <- c(ycoordinates,tempycoordinates[position2])
  }
  return(list(xcoordinates, ycoordinates))
}

importing_variables <- function(stationData, stationPlaces, typeCS, typeZenith, init_dates, Includedvariables, 
                                variableNamesCorr, leadTimes, leadTimeSizes, xDomainSizes, yDomainSizes, heights2, filename){
AllStationsvariableData <- array(0, c(length(init_dates)*length(stationPlaces), length(variableNamesCorr), length(leadTimes)))
for (place in 1:length(stationPlaces)){
  xcoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], xDomainSizes, yDomainSizes)[[1]]
  ycoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], xDomainSizes, yDomainSizes)[[2]]
  xcoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[1]]
  ycoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[2]]
  length1 <- length(xcoordinate)*length(ycoordinate)
  tempvariableData <- c()
for (j in 1:length(init_dates)){
  temperatureData <- readRDS(file = paste0("/nobackup/users/bakker/Data/temperaturevariables/",init_dates[j],".rds"))
  temperatureData <- filter(temperatureData, xcoor %in% xcoordinate, ycoor %in% ycoordinate)

  humidityData <- readRDS(file = paste0("/nobackup/users/bakker/Data/humidityvariables/",init_dates[j],".rds"))
  humidityData <- filter(humidityData, xcoor %in% xcoordinate, ycoor %in% ycoordinate)

  radiationData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",init_dates[j],".rds"))
  radiationData <- filter(radiationData, xcoor %in% xcoordinate, ycoor %in% ycoordinate)
  tempRadiationData <- (radiationData[(length1+1):length(radiationData[[1]]), 4:8] - 
                       radiationData[1:(length(radiationData[[1]]) - length1), 4:8])/3600
  radiationData <- rbind(data.frame(Global = array(0,c(length1)), Direct_SURF = array(0,c(length1)), Direct_TOA = array(0,c(length1)), 
                                        NCS_SURF = array(0,c(length1)), NCS_TOA = array(0,c(length1))), tempRadiationData)
    
  cloudData <- readRDS(file = paste0("/nobackup/users/bakker/Data/cloudvariables/",init_dates[j],".rds"))
  cloudData <- filter(cloudData, xcoor %in% xcoordinate, ycoor %in% ycoordinate)
  tempcloudData <- c(array(0,c(length1)), cloudData[(length1+1):length(cloudData[[1]]), 16] - cloudData[1:(length(cloudData[[1]]) - length1), 16])
  cloudData[,16] <- tempcloudData

  particlesData <- readRDS(file = paste0("/nobackup/users/bakker/Data/particlesvariables/",init_dates[j],".rds"))
  particlesData <- filter(particlesData, xcoor %in% xcoordinate2, ycoor %in% ycoordinate2)
  tempTimes <- unique(particlesData[[1]])
  tempTimes2 <- as.numeric(as.character(temperatureData[[1]]))
  closestleadTimes <- c(1:length(tempTimes2))
  for (k in 1:length(temperatureData[[1]])){
    closestleadTimes[k] <- which(min(abs(tempTimes2[k] - tempTimes)) == abs(tempTimes2[k] - tempTimes))
  }
  tempParticlesData <- data.frame(AOD = particlesData[closestleadTimes,4], ANG = particlesData[closestleadTimes,5], OZ = particlesData[closestleadTimes,6])

  tempData <- data.frame(Date = init_dates[j],temperatureData, humidityData[,4:16], radiationData[,1:5], cloudData[,4:18], tempParticlesData[,1:3])
  tempvariableData <- rbind(tempvariableData,tempData)
}

  variableData <- array(0, c(length(init_dates), length(variableNamesCorr), length(leadTimes)))
  for (time in 1:length(leadTimes)){
    leadTime <- leadTimes[time]
    
    CSR <- c(1:length(init_dates))
    if (leadTime <= 24){
      temp_init_dates <- init_dates
    } else {
      temp_init_dates <- c(init_dates[-1],20180317)
    }
    for (k in 1:length(init_dates)){
      CSR[k] <- ClearskyRadiation((leadTime %% 24),temp_init_dates[k],stationData[[4]][place],stationData[[3]][place],0, typeCS, typeZenith)
    }

    tempvariableData2 <- filter(tempvariableData, as.numeric(as.character(lt)) %in% c(leadTime + leadTimeSizes))

  for (p in 1:length(init_dates)){
    for (q in 1:5){
      variableData[p,q, time] <- round(mean(filter(tempvariableData2, Date == init_dates[p])[4 + Includedvariables[q]][[1]]), digits = 4)/CSR[p]
    }
    for (q in 6:length(Includedvariables)){
      variableData[p,q, time] <- round(mean(filter(tempvariableData2, Date == init_dates[p])[4 + Includedvariables[q]][[1]]), digits = 4)
    }
    variableData[p,length(Includedvariables) + 1, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[6:12]))*diff(heights2)[1:7])/sum(diff(heights2)[1:7]), digits = 4)
    variableData[p,length(Includedvariables) + 2, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[13:15]))*diff(heights2)[8:10])/sum(diff(heights2)[8:10]), digits = 4)
    variableData[p,length(Includedvariables) + 3, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[16:18]))*diff(heights2)[11:13])/sum(diff(heights2)[11:13]), digits = 4)
    variableData[p,length(Includedvariables) + 4, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[19:25]))*diff(heights2)[1:7])/sum(diff(heights2)[1:7]), digits = 4)
    variableData[p,length(Includedvariables) + 5, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[26:28]))*diff(heights2)[8:10])/sum(diff(heights2)[8:10]), digits = 4)
    variableData[p,length(Includedvariables) + 6, time] <- round(sum(as.numeric(colMeans(filter(tempvariableData2, Date == init_dates[p])[29:31]))*diff(heights2)[11:13])/sum(diff(heights2)[11:13]), digits = 4)
  }
  }
  AllStationsvariableData[(length(init_dates)*(place - 1) + 1):(length(init_dates)*place), , ] <- variableData
}
saveRDS(AllStationsvariableData, file  = paste0("/usr/people/bakker/kilianbakker/Data/",filename,".rds"))
}

importingclearskyhours <- function(Globalmin, Globalmax, Difmax, sdIntervalWidth, Stdevmax){
init_dates <- as.numeric(sub("-","",sub("-","",seq(as.Date('2016-04-01'),as.Date('2018-07-16'),by = "days"))))
ClearSkyHours <- array(0,c(24,length(init_dates)))

for (p in 1:length(init_dates)){
  for (q in 1:24){
    dateNumber <- init_dates[p]
    leadTime <- q
    
    observationData <- read.csv(paste0("/nobackup/users/bakker/Data/observations_min/",dateNumber,".csv"))
    
    Adjustment <- 0
    observations_minute <- observationData[[7]][(60*(leadTime-1)+1-Adjustment):(60*leadTime - Adjustment)]
    diffuse_minute <- observationData[[3]][(60*(leadTime-1)+1-Adjustment):(60*leadTime - Adjustment)]
    zenith_minute <- observationData[[2]][(60*(leadTime-1)+1-Adjustment):(60*leadTime - Adjustment)]
    zenith_noon <- max(observationData[[2]])
    
    TestCounter <- 0
    if (sum(as.numeric(is.na(observations_minute))) <= 10 & sum(as.numeric(is.na(diffuse_minute))) <= 10){
      MissingValues1 <- which(is.na(observations_minute))
      for (k in 1:length(MissingValues1)){
        observations_minute[MissingValues1[k]] <- mean(c(observations_minute[MissingValues1[k]-2],observations_minute[MissingValues1[k]-1],
                                                         observations_minute[MissingValues1[k]+1],observations_minute[MissingValues1[k]+2]), na.rm = T)
      }
      MissingValues2 <- which(is.na(diffuse_minute))
      for (k in 1:length(MissingValues1)){
        diffuse_minute[MissingValues2[k]] <- mean(c(diffuse_minute[MissingValues2[k]-2],diffuse_minute[MissingValues2[k]-1],
                                                    diffuse_minute[MissingValues2[k]+1],diffuse_minute[MissingValues2[k]+2]), na.rm = T)
      }
      if (sum(as.numeric(zenith_minute > 0)) == 60){
        #Test 1
        if (sum(as.numeric(observations_minute/zenith_minute >= Globalmin)) == 60){
          TestCounter <- TestCounter + 0.5
        }
        if (sum(as.numeric(observations_minute/zenith_minute <= Globalmax)) == 60){
          TestCounter <- TestCounter + 0.5
        }
        
        #Test 2
        if (sum(as.numeric(diffuse_minute/zenith_minute^0.5 <= Difmax)) == 60){
          TestCounter <- TestCounter + 1  
        }
        
        #Test 3
        observations_TOA <- 1365*zenith_minute
        DeltaT <- 1
        C <- 20
        R <- 10
        tempsum1 <- c(1:(60-2*DeltaT))
        tempsum2 <- c(1:(60-2*DeltaT))
        for (i in (1 + DeltaT):(60 - DeltaT)){
          Actual <- abs(observations_minute[i - DeltaT] - observations_minute[i + DeltaT])/(2*DeltaT)
          Max <- abs(observations_TOA[i - DeltaT] - observations_TOA[i + DeltaT])/(2*DeltaT) + C*zenith_minute[i]
          Min <- abs(observations_TOA[i - DeltaT] - observations_TOA[i + DeltaT])/(2*DeltaT) - R*(zenith_noon + 0.1)/zenith_minute[i]
          tempsum1[i - DeltaT] <- as.numeric(Actual <= Max)
          tempsum2[i - DeltaT] <- as.numeric(Actual >= Min)
        }
        
        if (sum(tempsum1) == (60 - 2*DeltaT)){
          TestCounter <- TestCounter + 0.5
        } 
        if (sum(tempsum2) == (60 - 2*DeltaT)){
          TestCounter <- TestCounter + 0.5
        }
        
        #Test 4
        StandardDeviations <- c(1:(60 - sdIntervalWidth + 1))
        for (i in 1:(60-sdIntervalWidth + 1)){
          StandardDeviations[i] <- sd(diffuse_minute[i:(i+sdIntervalWidth - 1)]/(zenith_minute[i:(i+sdIntervalWidth - 1)]*observations_minute[i:(i+sdIntervalWidth - 1)]))
        }
        
        if (sum(as.numeric(StandardDeviations <= Stdevmax)) == (60 - sdIntervalWidth + 1)){
          TestCounter <- TestCounter + 1
        }
      }
    }
    if (TestCounter == 4){
      ClearSkyHours[q,p] <- 1
    }
  }
}

ClearSkyHours <- rbind(ClearSkyHours[,1:length(ClearSkyHours[1,])-1],ClearSkyHours[,2:length(ClearSkyHours[1,])])
ClearSkyHours <- data.frame(ClearSkyHours)
colnames(ClearSkyHours) <- init_dates[-length(init_dates)]
saveRDS(ClearSkyHours, file = "/usr/people/bakker/kilianbakker/Data/clearskyhours.rds")
}

importinglargeErrorhours <- function(Errorbound, stationData, stationPlace, init_dates){
  LargeErrorHours <- array(0,c(48,length(init_dates)-1))
  observationData <-read.table("/usr/people/bakker/kilianbakker/Data/observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
  observationData[[4]] <- observationData[[4]]*(10000/3600)
  
  xcoordinate <- coordinatesLargeGrid(stationData, stationPlace, c(0), c(0))[1]
  ycoordinate <- coordinatesLargeGrid(stationData, stationPlace, c(0), c(0))[2]
  
  for (p in 1:(length(init_dates)-1)){
    observations <- filter(observationData, Station == stationData[[2]][[stationPlace]], Date %in% c(init_dates[p], init_dates[p+1]))[[4]]
    forecastData <- readRDS(file = paste0("/nobackup/users/bakker/Data/radiationvariables/",init_dates[p],".rds"))
    tempforecasts <- filter(forecastData, xcoor == xcoordinate, ycoor == ycoordinate)[[4]]
    forecastsDay1 <- c(0,0,0,c(tempforecasts[1],diff(tempforecasts[1:17]))/3600,0,0,0,0)
    forecastsDay2 <- c(0,0,0,c(tempforecasts[18],diff(tempforecasts[18:34]))/3600,0,0,0,0)
    forecasts <- c(forecastsDay1,forecastsDay2)
    for (q in 1:48){
      if (observations[q] == 0){
        LargeErrorHours[q,p] <- NA
      } else {
        if (abs(forecasts[q] - observations[q])/max(observations[q],forecasts[q]) > Errorbound){
          LargeErrorHours[q,p] <- 1
        }
      }
    }
  }
  
  LargeErrorHours <- data.frame(LargeErrorHours)
  colnames(LargeErrorHours) <- init_dates[-length(init_dates)]
  saveRDS(LargeErrorHours, file = "/usr/people/bakker/kilianbakker/Data/largeErrorhours.rds")
}

predict.multqrsvm.new <- function(model, newdata) {
  requireNamespace("kernlab", quietly = TRUE)
  
  if (class(model)=="multqrsvm"){
    prediction<-list()
    
    if (ncol(newdata) != ncol(model[[1]]$xtrain)) {
      cat("Newdata has different number of columns than xtrain please check consistency!",
          fill = TRUE)
    }
    for (i in 1:length(model)){
      xold <- model[[i]]$xtrain
      alpha <- model[[i]]$alpha
      kern <- model[[i]]$kernel
      b <- model[[i]]$b0
      
      prediction[[i]] <- kernelMult(kern, newdata, xold, alpha) + b
    }
  }
  
  return(prediction)
}

predict.qrsvm.new <- function(model, newdata) {
  requireNamespace("kernlab", quietly = TRUE)
  
  
  
  if (class(model)=="qrsvm"){
    xold <- model$xtrain
    
    if (ncol(newdata) != ncol(xold)) {
      cat("Newdata has different number of columns than xtrain please check consistency!",
          fill = TRUE)
    }
    alpha <- model$alpha
    kern <- model$kernel
    b <- model$b0
    prediction <- kernelMult(kern, newdata, xold, alpha) + b
  }
  return(prediction)
}

calculating_scores <- function(ContMethod, ProbMethod, observations, forecasts, DiscretizeWidth){
  if (ProbMethod == 1){
    observations <- observations[is.na(forecasts) == F]
    forecasts <- na.omit(forecasts)
  } else if (ProbMethod == 2){
    for (l in 1:length(observations)){
      for (q in 1:length(quants)){
        if (is.na(forecasts[l,q])){
          forecasts[l,] <- NA
        }
      }
    }
    observations <- observations[is.na(forecasts[,1]) == F]
    forecasts <- na.omit(forecasts)
  }
  if (length(observations) > 1){
    if (ContMethod == 1 & ProbMethod == 1){
      Score <- HeidkeScores(observations, forecasts, DiscretizeWidth)
    } else if (ContMethod == 2 & ProbMethod == 1){
      Score <- RMSE(observations, forecasts)[1]
    } else if (ContMethod == 1 & ProbMethod == 2){
      Score <- RPS(observations, forecasts, DiscretizeWidth)
    } else if (ContMethod == 2 & ProbMethod == 2){
      Score <- CRPS(observations, forecasts)[1]
    }
  } else {
    Score <- NA
  }
  return(Score)
}

saving_dailyclimatologies <- function(ProbMethod, quants, stationData, stationPlaces){
  observationData2 <-read.table("/usr/people/bakker/kilianbakker/Data/extended_observation_data.txt", sep=",", col.names=c("Station","Date","Time","Radiation"), fill=TRUE)
  observationData2[[4]] <- observationData2[[4]]*(10000/3600)
  Times <- c(1:24)

  if (ProbMethod == 1){
    dailyclimatologies <- array(NA, c(365,length(stationPlaces), length(Times)))
  } else if (ProbMethod == 2){
    dailyclimatologies <- array(NA, c(365,length(stationPlaces), length(Times), length(quants)))
  }
  for (time in 1:length(Times)){
  for (place in 1:length(stationPlaces)){
  stationPlace <- stationPlaces[place]
  tempobservationData2 <- filter(observationData2, Station %in% stationData[[2]][stationPlace], Time == (leadTime %% 24), !(Date %in% c(20080229,20120229,20160229)))

  yearNumbers <- c(2008:2017)
  observations <- array(0,c(365,length(yearNumbers)))
  for (i in 1:10){
    observations[,i] <- filter(tempobservationData2, floor(Date/10000) == yearNumbers[i])[[4]]
  }
  if (ProbMethod == 1){
    dailyclimatologies[,place,time] <- rowMeans(observations)
  } else if (ProbMethod == 2){
    for (d in 1:365){
      dailyclimatologies[d,place,time,] <- unname(quantile(observations[d,], quants))
    }
  }
  }
  }
  saveRDS(dailyclimatologies, file = paste0("/usr/people/bakker/kilianbakker/Data/dailyclimatologies_",ProbMethod,".rds"))
}

calculating_dailyclimatologies <- function(ContMethod, ProbMethod, leadTime, observations, init_dates, quants, stationPlace, DiscretizeWidth){
  tempdailyclimatologies <- readRDS(file = paste0("/usr/people/bakker/kilianbakker/Data/dailyclimatologies_",ProbMethod,".rds"))
  if (ProbMethod == 1){
    dailyclimatologies <- array(NA, c(length(init_dates)))
  } else if (ProbMethod == 2){
    dailyclimatologies <- array(NA, c(length(init_dates), length(quants)))
  }
  
    if (leadTime > 24){
      for (d in 1:length(init_dates)){
        tempDate <- paste0(substr(init_dates[d],1,4), "-", substr(init_dates[d],5,6), "-", substr(init_dates[d],7,8))
        init_dates[d] <- as.numeric(gsub("-","",as.Date(tempDate) + 1))
      }
    }
    leadTime <- leadTime %% 24
    
    for (k in 1:length(init_dates)){
      tempdate <- paste0(substr(init_dates[k],1,4), "-", substr(init_dates[k],5,6), "-", substr(init_dates[k],7,8))
      DoY <- as.numeric(strftime(tempdate, format = "%j"))
      if (as.numeric(substr(init_dates[k],1,4)) == 2016){
        DoY <- DoY - 1
      }
      if (ProbMethod == 1){
        dailyclimatologies[k] <- tempdailyclimatologies[DoY,stationPlace,leadTime]
      } else if (ProbMethod == 2){
        dailyclimatologies[k,] <- tempdailyclimatologies[DoY,stationPlace,,leadTime]
      }
    }

      if (ContMethod == 1){
    dailyclimatologies <- round(dailyclimatologies, digits = log(1/DiscretizeWidth)/log(10))
  }
  
  clim_scores <- calculating_scores(ContMethod, ProbMethod, observations, dailyclimatologies, DiscretizeWidth)
  return(clim_scores)
}

calculating_sampleclimatologies <- function(ContMethod, ProbMethod, observations, CSR, init_dates, quants, DiscretizeWidth){
    init_dates_temp <- init_dates[which(CSR > 20)]
    if (length(init_dates_temp) > 0){
    tempobservations <- observations[which(init_dates %in% init_dates_temp)]
    if (ProbMethod == 1){
      sampleclimatologies <- array(mean(tempobservations), c(length(init_dates)))
    } else if (ProbMethod == 2){
      sampleclimatologies <- array(unname(quantile(tempobservations, quants)), c(length(init_dates), length(quants)))
    }
    }
  if (ContMethod == 1){
    sampleclimatologies <- round(sampleclimatologies, digits = log(1/DiscretizeWidth)/log(10))
  }
  clim_scores <- calculating_scores(ContMethod, ProbMethod, observations, sampleclimatologies, DiscretizeWidth)
  return(clim_scores)
}

calculating_distances <- function(stationxcoor, stationycoor, Type){
Distances <- array(0, c(length(stationxcoor)))
if (Type == "Coast"){
  worldmap <- map("world", interior = F, xlim = c(3.4,7.1), ylim = c(50.8,53.5), plot = F)
  
  for (i in 1:length(worldmap$x)){
    if (is.na(worldmap$x[i]) == F & is.na(worldmap$y[i]) == F){
      if (worldmap$x[i] < 3.4 | worldmap$x[i] > 7.1 | worldmap$y[i] < 50.8 | worldmap$y[i] > 53.5){
        worldmap$x[i] <- NA
        worldmap$y[i] <- NA
      } else if (worldmap$x[i] > 4.9 & worldmap$x[i] < 6 & worldmap$y[i] > 52.2 & worldmap$y[i] < 52.9){
        worldmap$x[i] <- NA
        worldmap$y[i] <- NA
      }
    }
  }
  
  tempx <- na.omit(worldmap$x)
  tempy <- na.omit(worldmap$y)
  
  for (j in 1:length(stationxcoor)){
    tempNumber <- which(min(sqrt((stationxcoor[j] - tempx)^2 + (stationycoor[j] - tempy)^2)) == 
                                                          sqrt((stationxcoor[j] - tempx)^2 + (stationycoor[j] - tempy)^2))[1]
    Distances[j] <- sqrt((stationxcoor[j] - tempx[tempNumber])^2 + (stationycoor[j] - tempy[tempNumber])^2)
  }
} else if (Type == "Water"){
  worldmap <- map("world", interior = F, xlim = c(3.4,7.1), ylim = c(50.8,53.5), plot = F)
  
  for (i in 1:length(worldmap$x)){
    if (is.na(worldmap$x[i]) == F & is.na(worldmap$y[i]) == F){
      if (worldmap$x[i] < 3.4 | worldmap$x[i] > 7.1 | worldmap$y[i] < 50.8 | worldmap$y[i] > 53.5){
        worldmap$x[i] <- NA
        worldmap$y[i] <- NA
      }
    }
  }
  
  tempx <- na.omit(worldmap$x)
  tempy <- na.omit(worldmap$y)
  
  for (j in 1:length(stationxcoor)){
    tempNumber <- which(min(sqrt((stationxcoor[j] - tempx)^2 + (stationycoor[j] - tempy)^2)) == 
                          sqrt((stationxcoor[j] - tempx)^2 + (stationycoor[j] - tempy)^2))[1]
    Distances[j] <- sqrt((stationxcoor[j] - tempx[tempNumber])^2 + (stationycoor[j] - tempy[tempNumber])^2)
  }
} else if (Type == "Inland"){
  for (j in 1:length(stationxcoor)){
    ThreeCountriesPoint <- c(6.017, 50.75)
    Distances[j] <- sqrt((stationxcoor[j] - ThreeCountriesPoint[1])^2 + (stationycoor[j] - ThreeCountriesPoint[2])^2)
  }
}
return(Distances)
}