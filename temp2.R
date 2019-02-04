leadTimeSizes <- c(-1,0,1)
DomainSizes <- c(-3,-2,-1,0,1,2,3)
varSize <- "7x7gridbox_3lt"

tmprds <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_",varSize,".rds"))

tmprds[,,c(1,16),c(1:5)] <- tmprds[,,c(1,16),c(1:5)]/1.5

tempTimes <- seq(0,48,3)
  
for (date in 1:length(init_dates)){
  print(date)
  Sys.sleep(0.01)
  temperatureData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/temperaturevariables/",init_dates[date],".rds"))
  humidityData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/humidityvariables/",init_dates[date],".rds"))
  radiationData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/radiationvariables/",init_dates[date],".rds"))
  cloudData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/cloudvariables/",init_dates[date],".rds"))
  particlesData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/particlesvariables/",init_dates[date],".rds"))
  
  for (place in 1:length(stationPlaces)){
    xcoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], DomainSizes, DomainSizes, temperatureData)[[1]]
    ycoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], DomainSizes, DomainSizes, temperatureData)[[2]]
    xcoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[1]]
    ycoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[2]]
    LowIndices <- c(1:7)
    MediumIndices <- c(8:10)
    HighIndices <- c(11:13)
    
    for (time in 1:length(leadTimes)){
      leadTime <- leadTimes[time]
      closestLeadTime <- tempTimes[which(min(abs(tempTimes - leadTime)) == abs(tempTimes - leadTime))]
      particles <- as.numeric(unname(filter(particlesData, xcoor %in% xcoordinate2, ycoor %in% ycoordinate2, as.numeric(as.character(lt)) == closestLeadTime)[4:6]))
      tmprds[date,place,time,c(19:21)] <- particles
      
      
      leadTimesRegion <- leadTime + leadTimeSizes
      if (leadTime %in% c(5,29)){
        leadTimesRegion <- leadTime + leadTimeSizes[-1]
      }
      
      clouds <- colMeans(filter(cloudData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[8:11])
      tmprds[date,place,time,c(10:13)] <- as.numeric(unname(clouds))
      
      if (length(leadTimeSizes) == 1){
        if (time %in% c(1,16)){
        leadTimesRegion <- leadTime

        radiations <- colMeans(filter(radiationData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[4:8])
        previousradiations <- colMeans(filter(radiationData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion - 1))[4:8])
        
        tmprds[date,place,time,1:5] <- as.numeric(unname((radiations - previousradiations)/3600))
        
        clouds <- colMeans(filter(cloudData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[4:16])
        previousRain <- colMeans(filter(cloudData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion - 1))[16])
        clouds[length(clouds)] <- clouds[length(clouds)] - previousRain
        
        tmprds[date,place,time,6:18] <- as.numeric(unname(clouds))
        
        temperatures <- colMeans(filter(temperatureData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[5:17])
        humidities <- colMeans(filter(humidityData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[4:16])
        
        tmprds[date,place,time,22] <- sum(temperatures[LowIndices]*diff(heights2)[LowIndices], na.rm = T)/sum(diff(heights2)[LowIndices])
        tmprds[date,place,time,23] <- sum(temperatures[MediumIndices]*diff(heights2)[MediumIndices], na.rm = T)/sum(diff(heights2)[MediumIndices])
        tmprds[date,place,time,24] <- sum(temperatures[HighIndices]*diff(heights2)[HighIndices], na.rm = T)/sum(diff(heights2)[HighIndices])
        tmprds[date,place,time,25] <- sum(humidities[LowIndices]*diff(heights2)[LowIndices], na.rm = T)/sum(diff(heights2)[LowIndices])
        tmprds[date,place,time,26] <- sum(humidities[MediumIndices]*diff(heights2)[MediumIndices], na.rm = T)/sum(diff(heights2)[MediumIndices])
        tmprds[date,place,time,27] <- sum(humidities[HighIndices]*diff(heights2)[HighIndices], na.rm = T)/sum(diff(heights2)[HighIndices])
        }
      }
    }
  }
}

saveRDS(tmprds, file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_upd_",varSize,".rds"))
