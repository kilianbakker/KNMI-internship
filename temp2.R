tmprds <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_9x9gridbox_3lt.rds"))


tempTimes <- seq(0,48,3)
  
  for (date in 1:length(init_dates)){
    particlesData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/particlesvariables/",init_dates[date],".rds"))
    
    for (place in 1:length(stationPlaces)){
      xcoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[1]]
      ycoordinate2 <- coordinatesSmallGrid(stationData,stationPlaces[place])[[2]]

      for (time in 1:length(leadTimes)){
        leadTime <- leadTimes[time]
        closestLeadTime <- tempTimes[which(min(abs(tempTimes - leadTime)) == abs(tempTimes - leadTime))]
        particles <- as.numeric(unname(filter(particlesData, xcoor %in% xcoordinate2, ycoor %in% ycoordinate2, as.numeric(as.character(lt)) == closestLeadTime)[4:6]))
        tmprds[date,place,time,c(19:21)] <- particles
      }
    }
  }

  tmprds[,,c(1,16),c(1:5,6:18,22:27)] <- tmprds[,,c(1,16),c(1:5,6:18,22:27)]*1.5
  saveRDS(tmprds, file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_9x9gridbox_3lt.rds"))
  