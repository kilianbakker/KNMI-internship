tmprds <- readRDS(file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_9x9gridbox_3lt.rds"))
tmprds[,,c(1,16),c(1:5,6:18,22:27)] <- tmprds[,,c(1,16),c(1:5,6:18,22:27)]/1.5
tmprds[,,c(1,16),c(1:5)] <- tmprds[,,c(1,16),c(1:5)]/1.5


leadTimeSizes <- c(-1,0,1)
DomainSizes <- c(-4,-3,-2,-1,0,1,2,3,4)
  for (date in 1:length(init_dates)){
    print(date)
    Sys.sleep(0.01)
    cloudData <- readRDS(file = paste0("/nobackup/users/bakker/Data2/cloudvariables/",init_dates[date],".rds"))

    for (place in 1:length(stationPlaces)){
      xcoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], DomainSizes, DomainSizes, cloudData)[[1]]
      ycoordinate <- coordinatesLargeGrid(stationData,stationPlaces[place], DomainSizes, DomainSizes, cloudData)[[2]]

      for (time in 1:length(leadTimes)){
        leadTime <- leadTimes[time]

        leadTimesRegion <- leadTime + leadTimeSizes
        if (leadTime %in% c(5,29)){
          leadTimesRegion <- leadTime + leadTimeSizes[-1]
        }
        
        clouds <- colMeans(filter(cloudData, xcoor %in% xcoordinate, ycoor %in% ycoordinate, as.numeric(as.character(lt)) %in% (leadTimesRegion))[8:11])
        tmprds[date,place,time,c(10:13)] <- as.numeric(unname(clouds))
      }
    }
  }
  saveRDS(tmprds, file  = paste0("/usr/people/bakker/kilianbakker/Data/variables_9x9gridbox_3lt.rds"))
  