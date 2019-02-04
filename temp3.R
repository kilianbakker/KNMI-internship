#low amount of clouds, high amount of aerosols
dates1 <- c()
for (i in 1:670){
  if (length(which(AllvariableData[i,23,1:15,19] > 0.5)) > 7.5){
    dates1 <- c(dates1,i)  
  }
}
idates1 <- init_dates[sort(unique(dates1))]

dates2 <- c()
for (i in 1:670){
  if (length(which(AllvariableData[i,23,1:15,7] == 0)) > 7.5){
    dates2 <- c(dates2,i)  
  }
}
idates2 <- init_dates[sort(unique(dates2))]
intersect(idates1,idates2)


#Fully overcast, no rain
dates3 <- c()
for (i in 1:670){
  if ((length(which(AllvariableData[i,23,1:15,6] > 0.75)) == 15) && (length(which(AllvariableData[i,23,1:15,7] > 0.75)) == 15)){
    dates3 <- c(dates3,i)  
  }
}
idates3 <- init_dates[sort(unique(dates3))]

dates4 <- c()
for (i in 1:670){
  if (length(which(AllvariableData[i,23,1:15,18] == 0)) == 15){
    dates4 <- c(dates4,i)  
  }
}
idates4 <- init_dates[sort(unique(dates4))]
intersect(idates3,idates4)
