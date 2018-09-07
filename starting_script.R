# set library paths:
libpath = "/net/pc132059/nobackup_1/users/whan/R/x86_64-redhat-linux-gnu-library/3.3"
.libPaths(c(libpath, .libPaths()))
library(Rgrib2)
library(rasterVis)
library(ggplot2)
require(maps)
require(maptools)

# define a map for plotting later:
worldmap <- map("world", fill = TRUE, plot=FALSE)

yearNumber <- 2017
indir     <- paste0("/net/pc150388/nobackup_3/users/schmeits/HARM38/",yearNumber)
fnames    <- data.frame(files = list.files(path = indir, full.names = TRUE), stringsAsFactors = FALSE)
init_dates  <- unique(gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(fnames$files))))

####################################
#Checking installed packages
####################################

ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
print(ip, row.names=FALSE)

####################################
# An example importing a single file:
####################################

# make a data.frame with the grib information (same for all files):
grib_info <- Rgrib2::Gopen(fnames$files[1])

# open a single grib file:
gfile     <- Gdec(fnames$files[1], 65)
rast      <- flip(raster(t(gfile), 
                          xmn=attributes(gfile)$domain$SW[1], 
                          xmx=attributes(gfile)$domain$NE[1],
                          ymn=attributes(gfile)$domain$SW[2], 
                          ymx=attributes(gfile)$domain$NE[2]), 
                   direction="y")

# you can also plot a raster with a map:
gplot(rast) + 
  geom_tile(aes(fill = value)) +
  geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = NA, color = "black", size = 0.1, linetype = "dashed") +
  scale_fill_distiller("", palette = "Spectral", direction=1, na.value = "white") +
  coord_fixed(xlim = c(0, 11), ylim = c(49, 55))

# write a function that imports the raster:
import_grib <- function(afile, param, adate, alt){
  tmpi      <- Gopen(afile)
  
  tmpg      <- Gdec(afile, param)
  tmpr      <- flip(raster(t(tmpg), 
                            xmn=attributes(tmpg)$domain$SW[1], 
                            xmx=attributes(tmpg)$domain$NE[1],
                            ymn=attributes(tmpg)$domain$SW[2], 
                            ymx=attributes(tmpg)$domain$NE[2]), 
                     direction="y")
  names(tmpr) <- paste0(tmpi$shortName[which(tmpi$position == param)], "_", adate, "_", alt)
  return(tmpr)
}

# run the function like this:
newraster <- import_grib(afile=fnames[grepl("20170416", fnames$files), ][[13]], param = 275, adate = "20180110", alt = "01")
newraster <- crop(newraster, extent(3.4, 7.1, 50.8, 53.5))

gplot(newraster) + 
  geom_tile(aes(fill = value)) +
  geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = NA, color = "black", size = 0.1, linetype = "dashed") +
  scale_fill_distiller("", palette = "Spectral", direction=1, na.value = "white") +
  coord_fixed(xlim = c(0, 11), ylim = c(49, 55))

###################################################################################
# Use an 'apply' function to import all leadtimes for one date (similar to a loop):
###################################################################################

## functions to extract the date and leadtime
get_date <- function(afile){
  gsub("0000_[0-9]{5}_GB", "", gsub("HA38_N25_SOLFC_", "", basename(afile)))
}
get_lt <- function(afile){
  gsub("00_GB", "", gsub("HA38_N25_SOLFC_[0-9]{12}_", "", basename(afile)))
}

# change this value to import different dates:
whichdate   <- 1

# import the rasters and stack them into a 'raster stack'
daterast    <- do.call(stack, apply(as.matrix(fnames[grepl(init_dates[whichdate], fnames$files), ]), 1, 
                                    function(d) import_grib(afile = d, param = 65, adate = get_date(d), alt = get_lt(d))))

# now we can plot different lead times:
plot(daterast, 2)
plot(daterast, 49)


##### extract data for some coordinates:
# create a SpatialPoints object for some lon-lat data and use the 'extract' function to get the data for this point:
pointdat      <- extract(daterast, SpatialPoints(cbind(4.927, 51.971)))