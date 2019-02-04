library(ggplot2)
require(maps)
require(maptools)

worldmap <- map("world", fill = TRUE, plot=FALSE)

stationycoor <- c(52.933, 52.317, 52.65, 53.383, 52.5, 52.1, 52.9, 52.45, 53.217, 52.7, 52.05, 53.417, 52.433, 52.75, 53.117, 
                  52.067, 53.2, 52.267, 51.45, 51.217, 51.983, 51.967, 51.971, 51.567, 51.85, 51.45, 51.65, 51.2, 50.9, 51.5)
stationxcoor <- c(4.783, 4.783, 4.983, 5.35, 4.6, 5.183, 5.383, 5.517, 5.75, 5.883, 5.867, 6.2, 6.267, 6.567, 6.583, 6.65, 7.15, 
                  6.883, 3.6, 3.867, 4.117, 4.45, 4.927, 4.933, 5.15, 5.383, 5.7, 5.767, 5.767, 6.2)

rast <- data.frame(longitude = stationxcoor, latitude = stationycoor)

tempPlot <- ggplot(rast) + 
  geom_point(aes(x = longitude, y =  latitude)) +
  geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = NA, color = "black", size = 0.1, linetype = "dashed") +
  scale_fill_distiller("", palette = "Spectral", direction=1, na.value = "white") + theme_bw() +
  coord_fixed(ratio = 1.3, xlim = c(3.25, 7.25), ylim = c(50.75, 53.5))
print(tempPlot)


ggsave(paste0("stations_NL2.pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/",  width = 12.5, height = 10, units = "cm")


obs <- function(x){if (x > 0.6){return(1)} else {return(0)}}
forec <- function(x){return(0.5*(1 + tanh((x-0.5)*6)))}

x <- seq(0,1,0.01)
observations <- array(NA, c(length(x)))
forecasts <- array(NA, c(length(x)))
for (i in 1:length(x)){
observations[i] <- obs(x[i])
forecasts[i] <- forec(x[i])
}
tempPlot <- ggplot(data.frame(xaxis = x, data = c(observations, forecasts), Type = rep(c("observations", "forecasts"), each = length(x)))) +
  geom_line(aes(x = xaxis, y = data, color = Type)) +
  scale_colour_manual(values = c("black","blue")) +
  geom_ribbon(data = data.frame(xaxis = x, observations, forecasts), mapping = aes(x = xaxis, ymin = observations, ymax = forecasts), fill="gray", alpha="0.5") +
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) + 
  ylab("Cumulative density function")
print(tempPlot)


ggsave(paste0("CRPS_illustration.pdf"), plot = tempPlot, device = "pdf", path = "/usr/people/bakker/kilianbakker/plots/",  width = 12.5, height = 10, units = "cm")
