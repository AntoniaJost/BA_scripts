rm(list=ls())

require(SIDFEx)
require(spheRlab)
require(RColorBrewer)

# for whatever reason 300234065498190 cannot be plotted for 2023-80-86 or 80-112

tid = c("900126")
# change tid to tids to use for loop
# for(i.single in 1:length(tids)){ # quick & dirty way to not enter every tid by hand
#   tid = tids[i.single]
#   print(tid)

# set your time preferences
d.start = 80
y.start = 2023
d.end = 112
y.end = 2023
# I recommend using "arctic" for plotting several buoys and "buoy" for single buoys
zoom = "buoy" # "buoy"

obs = sidfex.read.obs(TargetID = tid)

if(length(tid) == 1) {
  list_obs = list(obs)
  curr_obs = list_obs[[1]]
}

lon.obs.all = c()
lat.obs.all = c()

for(i in 1:length(tid)) {
 if(length(tid) > 1) {
   curr_obs = obs[[i]]
 }
  
  first = which(curr_obs$data$Year == y.start & trunc(curr_obs$data$POS_DOY) == d.start)[1] # get index of first line that matches date input
  last_temp = which(curr_obs$data$Year == y.end & trunc(curr_obs$data$POS_DOY) == d.end)
  last = tail(last_temp, n = 1) # get index of last line
  
  obs_red = curr_obs
  obs_red$data = obs_red$data[first:last,] # reduced observational data set, only for selected time frame
  
  lon.obs = obs_red$data$Lon[!is.na(obs_red$data$Lon)] # extract lon
  lat.obs = obs_red$data$Lat[!is.na(obs_red$data$Lat)] # extract lat
  lon.obs.all = c(lon.obs, lon.obs.all)
  lat.obs.all = c(lat.obs, lat.obs.all)
}

bc = sl.boundingcircle(lon.obs.all,lat.obs.all) # returns parameters for smallest region which contains every coordinate
if(zoom == "buoy") {
  # plot centered on bouy
  lon.center = bc$center_lon # new north pole longitude (north pole = center of new/zoomed plot)
  lat.center = bc$center_lat # new north pole latitude
  rot.center = 0 # rotate around new north pole by this amount of degrees
  bounding.lat = 90-1.2*bc$radius # bounding latitude of zoomed plot: radius of boundingcircle, (normally multiplied by a factor (1.2) to make sure that it is actually a bit larger than the trajectory region, now deactivated bc polar.latbound must be within [0,90])
  if (bounding.lat < 0) {
    bounding.lat = 0
  }
} else if (zoom == "arctic") {
  # Arctic-wide plot
  lon.center = 0 
  lat.center = 90 
  rot.center = 0 
  bounding.lat = 65
}

# creating structure of map
if(length(tid) == 1) {
  file.name = paste0(tid, "_", y.start, ":", d.start, "-", d.end, "_obs.png")
} else {
  file.name = paste0("obs_", y.start, ":", d.start, "-", d.end, ".png")
}
title = paste0("\n", file.name)
plot.dir = "/home/anjost001/Downloads"#"/home/anjost001/Documents/AWI/Bachelorarbeit/data_analysis/obs"

if (file.exists(file.path(plot.dir, file.name))) {
  stop(paste0("File '", file.name, "' already exists. Stopped.")) # skip if file already exists
} else {
  par(oma=c(0,0,2,0))
  pir = sl.plot.init(projection = "polar", polar.latbound = bounding.lat, polar.lonlatrot = c(lon.center,lat.center,rot.center), do.init = T, file.name = file.path(plot.dir, file.name), main = title, device = "png") # build map
  sl.plot.naturalearth(pir, what = "land", resolution = "coarse", lwd = 0.5, lines.col = "black") # plot coastlines
  sl.plot.lonlatgrid(pir, pole.hole = TRUE, labels = TRUE, col = "grey", labels.col = "black") # plot grid lines
  
  label = c()
  colorLegend = c()
  
  for(i.obs in 1:length(tid)) {
    if(length(tid) > 1) {
      curr_obs = obs[[i.obs]]
    }
    first = which(curr_obs$data$Year == y.start & trunc(curr_obs$data$POS_DOY) == d.start)[1]
    last_temp = which(curr_obs$data$Year == y.end & trunc(curr_obs$data$POS_DOY) == d.end)
    last = tail(last_temp, n = 1)
    
    obs_red = curr_obs
    obs_red$data = obs_red$data[first:last,] # reduced observational data set, only for selected time frame
    
    lon.obs = obs_red$data$Lon[!is.na(obs_red$data$Lon)]
    lat.obs = obs_red$data$Lat[!is.na(obs_red$data$Lat)]
  
    ### plot trajectories
    color = colorRampPalette(brewer.pal(8, "RdYlBu"))(length(tid)) # colorblind-friendly color scale
    sl.plot.lines(pir, lon = lon.obs, lat = lat.obs, col = color[i.obs]) # plot trajectory
    sl.plot.points(pir,lon = lon.obs[1],lat = lat.obs[1], pch = 4) # plot square at last point of trajectory
    labelTid = tid[[i.obs]]         # save TargetID to later display in legend
    label = c(labelTid, label)               # save as array, accessible outside of loop
    saveColor = color[i.obs]                 # save color to later display in legend
    colorLegend = c(saveColor, colorLegend)  # save as array, accessible outside of loop
    
  }
  
  sl.plot.end(pir, do.close.device = F) # ending plot before plotting the legend to prevent overlapping
  legend("bottomright", legend = c(label, "start"), col = c(colorLegend, "black"), pch = c(rep(NA, length(tid)), 4), lty = c(rep(1, length(tid)), 0), bg = "white") # legend
  
  dev.off()
}
#}