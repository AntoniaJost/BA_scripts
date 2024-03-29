rm(list=ls())

require(SIDFEx)
require(spheRlab)
require(RColorBrewer)

# for whatever reason 300234065498190 cannot be plotted for 2023-80-86 or 80-112

tid = c("900120")#, "900121", "900126", "900128", "300234065498190", "300234067527540")#, "300534063809090", "300534063807110")
# change tid to tids to use for loop
# for(i.single in 1:length(tids)){ # quick & dirty way to not enter every tid by hand
#   tid = tids[i.single]
#   print(tid)

# set your time preferences
d.start = 46
y.start = 2023
d.end = 60
y.end = 2023
# I recommend using "arctic" for plotting several buoys and "buoy" for single buoys
zoom = "buoy" # "buoy" or "arctic"

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

firsty = obs_red$data$Year[1]
firstdoy = obs_red$data$DOY[1]
lasty = obs_red$data$Year[length(obs_red$data$Year)]
lastdoy = obs_red$data$DOY[length(obs_red$data$DOY)]
daysPassed = sidfex.ydoy2reltime(Year = lasty, DayOfYear = lastdoy, RefYear = firsty, RefDayOfYear = firstdoy)
days = c(0:daysPassed)
ydoy = sidfex.reltime2ydoy(reltime = days, RefYear = firsty, RefDayOfYear = firstdoy)
obs.remap = sidfex.remaptime.obs(obs_red, newtime.YearDayOfYear = ydoy)
obs.daily = obs.remap
obs.daily$data = obs.daily$data[-1,] # remove first data point so it doesn't overlap with cross for first day

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
  bounding.lat = 73
}

# creating structure of map
if(length(tid) == 1) {
  file.name = paste0(tid, "_", y.start, ":", d.start, "-", d.end, "_obs.pdf")
} else {
  file.name = paste0("obs_", y.start, ":", d.start, "-", d.end, ".pdf")
}
title = paste0("\n", file.name)
plot.dir = "/home/anjost001/Documents/AWI/Bachelorarbeit/data_analysis/obs/2023:46-74/"

if (file.exists(file.path(plot.dir, file.name))) {
  stop(paste0("File '", file.name, "' already exists. Stopped.")) # skip if file already exists
} else {
  pir = sl.plot.init(projection = "polar", polar.latbound = bounding.lat, polar.lonlatrot = c(lon.center,lat.center,rot.center), do.init = T, main = title, file.name = file.path(plot.dir, file.name), device = "pdf") # build map
  sl.plot.naturalearth(pir, what = "land", resolution = "medium", lwd = 2, lines.col = "black") # plot coastlines
  sl.plot.lonlatgrid(pir, pole.hole = TRUE, labels = TRUE, col = "grey", labels.col = "black", lwd = 3, labels.cex = 2) # plot grid lines
  
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
    sl.plot.lines(pir, lon = lon.obs, lat = lat.obs, col = color[i.obs], lwd = 4) # plot trajectory
    if (zoom == "buoy") {
      cross = 6
    } else if (zoom == "arctic") {
      cross = 1
    }
    sl.plot.points(pir,lon = lon.obs[1],lat = lat.obs[1], pch = 4, cex = cross) # 6 for buoy # plot cross at first point of trajectory
    if(zoom == "buoy") {
      sl.plot.points(pir, lon = obs.daily$data$Lon, lat = obs.daily$data$Lat, col = color[i.obs], pch = 19, cex = 2) # point for every new day
    }
    labelTid = tid[[i.obs]]         # save TargetID to later display in legend
    label = c(labelTid, label)               # save as array, accessible outside of loop
    saveColor = color[i.obs]                 # save color to later display in legend
    colorLegend = c(saveColor, colorLegend)  # save as array, accessible outside of loop
    
  }
  
  sl.plot.end(pir, do.close.device = F) # ending plot before plotting the legend to prevent overlapping
  if(zoom == "arctic") {
    legend("bottomright", legend = c(label, "start"), col = c(colorLegend, "black"), pch = c(rep(NA, length(tid)), 4), lty = c(rep(1, length(tid)), 0), bg = "white", cex = 1.8) # legend
  } else if(zoom == "buoy") {
    legend("bottomright", legend = c(label, "new day", "start"), col = c(colorLegend, "red", "black"), pch = c(rep(NA, length(tid)), 16, 4), lty = c(rep(1, length(tid)), 0, 0), bg = "white", cex = 1.8) # legend
  }
  
  dev.off()
}
#}