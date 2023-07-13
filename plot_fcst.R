## Load modules
rm(list=ls())
require(SIDFEx)
require(spheRlab)

plot.dir = "/home/anjost001/Documents/AWI/Bachelorarbeit"

index = sidfex.load.index()
tid = c("900120")
iy = 2022
idoy = 233
subind.eccc = sidfex.fcst.search.extractFromTable(index = index, tid = tid, gid = "eccc001")
plot.dir.tid = file.path(plot.dir, tid)
if (!dir.exists(plot.dir.tid)){ # create directory if it does not exist yet
  dir.create(plot.dir.tid)
}
obs = sidfex.read.obs(TargetID = tid)
subind.tid = subind.eccc[subind.eccc$TargetID == tid,] # get index for current TargetID
subind.tid2 = subind.tid[subind.tid$InitYear == iy,]
subind.tid3 = subind.tid2[subind.tid2$InitDayOfYear == idoy,] # only wanted fcst

for(i.fcst in 1:length(subind.tid3$File)) {
  fcst = sidfex.read.fcst(files = subind.tid3$File[i.fcst], ens.merge = F)
  file.name = paste0("trajectory_" ,subind.tid3$File[i.fcst], ".pdf") # filename for saving
  if (file.exists(file.path(plot.dir.tid,file.name))){ # skip if file exists already
    print(paste0(file.path(plot.dir.tid,file.name), " exists already!"))
    next
  }
  obs.remap = sidfex.remaptime.obs2fcst(obs = obs, fcst = fcst)
  if (any(is.na(obs.remap$res.list[[1]]$data$Lon))){
    print("No obervations during forecast time, so no remapping possible. Skipping...")
    next
  }
  fcst = sidfex.rot.fcst(fcst = fcst) # rotate forecast so that initial position matches
  if (all(is.na(fcst$res.list[[1]]$data$Lat)) == TRUE){ # skip if there were only NAs
    next
  }
  fcst.reltime.max = fcst$res.list[[1]]$DaysForecastLength # last day
  fcst.daily = sidfex.remaptime.fcst(fcst = fcst, newtime.DaysLeadTime = 1:fcst.reltime.max) # forecast with daily resolution
  obs.remap.daily = sidfex.remaptime.obs2fcst(obs = obs, fcst = fcst.daily) # observation with dailyy resolution
  
  bc = sl.boundingcircle(lon = fcst$res.list[[1]]$data$Lon, lat =  fcst$res.list[[1]]$data$Lat) # circle enclosing all longitudes and latitudes
  lon.center = bc$center_lon # new north pole longitude (north pole = center of new/zoomed plot)
  lat.center = bc$center_lat # new north pole latitude
  rot.center = 0 # rotate around new north pole by this amount of degrees
  bounding.lat = 90-1.2*(bc$radius)
  
  pir = sl.plot.init(projection = "polar", polar.latbound = bounding.lat, polar.lonlatrot = c(lon.center,lat.center,rot.center), device = "pdf",width = 10, file.name = file.path(plot.dir.tid,file.name),main = paste0("\n",tid," ", fcst$res.list[[1]]$InitYear,"-",fcst$res.list[[1]]$InitDayOfYear, " lead time: ",fcst$res.list[[1]]$DaysForecastLength, " days")) # start plot
  sl.plot.naturalearth(pir, what = "land", resolution = "medium", lwd = 2, lines.col = "black") # plot coastline
  sl.plot.lonlatgrid(pir,pole.hole = TRUE, labels = TRUE, col = "grey", lwd = 2, labels.cex = 1.5) # plot grid
  # obs
  sl.plot.lines(pir, lon = obs.remap$res.list[[1]]$data$Lon, lat =  obs.remap$res.list[[1]]$data$Lat, col = "red", lwd=4) # plot observations
  sl.plot.points(pir, lon = obs.remap.daily$res.list[[1]]$data$Lon, lat =  obs.remap.daily$res.list[[1]]$data$Lat, pch = 19, col = "red") # circle for every day of obs
  sl.plot.points(pir, lon = obs.remap$res.list[[1]]$data$Lon[1], lat =  obs.remap$res.list[[1]]$data$Lat[1], col = "black", pch = 4, cex = 6) # cross for first day
  # fcst
  col.fcst = "blue"
  sl.plot.lines(pir, lon = fcst$res.list[[1]]$data$Lon, lat =  fcst$res.list[[1]]$data$Lat, col = col.fcst, lwd=4) # plot forecast trajectory
  sl.plot.points(pir, lon = fcst.daily$res.list[[1]]$data$Lon[2:(fcst.reltime.max)], lat =  fcst.daily$res.list[[1]]$data$Lat[2:(fcst.reltime.max)], col = col.fcst, pch = 19) # cross at each day of forecast
  legend("topleft",legend = c("Obs", fcst$res.list[[1]]$MethodID,"New Day","Start"),col = c("red",col.fcst,rep("black",2)), pch = c(NA,NA,16,4),lty = c(rep("solid",2),rep(NA,2)), lwd = c(3,3,NA,NA), cex = 1.8) # legend
  sl.plot.end(pir) # stop plot
  
}
