# library(devtools)
# install_github("helgegoessling/SIDFEx")
##########
# Rscript to demonstrate how to rotate a SIDFEx trajectory to the North pole and to align
# the 0h and 24h points of the forecast to match the observed positions

# if desired, clear workspace
rm(list=ls())

# load packages (if not already loaded)
require(spheRlab)
require(SIDFEx)
library(gridExtra)

# update fcst and obs data (set TRUE if desired)
if (FALSE) {
  res = sidfex.download.fcst()
  res2 = sidfex.download.obs()
}

# load SIDFEx index and read data
indx = sidfex.load.index()
#fcst = sidfex.read.fcst(files = indx[indx$File=="eccc001_giops_900128_2023-113_001",])
fcst = sidfex.read.fcst(files = indx[indx$File=="eccc001_giops_900120_2021-033_001",])
obs = sidfex.remaptime.obs2fcst(fcst = fcst)

###################### TESTING ZONE ##########################
# Spalte "DaysLeadTime" extrahieren
daysLeadTime <- fcst$res.list[[1]]$data$DaysLeadTime

# Ganze Zahlen in der Spalte "DaysLeadTime" identifizieren
wholeNumbers <- unique(floor(daysLeadTime))

# Liste für die Teildatensätze erstellen
teilDatensaetze <- list()

# Schleife über die ganzen Zahlen
for (i in 1:length(wholeNumbers)) {
  # Ganze Zahl
  currentNumber <- wholeNumbers[i]
  
  # Index für den Beginn des Teildatensatzes
  startIdx <- which(floor(daysLeadTime) == currentNumber)
    
  if(i == length(wholeNumbers)-1) {
    startIdx = c(startIdx, startIdx[length(startIdx)]+1)
    teilDatensatz <- fcst$res.list[[1]]$data[startIdx[1]:startIdx[length(startIdx)],] #[startIdx:endIdx, ] #startIdx[1]:startIdx[length(startIdx)]
    teilDatensaetze[[i]] <- teilDatensatz
    break
  }
  # Teildatensatz erstellen
  teilDatensatz <- fcst$res.list[[1]]$data[startIdx[1]:startIdx[length(startIdx)],] #[startIdx:endIdx, ] #startIdx[1]:startIdx[length(startIdx)]

  # Teildatensatz der Liste hinzufügen
  teilDatensaetze[[i]] <- teilDatensatz
}

##############################################################


fcst.init = sidfex.remaptime.fcst(fcst = fcst,newtime.DaysLeadTime = 0)
# in the previous line one could also use InitLon & InitLat directly, but this is to make sure all is consistent
obs.init = sidfex.remaptime.obs2fcst(fcst = fcst.init)

# rotate to north pole
# determine fcst rotation parameters from initial position, then rotate
fcst.abg = sl.lonlatrot2abg(lonlatrot = c(fcst.init$res.list[[1]]$data$Lon,fcst.init$res.list[[1]]$data$Lat,0))
fcst.rot = sl.rot(lon = fcst$res.list[[1]]$data$Lon, lat = fcst$res.list[[1]]$data$Lat,
                  alpha = fcst.abg[1], beta = fcst.abg[2], gamma = fcst.abg[3])
# determine obs rotation parameters from obs position at initial time, then rotate
obs.abg = sl.lonlatrot2abg(lonlatrot = c(obs.init$res.list[[1]]$data$Lon,obs.init$res.list[[1]]$data$Lat,0))
obs.rot = sl.rot(lon = obs$res.list[[1]]$data$Lon, lat = obs$res.list[[1]]$data$Lat,
                 alpha = obs.abg[1], beta = obs.abg[2], gamma = obs.abg[3])

# rotate so that the position after one day is at the Greenwich meridian
for (i in 1:obs$res.list[[1]]$data$DaysLeadTime[length(obs$res.list[[1]]$data$DaysLeadTime)]) {
  
}
# first the obs
obs.rot.1d = sl.trajectory.remaptime(oldtime = obs$res.list[[1]]$data$DaysLeadTime, newtime = i, # hier muesste i statt 1 stehen fuer dynamisch
                                     oldlat = obs.rot$lat, oldlon = obs.rot$lon) # extrahiert Wertepaar nach genau 1 Tag Leadtime
obs.rot$lon = obs.rot$lon - obs.rot.1d$Lon # damit latitude nach 1 Tag genau wieder 0 ist
# now the fcst
fcst.rot.1d = sl.trajectory.remaptime(oldtime = fcst$res.list[[1]]$data$DaysLeadTime, newtime = i,
                                      oldlat = fcst.rot$lat, oldlon = fcst.rot$lon)
fcst.rot$lon = fcst.rot$lon - fcst.rot.1d$Lon
# now stretch or squeeze the fcst to match the observed latitude after one day
fcst.rot$lat = 90 - (((90 - obs.rot.1d$Lat) / (90 - fcst.rot.1d$Lat)) * (90 - fcst.rot$lat))

# create a "standard" fcst and obs object again with the adjusted trajectories so that sidfex.evaluate() can be used
# note: the header information is not adjusted, just the trajectory in the data object!
fcst.adj = fcst
fcst.adj$res.list[[1]]$data$Lat = fcst.rot$lat
fcst.adj$res.list[[1]]$data$Lon = fcst.rot$lon
obs.adj = obs
obs.adj$res.list[[1]]$data$Lat = obs.rot$lat
obs.adj$res.list[[1]]$data$Lon = obs.rot$lon

# create a "reference" forecast that drifts linearly from 0h to 24h (and beyond)
fcst.2points = sidfex.remaptime.fcst(fcst = fcst.adj, newtime.DaysLeadTime = c(0,1))
fcst.lin = sidfex.remaptime.fcst(fcst = fcst.2points, newtime.DaysLeadTime = fcst.adj$res.list[[1]]$data$DaysLeadTime, extrapolate = TRUE)

# evaluation
fcst.adj.eval = sidfex.evaluate(obs = obs.adj, fcst = fcst.adj, do.speedangle = TRUE, verbose=FALSE)
fcst.lin.eval = sidfex.evaluate(obs = obs.adj, fcst = fcst.lin, do.speedangle = TRUE, verbose=FALSE)

# some plotting
par(mfrow=c(2,2))

# plot trajectories in x-y-space
fcst.adj.xyz = sl.lonlat2xyz(lon=fcst.adj$res.list[[1]]$data$Lon, lat=fcst.adj$res.list[[1]]$data$Lat)
fcst.lin.xyz = sl.lonlat2xyz(lon=fcst.lin$res.list[[1]]$data$Lon, lat=fcst.lin$res.list[[1]]$data$Lat)
obs.adj.xyz = sl.lonlat2xyz(lon=obs.adj$res.list[[1]]$data$Lon, lat=obs.adj$res.list[[1]]$data$Lat)
xyz.1d = sl.lonlat2xyz(lon=obs.rot.1d$Lon, lat=obs.rot.1d$Lat)  # position after 1 day, so that the plot can be scaled reasonably
dist.1d = sqrt(xyz.1d$x^2 + xyz.1d$y^2)

plot(NA,xlim=c(-1,1)*dist.1d*1.5, ylim=c(-1,1)*dist.1d*1.5,xlab="x / Earth radius",ylab="y / Earth radius")
#abline(h=0,v=0,col="grey",lty=3)
points(x=obs.adj.xyz$x, y=obs.adj.xyz$y, col="black",cex=0.3)
points(x=fcst.lin.xyz$x, y=fcst.lin.xyz$y, col="grey",cex=0.3)
points(x=fcst.adj.xyz$x, y=fcst.adj.xyz$y, col="red",cex=0.3)
points(x=dist.1d, y=0, pch="+", col="orange",cex=2)
legend("topleft", c("obs adjusted", "fcst linear", "fcst adjusted", "1 day"), col = c("black","grey","red","orange"), pch = c(1,1,1,3), bty = "n", cex = 0.7, y.intersp = 0.5)

# plot evaluation results (first 2 days only; note that the second day would be treated differently though...)
plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.gc.dist, xlab="days lead time", ylab = "great-circle distance / km",
     xlim=c(0,2), ylim=range(fcst.lin.eval$res.list[[1]]$ens.mean.gc.dist[fcst.lin$res.list[[1]]$data$DaysLeadTime<=2]), col="grey")
abline(h=0,v=1,col="grey",lty=3)
points(x=1, y=0, pch="+", col="orange",cex=2)
points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist, col="red")
legend("topleft", c("Error betw. obs & daily fcst", "Error betw. obs & hh fcst"), col = c("grey", "red"), pch = c(1,1), bty = "n", cex = 0.7, y.intersp = 0.5)

plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.relspeed, xlab="days lead time", ylab = "relative speed",
     xlim=c(0,2), ylim=range(fcst.lin.eval$res.list[[1]]$ens.mean.relspeed[fcst.lin$res.list[[1]]$data$DaysLeadTime<=2],na.rm = TRUE), col="grey")
points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.relspeed, xlim=c(0,2), col="red")
abline(h=1,v=1,col="grey",lty=3)
points(x=1, y=1, pch="+", col="orange",cex=2)

plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.angle, xlab="days lead time", ylab = "relative angle / degree left",
     xlim=c(0,2), ylim=range(fcst.lin.eval$res.list[[1]]$ens.mean.angle[fcst.lin$res.list[[1]]$data$DaysLeadTime<=2],na.rm = TRUE), col="grey")
abline(h=0,v=1,col="grey",lty=3)
points(x=1, y=0, pch="+", col="orange",cex=2)
points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.angle, xlim=c(0,2), col="red")
##########
par(oma = c(0, 0, 3, 0))
title(main=paste0("Analysis subdaily resolution \n ", fcst$res.list[[1]]$TargetID, "_", fcst$res.list[[1]]$GroupID, "_", fcst$res.list[[1]]$MethodID, "_", fcst$res.list[[1]]$InitYear, "-", fcst$res.list[[1]]$InitDayOfYear), outer = T)
