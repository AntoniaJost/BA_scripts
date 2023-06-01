# library(devtools)
# install_github("helgegoessling/SIDFEx")
##########
# Rscript to demonstrate how to rotate a SIDFEx trajectory to the North pole and to align
# the 0h and 24h points of the forecast to match the observed positions

# if desired, clear workspace
rm(list=ls())
load("~/Downloads/sidfex.evaluate.R")

# load packages (if not already loaded)
require(spheRlab)
require(SIDFEx)

# setting wd
if(Sys.info() ["user"] == "anjost001") {
  wd = setwd("/home/anjost001/Documents/AWI/Bachelorarbeit/data_analysis/subdaily_res") # wd on Antonia's machine
} else { # wd on everyone else's machine
  work_dir = getwd() # get current working directory
  save_dir = "subdaily_res" # name of (new) folder where all data is going to be stored
  if (!file.exists(save_dir)) { # new folder is only being created if doesn't exist already
    dir.create(save_dir)
    print(paste0("new directory created:", work_dir, "/", save_dir))
  }
  wd = setwd(save_dir) # set new wd
}

# input
tids = c("900120", "900121")
gid = "eccc001"
mid = "giops"
# if you don't want any date specification, please remove parameters iy and doy completely from subind.fcst (doesn't work with empty string)
iy = "2021"
idoy = "50" 
plot = FALSE # if plots are wanted (might take quite long)
matrix = TRUE # if matrix with errors is wanted

# update fcst and obs data (set TRUE if desired)
if (FALSE) {
  res = sidfex.download.fcst()
  res2 = sidfex.download.obs()
}

# load SIDFEx index and read data
indx = sidfex.load.index()

# Loop over TargetIDs
for (i.tid in 1:length(tids)) {
  tid = tids[i.tid]
  
  subind.fcst = sidfex.fcst.search.extractFromTable(index = indx, tid = tid, gid = gid, mid = mid) #, iy = iy, idoy = idoy) # delete iy and idoy if wanted
  fcst = sidfex.read.fcst(subind.fcst)

  # so sieht aus, wenn iy und idoy nicht angegeben
  # subind.all = sidfex.fcst.search.extractFromTable(index = indx, tid = tid, gid = gid, mid = mid)
  # fcst.all = sidfex.read.fcst(subind.all)
  
  # obs1 = sidfex.read.obs(index = indx, TargetID = tid)
  obs = sidfex.remaptime.obs2fcst(fcst = fcst)

  if (length(fcst$res.list) != length(obs$res.list)) { # quick sanity check
    stop("check fcst and obs, they don't have same length")
  }
  ### Loop over every day in res.list
  for (i in 1:length(fcst$res.list)) {
    
    # dividing obs and fcst datasets of 10 days leadtime into 10 separate datasets for each day
    daysLeadTime = fcst$res.list[[i]]$data$DaysLeadTime # extract Leadtime
    integers = unique(floor(daysLeadTime)) # reduce to integers
    
    # function to create subset
    subdiv_dataset = function(dataset) { 
    
      teilDatensaetze = list() # empty list to store sub datasets
      for (i.num in 1:length(integers)) { # loop over integers (=every day)
        currentNumber = integers[i.num]
        startIdx = which(floor(daysLeadTime) == currentNumber) # index for begin of dataset
          
        if(i.num == length(integers)-1) { # if last set (day 9 to 10) add last single value of 10.000 to dataset of day 9
          startIdx = c(startIdx, startIdx[length(startIdx)]+1)
          teilDatensatz = dataset$res.list[[i]]$data[startIdx[1]:startIdx[length(startIdx)],] 
          teilDatensaetze[[i.num]] = teilDatensatz 
          break # stop loop
        }
        
        teilDatensatz = dataset$res.list[[i]]$data[startIdx[1]:startIdx[length(startIdx)],] # extract
        teilDatensaetze[[i.num]] <- teilDatensatz # save
      }
      return(teilDatensaetze)
    }
    
    obs_subdiv = subdiv_dataset(obs) # subdivided obs dataset
    fcst_subdiv = subdiv_dataset(fcst) # subdivided fcst dataset
    
    # loop over every of the sub data sets
    for (i.sub in 1:length(obs_subdiv)) {
      if (i == 49 && i.sub == 9 || i == 234 && i.sub == 10) {
        i.sub=- i.sub + 1
        next  # skip error plot (too close to north pole, creating NaNs and Infs)
      }
      # bringing subdivided fcst and obs datasets back to "normal" format to use SIDFEx and spherelab functions
      obs_subdiv_temp = obs # copy normal obs for usual format
      obs_subdiv_temp$res.list = list(obs_subdiv_temp$res.list[[i]])
      obs_subdiv_temp$res.list[[1]]$data = NULL # remove any data
      obs_subdiv_temp$res.list[[1]]$data = as.data.frame(append(obs_subdiv_temp$res.list[[1]]$data, obs_subdiv[[i.sub]])) # "overwrite" with reduced dataset of just one day leadtime
      # same for fcst
      fcst_subdiv_temp = fcst # hier beginnt Fehler! Darf nicht einfach res.list[[1]] Sachen (vor Data) von fcst uebernehmen, denn das wird immer das Alte sein
      fcst_subdiv_temp$res.list = list(fcst_subdiv_temp$res.list[[i]]) # hier mit fcst$res.list[[i]] ueberschreiben?
      fcst_subdiv_temp$res.list[[1]]$data = NULL
      fcst_subdiv_temp$res.list[[1]]$data = as.data.frame(append(fcst_subdiv_temp$res.list[[1]]$data, fcst_subdiv[[i.sub]]))
    
      # determine fcst & obs rotation parameters from initial position, then rotate
      fcst.init = as.data.frame(lapply(fcst_subdiv_temp$res.list[[1]]$data, function(x) x[1]))
      obs.init  = as.data.frame(lapply(obs_subdiv_temp$res.list[[1]]$data, function(x) x[1]))
      
      # rotate to north pole
      fcst.abg = sl.lonlatrot2abg(lonlatrot = c(fcst.init$Lon,fcst.init$Lat,0))
      fcst.rot = sl.rot(lon = fcst_subdiv_temp$res.list[[1]]$data$Lon, lat = fcst_subdiv_temp$res.list[[1]]$data$Lat,
                        alpha = fcst.abg[1], beta = fcst.abg[2], gamma = fcst.abg[3])
      # determine obs rotation parameters from obs position at initial time, then rotate
      obs.abg = sl.lonlatrot2abg(lonlatrot = c(obs.init$Lon,obs.init$Lat,0))
      obs.rot = sl.rot(lon = obs_subdiv_temp$res.list[[1]]$data$Lon, lat = obs_subdiv_temp$res.list[[1]]$data$Lat,
                       alpha = obs.abg[1], beta = obs.abg[2], gamma = obs.abg[3])
      
      # rotate so that the position after one day is at the Greenwich meridian
      # first the obs
      obs.rot.1d = list(Lat = numeric(), Lon = numeric())
      obs.rot.1d$Lat = obs.rot$lat[length(obs.rot$lat)] # get last data pair (lon,lat) of the day
      obs.rot.1d$Lon = obs.rot$lon[length(obs.rot$lon)]
      obs.rot$lon = obs.rot$lon - obs.rot.1d$Lon # substract previously extracted lon after 1 day to be at 0 (=Greenwich meridian)
      
      # now the fcst
      fcst.rot.1d = list(Lat = numeric(), Lon = numeric())
      fcst.rot.1d$Lat = fcst.rot$lat[length(fcst.rot$lat)]
      fcst.rot.1d$Lon = fcst.rot$lon[length(fcst.rot$lon)]
      fcst.rot$lon = fcst.rot$lon - fcst.rot.1d$Lon
      
      # now stretch or squeeze the fcst to match the observed latitude after one day
      fcst.rot$lat = 90 - (((90 - obs.rot.1d$Lat) / (90 - fcst.rot.1d$Lat)) * (90 - fcst.rot$lat))
      
      # create a "standard" fcst and obs object again with the adjusted trajectories so that sidfex.evaluate() can be used
      # note: the header information is not adjusted, just the trajectory in the data object!
      fcst.adj = fcst_subdiv_temp
      fcst.adj$res.list[[1]]$data$Lat = fcst.rot$lat
      fcst.adj$res.list[[1]]$data$Lon = fcst.rot$lon
      obs.adj = obs_subdiv_temp
      obs.adj$res.list[[1]]$data$Lat = obs.rot$lat
      obs.adj$res.list[[1]]$data$Lon = obs.rot$lon
      
      # create a "reference" forecast that drifts linearly from 0h to 24h (and beyond)
      fcst.2points = fcst.adj
      fcst.2points$res.list[[1]]$data = fcst.2points$res.list[[1]]$data[c(1, nrow(fcst.2points$res.list[[1]]$data)), ] # get first and last data point of the day (= 0-23.5.. h)
      
      var = sidfex.ydoy2reltime(Year = fcst.2points$res.list[[1]]$data$Year[1], DayOfYear = fcst.2points$res.list[[1]]$data$DayOfYear[1], RefYear = fcst.2points$res.list[[1]]$InitYear, RefDayOfYear = fcst.2points$res.list[[1]]$InitDayOfYear) # variable to store difference (in days) between InitDOY and DOY
      # some printing for error hunting
      print(paste0("i:", i))
      print(paste0("i.sub:", i.sub))
      print(var)
      print(paste0("DOY:", fcst.2points$res.list[[1]]$data$DayOfYear[1]))
      print(paste0("DOY_Init:", fcst.2points$res.list[[1]]$InitDayOfYear))
      
      fcst.lin = sidfex.remaptime.fcst(fcst = fcst.2points, newtime.DaysLeadTime = (fcst.adj$res.list[[1]]$data$DaysLeadTime), extrapolate = TRUE) # create 48 corresponding data points in between for better comparison
      # important! newtime.DaysLeadTime NEEDS parameter var because it always takes Init Value of fcst, not res.list[[1]]$data value!!
      
      # evaluation
      fcst.adj.eval = sidfex.evaluate(obs = obs.adj, fcst = fcst.adj, do.speedangle = TRUE, verbose=FALSE)
      fcst.lin.eval = sidfex.evaluate(obs = obs.adj, fcst = fcst.lin, do.speedangle = TRUE, verbose=FALSE)
      
      # some plotting
      if (plot = TRUE) {
        
        plot.title = paste0(fcst$res.list[[i]]$TargetID, "_", fcst$res.list[[i]]$GroupID, "_", fcst$res.list[[i]]$MethodID, "_", fcst$res.list[[i]]$InitYear, "-", fcst$res.list[[i]]$InitDayOfYear, "_ld:", i.sub, ".png")
        if (!file.exists(file.path(wd, plot.title))) {
          png(filename = plot.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
        } else {
          print(paste0("File '", plot.title, "' already exists. Skipped it.")) # skip if file already exists
          next
        }
        
        par(mfrow=c(2,2))
        
        # calculating ylim for plots
        range_calc = function(eval_value) {
          ylim = c()
          is_naN_inf = is.na(eval_value) | is.nan(eval_value) | is.infinite(eval_value)
          ylim = range(eval_value[!is_naN_inf])
          return(ylim)
        }
        
        ylim.lin.speed = range_calc(fcst.lin.eval$res.list[[1]]$ens.mean.relspeed)
        ylim.lin.angle = range_calc(fcst.lin.eval$res.list[[1]]$ens.mean.angle)
        ylim.lin.gc = range_calc(fcst.lin.eval$res.list[[1]]$ens.mean.gc.dist)
        ylim.adj.speed = range_calc(fcst.adj.eval$res.list[[1]]$ens.mean.relspeed)
        ylim.adj.angle = range_calc(fcst.adj.eval$res.list[[1]]$ens.mean.angle)
        ylim.adj.gc = range_calc(fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist)
        
        range_gc = range(ylim.lin.gc, ylim.adj.gc)
        range_speed = range(ylim.lin.speed, ylim.adj.speed)
        range_angle = range(ylim.lin.angle, ylim.adj.angle)
        
        # plot trajectories in x-y-space
        fcst.adj.xyz = sl.lonlat2xyz(lon=fcst.adj$res.list[[1]]$data$Lon, lat=fcst.adj$res.list[[1]]$data$Lat)
        fcst.lin.xyz = sl.lonlat2xyz(lon=fcst.lin$res.list[[1]]$data$Lon, lat=fcst.lin$res.list[[1]]$data$Lat)
        obs.adj.xyz = sl.lonlat2xyz(lon=obs.adj$res.list[[1]]$data$Lon, lat=obs.adj$res.list[[1]]$data$Lat)
        xyz.1d = sl.lonlat2xyz(lon=obs.rot.1d$Lon, lat=obs.rot.1d$Lat)  # position after 1 day, so that the plot can be scaled reasonably
        dist.1d = sqrt(xyz.1d$x^2 + xyz.1d$y^2)
        
        plot(NA,xlim=c(-0.3,1)*dist.1d*1.1, ylim=c(-1,1)*dist.1d*1.1,xlab="x / Earth radius",ylab="y / Earth radius")
        #abline(h=0,v=0,col="grey",lty=3)
        points(x=obs.adj.xyz$x, y=obs.adj.xyz$y, col="black",cex=0.3)
        points(x=fcst.lin.xyz$x, y=fcst.lin.xyz$y, col="grey",cex=0.3)
        points(x=fcst.adj.xyz$x, y=fcst.adj.xyz$y, col="red",cex=0.3)
        points(x=dist.1d, y=0, pch="+", col="orange",cex=2)
        legend("topleft", c("obs adjusted", "fcst linear", "fcst adjusted", "1 day"), col = c("black","grey","red","orange"), pch = c(1,1,1,3), bty = "n", cex = 0.7)
        
        # plot evaluation results
        # great-circle distance
        plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.gc.dist, xlab="days lead time", ylab = "great-circle distance / m",
             xlim=c(i.sub-1.05,i.sub+0.05), ylim = range_gc, col = "grey")
        abline(h=0,v=i.sub,col="grey",lty=3)
        points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist, col="red")
        legend("topleft", c("Error betw. obs & daily fcst", "Error betw. obs & hh fcst"), col = c("grey", "red"), pch = c(1,1), bty = "n", cex = 0.7)
        
        # relative speed
        plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.relspeed, xlab="days lead time", ylab = "relative speed",
             xlim=c(i.sub-1.05,i.sub+0.05), ylim=range_speed, col="grey")
        points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.relspeed, col="red")
        abline(h=1,v=i.sub,col="grey",lty=3)
        
        # relative angle
        plot(x=fcst.lin$res.list[[1]]$data$DaysLeadTime, y=fcst.lin.eval$res.list[[1]]$ens.mean.angle, xlab="days lead time", ylab = "relative angle / degree left",
             xlim=c(i.sub-1.05,i.sub+0.05), ylim=range_angle, col="grey")
        abline(h=0,v=i.sub,col="grey",lty=3)
        points(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.angle, col="red")
        
        # add title
        par(oma = c(0, 0, 3, 0))
        title(main=paste0("Analysis subdaily resolution \n ", fcst$res.list[[1]]$TargetID, "_", fcst$res.list[[1]]$GroupID, "_", fcst$res.list[[1]]$MethodID, "_", fcst$res.list[[1]]$InitYear, "-", fcst$res.list[[1]]$InitDayOfYear), outer = T)
        
        dev.off()
      } # end if-clause for plotting
      
      # some statistics
      if(matrix = TRUE) {
        
      }
      
    } # end for-loop sub dataset 

  } # end for-loop days in res.list
  
} # end for-loop TargetID

