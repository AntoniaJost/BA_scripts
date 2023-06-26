# library(devtools)
# install_github("helgegoessling/SIDFEx")
##########
# Rscript to demonstrate how to rotate a SIDFEx trajectory to the North pole and to align
# the 0h and 24h points of the forecast to match the observed positions

# if desired, clear workspace
rm(list=ls())
source("/home/anjost001/Documents/BA_scripts/sidfex.evaluate.R")
source("/home/anjost001/Documents/BA_scripts/subset.R")
source("/home/anjost001/Documents/BA_scripts/range.calc.R")
source("/home/anjost001/Documents/BA_scripts/matrix.R")
source("/home/anjost001/Documents/BA_scripts/calc.mean.R")
source("/home/anjost001/Documents/BA_scripts/split.daily.R")

# load packages (if not already loaded)
require(spheRlab)
require(SIDFEx)
require(viridis)
require(matrixStats)

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
  wd = setwd(paste0(work_dir, "/" , save_dir)) # set new wd
}

# input
tids = c("300234065498190")
gid = "eccc001"
mid = "giops"
# if you don't want any date specification, please remove parameters iy and doy completely from subind.fcst (doesn't work with empty string)
iy = "2023"
idoy = c(as.character(46:74))

update = FALSE # if updated fcst and obs are wanted
plot = TRUE # if plots are wanted (might take quite long)
matrixs = FALSE # if matrices with errors are wanted
# statistics & boxplots can only be TRUE if matrixs is TRUE
statistics = FALSE # if statistics are wanted
boxplots = FALSE # if boxplots are wanted
# specify time period for statistic & boxplot analysis (not more than 3 months recommended)
# VERY IMPORTANT: Please make sure, that a forecast exists for those two dates! 
#                 If you don't know that, don't worry, the possible error will be caught anyway.
#                 Also, please make sure to enter some dates here anyway, even if you don't want statistics. 
time.start = "46_2023"
time.end = "74_2023"

# update fcst and obs data (set TRUE if desired)
if (update == TRUE) {
  res = sidfex.download.fcst()
  res2 = sidfex.download.obs()
}

# load SIDFEx index and read data
indx = sidfex.load.index()
  
# Loop over TargetIDs
for (i.tid in 1:length(tids)) {
  tid = tids[i.tid]
  
  subind.fcst = sidfex.fcst.search.extractFromTable(index = indx, tid = tid, gid = gid, mid = mid, iy = iy, idoy = idoy) # delete iy and idoy if wanted
  fcst = sidfex.read.fcst(subind.fcst)
  obs = sidfex.remaptime.obs2fcst(fcst = fcst)

  if (length(fcst$res.list) != length(obs$res.list)) { # quick sanity check
    stop("check fcst and obs, they don't have same length")
  }
  
  # global matrix pre-settings
  ncol = (10*48)+1 # only valuable for giops (10 days a 48 fcst + last day (10.000))
  nrow = length(fcst$res.list) # number of fcst available in total 
  matrix = matrix(nrow = nrow, ncol = ncol) # create empty generic matrix
  
  col = viridis(length(fcst$res.list)) # color scheme for separate gc plot
  
  ### Loop over every day in res.list
  for (i in 1:length(fcst$res.list)) {
    
    # dividing obs and fcst datasets of 10 days leadtime into 10 separate datasets for each day
    daysLeadTime = fcst$res.list[[i]]$data$DaysLeadTime # extract Leadtime
    integers = unique(floor(daysLeadTime)) # reduce to integers
    
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
      fcst_subdiv_temp = fcst 
      fcst_subdiv_temp$res.list = list(fcst_subdiv_temp$res.list[[i]])
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
      if (plot == TRUE) {
        
        plot.title = paste0(fcst$res.list[[i]]$TargetID, "_", fcst$res.list[[i]]$GroupID, "_", fcst$res.list[[i]]$MethodID, "_", fcst$res.list[[i]]$InitYear, "-", fcst$res.list[[i]]$InitDayOfYear, "_ld:", i.sub, ".png")
        if (!file.exists(file.path(wd, plot.title))) {
          png(filename = plot.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
        } else {
          print(paste0("File '", plot.title, "' already exists. Skipped it.")) # skip if file already exists
          next
        }
        
        par(mfrow=c(2,2))
        
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
        title(main=paste0("Analysis subdaily resolution \n ", fcst$res.list[[1]]$TargetID, "_", fcst$res.list[[1]]$GroupID, "_", fcst$res.list[[1]]$MethodID, "_", fcst$res.list[[1]]$InitYear, "-", fcst$res.list[[i]]$InitDayOfYear), outer = T) 
        
        dev.off()
        
        # control plots for gc dist. all errors over entire time for single lead time days (i.sub needs to be set to 1:1){# instead of 1:length... for every day separately) 
        # if(i == 1 && i.sub >= 1 && i.sub <= 10) {
        #   plot(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist, col=col[i], xlab="days lead time", ylab = "great-circle distance / m",
        #        xlim=c(i.sub-1.05,i.sub+0.05), ylim = c(range_gc[1], range_gc[2]+8000), type = "l", main = paste0("Analysis subdaily resolution \n ", fcst$res.list[[1]]$TargetID, "_", fcst$res.list[[1]]$GroupID, "_", fcst$res.list[[1]]$MethodID, "_", fcst$res.list[[1]]$InitYear, "-", idoy[1], ":", idoy[length(idoy)]))
        #   abline(h=0,v=i.sub,col="grey",lty=3)
        #   legend("topleft", c( "Error betw. obs & hh fcst"), lty = 1, col = "black", bty = "n", cex = 0.7)
        # } else {
        #   lines(x=fcst.adj$res.list[[1]]$data$DaysLeadTime, y=fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist, col=col[i])
        # }
      } # end if-clause for plotting
      
      # some matrices - calc 6 matrices (3 for gc_dist, speed & angle; twice for linear & normal (half-hourly))
      if(matrixs == TRUE) {
        
        if(i == 1 && i.sub == 1) {
          matr_dis = matrix
          matr_dis_lin = matrix
          matr_speed = matrix
          matr_speed_lin = matrix
          matr_angle = matrix
          matr_angle_lin = matrix
        }
        # matrices for great-circle distance
        matr_dis = matrix_calc(matr_dis, fcst.adj.eval$res.list[[1]]$ens.mean.gc.dist) # half-hourly
        matr_dis_lin = matrix_calc(matr_dis_lin, fcst.lin.eval$res.list[[1]]$ens.mean.gc.dist) # linear
        # removing outliners
        if(tid == "900128" && time.start == "80_2023" && time.end =="112_2023") {
          matr_dis[205,] = NA # remove day 2023-107
          matr_dis[206,] = NA # remove day 2023-108
        }
        if(tid == "900126" && time.start == "80_2023" && time.end =="112_2023") {
          matr_dis[201,] = NA # remove day 2023-110
        }
        
        # matrices for relative speed
        matr_speed = matrix_calc(matr_speed, fcst.adj.eval$res.list[[1]]$ens.mean.relspeed) # half-hourly
        matr_speed_lin = matrix_calc(matr_speed_lin, fcst.lin.eval$res.list[[1]]$ens.mean.relspeed) # linear
        # removing outliners
        if(tid == "900128" && time.start == "80_2023" && time.end =="112_2023") {
          matr_speed[205,] = NA # remove day 2023-107
          matr_speed[206,] = NA # remove day 2023-108
        }
        if(tid == "900126" && time.start == "80_2023" && time.end =="112_2023") {
          matr_speed[201,] = NA # remove day 2023-110
        }
        
        # matrices for angle
        matr_angle = matrix_calc(matr_angle, fcst.adj.eval$res.list[[1]]$ens.mean.angle) # half-hourly
        matr_angle_lin = matrix_calc(matr_angle_lin, fcst.lin.eval$res.list[[1]]$ens.mean.angle) # linear
        # removing outliners
        if(tid == "900128" && time.start == "80_2023" && time.end =="112_2023") {
          matr_angle[205,] = NA # remove day 2023-107
          matr_angle[206,] = NA # remove day 2023-108
        }
        if(tid == "900126" && time.start == "80_2023" && time.end =="112_2023") {
          matr_angle[201,] = NA # remove day 2023-110
        }
      } # end if matrix calc
      
    } # end for-loop sub dataset
    
  } # end for-loop days in res.list    
  

    
  # some statistics
  if(matrixs == TRUE && statistics == TRUE){
    
    # create reference frame to know at which index position which date is located in matrices
    reference.frame = data.frame(index = 1:length(fcst$res.list), y_doy = NA)
    for(i.doy in 1:length(fcst$res.list)){
      reference.frame$y_doy[i.doy] = paste0(fcst$res.list[[i.doy]]$InitDayOfYear, "_", fcst$res.list[[i.doy]]$InitYear)
    }
    ref.date.start = which(reference.frame$y_doy == time.start)
    if(identical(ref.date.start, integer(0))) {
      stop("No forecast exists for that starting date. Please check '~/SIDFEx/data/fcst/eccc001' for 
         'ls eccc001_giops_tid_20xx*' (add your tid & year instead of xx) to see if a forecast exists.")
      print(tid)
    }
    ref.date.end = which(reference.frame$y_doy == time.end)
    if(identical(ref.date.end, integer(0))) {
      stop("No forecast exists for that ending date. Please check '~/SIDFEx/data/fcst/eccc001' for 
         'ls eccc001_giops_tid_20xx*' (add your tid & year instead of xx) to see if a forecast exists.")
      print(tid)
    } 
    time.sel = c(ref.date.start:ref.date.end) 
    
    # mean value for every step of the highly resolved leadtime
    colmeans_dis = colMeans(matr_dis[time.sel,], na.rm = T)
    colmeans_dis_lin = colMeans(matr_dis_lin[time.sel,], na.rm = T)
    colmeans_speed = colMeans(matr_speed[time.sel,], na.rm = T)
    colmeans_speed_lin = colMeans(matr_speed_lin[time.sel,], na.rm = T)
    colmeans_angle = colMeans(matr_angle[time.sel,], na.rm = T)
    colmeans_angle_lin = colMeans(matr_angle_lin[time.sel,], na.rm = T)
    
    # One mean value for every day of the leadtime 
    # taking only the values of colmeans
    dis_1d_mean = calc.mean(colmeans_dis)
    dis_lin_1d_mean = calc.mean(colmeans_dis_lin)
    speed_1d_mean = calc.mean(colmeans_speed)
    speed_lin_1d_mean = calc.mean(colmeans_speed_lin)
    angle_1d_mean = calc.mean(colmeans_angle)
    angle_lin_1d_mean = calc.mean(colmeans_angle_lin)
    
    # standard deviation
    sd_dis = colSds(matr_dis[time.sel,], na.rm = T)
    sd_disP = colmeans_dis + sd_dis
    sd_disM = colmeans_dis - sd_dis
    sd_dis_lin = colSds(matr_dis_lin[time.sel,], na.rm = T)
    sd_disP_lin = colmeans_dis_lin + sd_dis_lin
    sd_disM_lin = colmeans_dis_lin - sd_dis_lin
    
    sd_speed = colSds(matr_speed[time.sel,], na.rm = T)
    sd_speedP = colmeans_speed + sd_speed
    sd_speedM = colmeans_speed - sd_speed
    sd_speed_lin = colSds(matr_speed_lin[time.sel,], na.rm = T)
    sd_speedP_lin = colmeans_speed_lin + sd_speed_lin
    sd_speedM_lin = colmeans_speed_lin - sd_speed_lin
    
    sd_angle = colSds(matr_angle[time.sel,], na.rm = T)
    sd_angleP = colmeans_angle + sd_angle
    sd_angleM = colmeans_angle - sd_angle
    sd_angle_lin = colSds(matr_angle_lin[time.sel,], na.rm = T)
    sd_angleP_lin = colmeans_angle_lin + sd_angle_lin
    sd_angleM_lin = colmeans_angle_lin - sd_angle_lin

    # some statistical plotting
    # function for plotting
    stat.plot = function(x, high.res, lin, sdP, sdM, sdP_lin, sdM_lin, ylab, title, pos.lg) {
      ylim1 = range(range_calc(high.res), range_calc(lin), range_calc(sdP), range_calc(sdM))
      plot(x = x, y = sdM, xlab="days lead time", ylab = ylab,
          main = title, ylim = ylim1, col = "lightgrey", type = "l", lty = "dashed")
      lines(x = x, y = sdP, col = "lightgrey", lty = "dashed")
      lines(x = x, y = sdP_lin, col = "lightblue", lty = "dashed")
      lines(x = x, y = sdM_lin, col = "lightblue", lty = "dashed")
      # polygon(c(x, rev(x)), c(sdP, rev(sdM)), col = "lightyellow", border = NA)
      abline(h = 0,v = c(1:10),col = "grey",lty = 3)
      lines(x = x, y = high.res, col = "black")
      lines(x = x, y = lin, col = "blue")
      legend(pos.lg, c("High, subdaily resolution", "Daily resolution", "Standard Deviation (High res)", "Standard Deviation (Daily res)"), col = c("black", "blue", "grey", "lightblue"), lty = c(1,1,2,2), bty = "n", cex = 0.7)
    }
    
    # plots for half-hourly steps
    # different x-axis options, depending on wanted Leadtime
    x = fcst$res.list[[1]]$data$DaysLeadTime # entire leadtime
    x_1d = fcst$res.list[[1]]$data$DaysLeadTime[1:49] # only day 1
    x_2to10 = fcst$res.list[[1]]$data$DaysLeadTime[-(1:48)] # everything except day 1
    
    year = strsplit(time.start, '_')[[1]][2] # extract year
    first = as.numeric(strsplit(time.start, '_')[[1]][1]) # extract first day, remove 0 if there
    last = as.numeric(strsplit(time.end, '_')[[1]][1]) # extract last day, remove 0 if there
    date = paste0(year, ":", first, "-", last) # total selected time period
    
    myDir = unlist(strsplit(getwd(), '/'))
    stat.Dir = paste0(myDir[-length(myDir)], collapse = '/') # directory for statistical plots
    wd2 = setwd(stat.Dir) # change directory
    wd2 = setwd(stat.Dir) # change directory
    print(paste0("Directory has been changed to ", wd2))
    
    dist.title = paste0(tid, "_", date, "_gc.dist.png")
    speed.title = paste0(tid, "_", date, "_speed.png")
    angle.title = paste0(tid, "_", date, "_angle.png")
    daily.title = paste0(tid, "_", date, "_daily.png")

    # gc-dist
    if (!file.exists(file.path(wd2, dist.title))) {
      png(filename = dist.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
    } else {
      print(paste0("File '", dist.title, "' already exists. Skipped it.")) # skip if file already exists
    }
    title_gc1 = paste0("Great-circle distance - ", tid, " - ", date)
    ylab_gc1 = "mean error / m"
    stat.plot(x, colmeans_dis, colmeans_dis_lin, sd_disP, sd_disM, sd_disP_lin, sd_disM_lin, ylab_gc1, title_gc1, "topleft")
    dev.off()
    
    # speed (entire leadtime)
    if (!file.exists(file.path(wd2, speed.title))) {
      png(filename = speed.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
    } else {
      print(paste0("File '", speed.title, "' already exists. Skipped it.")) # skip if file already exists
    }
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) # layout for 3 plots on one page
    
    title_sp1 = paste0("Relative Speed - ", tid, " - ", date)
    ylab_sp1 = "mean error"
    stat.plot(x, colmeans_speed, colmeans_speed_lin, sd_speedP, sd_speedM, sd_speedP_lin, sd_speedM_lin, ylab_sp1, title_sp1, "topright")
    # leadtime 1
    title_sp2 = paste0("Relative Speed - ", tid," - ", date, "; Leadtime 0-1")
    stat.plot(x_1d, colmeans_speed[1:49], colmeans_speed_lin[1:49], sd_speedP[1:49], sd_speedM[1:49], sd_speedP_lin[1:49], sd_speedM_lin[1:49], ylab_sp1, title_sp2, "topright")
    # leadtime 2:10
    title_sp3 = paste0("Relative Speed - ", tid, " - ", date, "; Leadtime 1-10")
    stat.plot(x_2to10, colmeans_speed[49:(length(colmeans_speed))], colmeans_speed_lin[49:(length(colmeans_speed_lin))], sd_speedP[49:(length(colmeans_speed_lin))], sd_speedM[49:(length(colmeans_speed_lin))], sd_speedP_lin[49:(length(colmeans_speed))], sd_speedM_lin[49:(length(colmeans_speed))], ylab_sp1, title_sp3, "topleft")
    dev.off()
    
    # angle (entire leadtime)
    if (!file.exists(file.path(wd2, angle.title))) {
      png(filename = angle.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
    } else {
      print(paste0("File '", angle.title, "' already exists. Skipped it.")) # skip if file already exists
    }
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))# layout for 3 plots on one page
    
    title_an1 = paste0("Relative Angle - ", tid, " - ", date)
    ylab_an1 = "mean error / degree left"
    stat.plot(x, colmeans_angle, colmeans_angle_lin, sd_angleP, sd_angleM, sd_angleP_lin, sd_angleM_lin, ylab_an1, title_an1, "topright")
    # leadtime 1
    title_an2 = paste0("Relative Angle - ", tid, " - ", date, "; Leadtime 0-1")
    stat.plot(x_1d, colmeans_angle[1:49], colmeans_angle_lin[1:49], sd_angleP[1:49], sd_angleM[1:49], sd_angleP_lin[1:49], sd_angleM_lin[1:49], ylab_an1, title_an2, "topright")
    # leadtime 2:10
    title_an3 = paste0("Relative Angle - ", tid, " - ", date, "; Leadtime 1-10")
    stat.plot(x_2to10, colmeans_angle[49:(length(colmeans_angle))], colmeans_angle_lin[49:(length(colmeans_angle_lin))], sd_angleP[49:(length(colmeans_angle_lin))], sd_angleM[49:(length(colmeans_angle_lin))], sd_angleP_lin[49:(length(colmeans_angle_lin))], sd_angleM_lin[49:(length(colmeans_angle_lin))], ylab_an1, title_an3, "topleft")
    dev.off()
    
    # plots for daily steps (same as above, only with type points instead of lines)
    stat.points = function(x, high.res, lin, ylab, title, pos.lg) {
      ylim1 = range(range_calc(high.res), range_calc(lin))
      plot(x = x, y = high.res, xlab="days lead time", ylab = ylab,
           main = title, ylim = ylim1, col = "black", type = "p")
      abline(h = 0,v = c(1:10),col = "grey",lty = 3)
      points(x = x, y=lin, col="blue")
      legend(pos.lg, c("High, subdaily resolution", "Daily resolution"), col = c("black", "blue"), lty = c(1,1), bty = "n", cex = 0.7)
    }
    
    if (!file.exists(file.path(wd2, daily.title))) {
      png(filename = daily.title, height = 21, width = 29.7, units = "cm", res = 300) # save in DinA4 format (open png)
    } else {
      print(paste0("File '", daily.title, "' already exists. Skipped it.")) # skip if file already exists
    }
    par(mfrow=c(3,1))
    stat.points(1:10, dis_1d_mean, dis_lin_1d_mean, ylab_gc1, title_gc1, "topleft")
    stat.points(1:10, speed_1d_mean, speed_lin_1d_mean, ylab_sp1, title_sp1, "topleft")
    stat.points(1:10, angle_1d_mean, angle_lin_1d_mean, ylab_an1, title_an1, "topright")
    dev.off()
    
    setwd(wd)
    print((paste0("Directory has been changed back to ", wd)))
  } # end if-clause for statistics
  
  if(matrixs == TRUE && statistics == TRUE && boxplots == TRUE) {
    
    myDir = unlist(strsplit(getwd(), '/'))
    stat.Dir = paste0(myDir[-length(myDir)], collapse = '/') # directory for statistical plots
    wd3 = setwd(stat.Dir) # change directory
    wd3 = setwd(stat.Dir) # change directory
    print(paste0("Directory has been changed to ", wd3))
    save_dir2 = date
    if (!file.exists(save_dir2)) { # new folder is only being created if doesn't exist already
      dir.create(save_dir2)
      print(paste0("new directory created:", wd3, "/", save_dir2))
    }
    wd4 = setwd(paste0(wd3, "/", save_dir2 )) # set new wd 
    wd4 = setwd(paste0(wd3, "/", save_dir2 )) # set new wd
    print(paste0("Directory has been changed to ", wd4))
    
    # getting daily data sets for each parameter
    dis_perDay = split.daily(matr_dis[time.sel,])
    dis_lin_perDay = split.daily(matr_dis_lin[time.sel,])
    speed_perDay = split.daily(matr_speed[time.sel,])
    speed_lin_perDay = split.daily(matr_speed_lin[time.sel,])
    angle_perDay = split.daily(matr_angle[time.sel,])
    angle_lin_perDay = split.daily(matr_angle_lin[time.sel,])
    
    # some boxplotting
    for(i.bplot in 1:length(dis_perDay)) {
      bxp.dis = paste0("bxp_", tid, "_", date, "_ld:", i.bplot-1, "-", i.bplot, "_gc.dist.png") # titles for saving plots
      bxp.speed = paste0("bxp_", tid, "_", date, "_ld:", i.bplot-1, "-", i.bplot, "_speed.png")
      bxp.angle = paste0("bxp_", tid, "_", date, "_ld:", i.bplot-1, "-", i.bplot, "_angle.png")
        
      # gc-dist boxplot
      if (!file.exists(file.path(wd4, bxp.dis))) {
        png(filename = bxp.dis, height = 29.7, width = 21, units = "cm", res = 300) # save in DinA4 format (open png)
        par(mfrow = c(2,1))
        boxplot(dis_perDay[[i.bplot]], 
              main = paste0("Great-circle distance - ", tid, " - ", date),
              xlab="days lead time", 
              ylab = "great-circle distance / m"
              )
        grid(lty = "dotted")
        boxplot(dis_lin_perDay[[i.bplot]], 
                main = paste0("Great-circle distance (linear) - ", tid, " - ", date),
                xlab="days lead time", 
                ylab = "great-circle distance / m"
        )
        grid(lty = "dotted")
        dev.off()
      } else {
        print(paste0("File '", bxp.dis, "' already exists. Skipped it.")) # skip if file already exists
      }
      
      # speed boxplot
      if (!file.exists(file.path(wd4, bxp.speed))) {
        png(filename = bxp.speed, height = 29.7, width = 21, units = "cm", res = 300) # save in DinA4 format (open png)
        par(mfrow = c(2,1))
        boxplot(speed_perDay[[i.bplot]], 
                main = paste0("Relative Speed - ", tid, " - ", date),
                xlab="days lead time",
                ylab = "relative speed"
        )
        grid(lty = "dotted")
        boxplot(speed_lin_perDay[[i.bplot]], 
                main = paste0("Relative Speed (linear) - ", tid, " - ", date),
                xlab="days lead time",
                ylab = "relative speed"
        )
        grid(lty = "dotted")
        dev.off()
      } else {
        print(paste0("File '", bxp.speed, "' already exists. Skipped it.")) # skip if file already exists
      }
      
      # angle boxplot
      if (!file.exists(file.path(wd4, bxp.angle))) {
        png(filename = bxp.angle, height = 29.7, width = 21, units = "cm", res = 300) # save in DinA4 format (open png)
        par(mfrow = c(2,1))
      boxplot(angle_perDay[[i.bplot]], 
              main = paste0("Relative Angle - ", tid, " - ", date),
              xlab="days lead time", 
              ylab = "relative angle / degree left"
      )
      grid(lty = "dotted")
      boxplot(angle_lin_perDay[[i.bplot]], 
              main = paste0("Relative Angle (linear) - ", tid, " - ", date),
              xlab="days lead time", 
              ylab = "relative angle / degree left"
      )
      grid(lty = "dotted")
      dev.off()
      } else {
        print(paste0("File '", bxp.angle, "' already exists. Skipped it.")) # skip if file already exists
      }
      
    } # end plotting loop
    
    setwd(wd)
    print((paste0("Directory has been changed back to ", wd)))
  } # end if-clause for boxplots

} # end for-loop TargetID
