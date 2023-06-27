## TEMPERATURE ##
## Load modules
rm(list=ls())
require(SIDFEx)
require(spheRlab)
require(viridis)
tid = "900128"
dir = "/home/anjost001/Documents/AWI/Bachelorarbeit/data_analysis/SIC/"
obs.remap = sidfex.read.obs(TargetID = tid, data.path = dir)

plot.dir = dir
options(repr.plot.width=16, repr.plot.height=14)
bc = sl.boundingcircle(lon = obs.remap$data$Lon,lat = obs.remap$data$Lat)
pir = sl.plot.init(projection = "polar",polar.latbound = 90-2*bc$radius, polar.lonlatrot = c(bc$center_lon, bc$center_lat, 0), file.name = file.path(plot.dir,paste0("map_",tid,"_sic.pdf")),device = "pdf",do.init = T,)
sl.plot.naturalearth(pir, what = "land", fill.col = "grey", resolution = "medium")
sl.plot.lonlatgrid(pir, labels = T)
cb.brks =seq(0,100,length.out = 20)
cb = sl.colbar(cols=viridis(length(cb.brks)),N=length(cb.brks)+1)
sic_num2cb = sl.num2colbar(obs.remap$data$sic_osisaf,breaks = cb.brks,colbar = cb)
#sl.plot.points(pir,lon = obs.remap$data$Lon,lat = obs.remap$data$Lat,pch=16,)
#sl.plot.lines(pir,lon = obs.remap$data$Lon,lat = obs.remap$data$Lat,lty = 1)
sl.plot.points(pir,lon = obs.remap$data$Lon,lat = obs.remap$data$Lat,col=unlist(cb[sic_num2cb$colour.index]),pch=16)
sl.plot.colbar(sic_num2cb,do.init=F,xshift=0.28*(pir$xlim[2]-pir$xlim[1]),yshift=0*(pir$ylim[2]-pir$ylim[1]),len = 0.2*(pir$ylim[2]-pir$ylim[1]),triag.ends = T,labels.cex = 1.2,units = "SIC")
sl.plot.end(pir)
