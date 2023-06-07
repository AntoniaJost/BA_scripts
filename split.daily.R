# function to calculate mean daily error (calcMean = T) & to divide matrix into subsets for each day of the leadtime (calcMean = F)

# unneccesarily complex way of trying to figure of full day of "numeric"
# fullday = which(colmeans_dis == 0)
# lastDay = fullday[(length(fullday))] # extract last day to attach later
# fullday = fullday[-(length(fullday))] # remove last day
# next_value = fullday + 1 # real full days, !! don't get irritated: index is being manipulated!! Not aligning with real data, only ugly work aroung
# fullday = c(colmeans_dis[1], next_value, lastDay) # bringing it all together

fullday = c(1, 49, 97, 145, 193, 241, 289, 337, 385, 433, 481) # manually setting full days
sep.day = list() # initialize empty vector to store data later

split.daily = function(matrx) {
  for(i.spd in 1:(length(fullday)-1)) { 
    if(i.spd < (length(fullday)-1)) {
      df_1day = as.data.frame(matrx[,fullday[i.spd]:(fullday[i.spd+1]-1)])
    } else { 
      df_1day = as.data.frame(matrx[,fullday[i.spd]:fullday[i.spd+1]]) # add 10th value (for mean calc)
    } 
    sep.day[[i.spd]] = df_1day # fill in list with daily subsets 
  }
  return(sep.day)
}