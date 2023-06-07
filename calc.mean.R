# function to calculate mean daily error (calcMean = T) & to divide matrix into subsets for each day of the leadtime (calcMean = F)

# unneccesarily complex way of trying to figure of full day of "numeric"
# fullday = which(colmeans_dis == 0)
# lastDay = fullday[(length(fullday))] # extract last day to attach later
# fullday = fullday[-(length(fullday))] # remove last day
# next_value = fullday + 1 # real full days, !! don't get irritated: index is being manipulated!! Not aligning with real data, only ugly work aroung
# fullday = c(colmeans_dis[1], next_value, lastDay) # bringing it all together

fullday = c(1, 49, 97, 145, 193, 241, 289, 337, 385, 433, 481) # manually setting full days

day_mean = c() # initialize empty vector to store data later

calc.mean = function(matrx) {
  matrx = as.data.frame(matrx) # bring into data.frame format to access via indices
  for(i.day in 1:(length(fullday)-1)) { 
    if(i.day < (length(fullday)-1)) {
      df_1d = as.data.frame(matrx[fullday[i.day]:(fullday[i.day+1]-1),])
    } else { 
      df_1d = as.data.frame(matrx[fullday[i.day]:fullday[i.day+1],]) # add 10th value (for mean calc)
    } 
    colnames(df_1d) = paste0("day_", (i.day-1), "-", (i.day)) # adjust colnames to days
    daily_mean = colMeans(df_1d) # calculate col mean
    day_mean = c(day_mean, daily_mean) # fill in vector with mean values
  }
  return(day_mean)
}