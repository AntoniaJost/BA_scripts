# function to create matrices with evaluation (error) data

matrix_calc = function(matrix, input.data) {
  if(ncol == length(fcst$res.list[[1]]$data$DaysLeadTime)) {
    colnames(matrix) = round(fcst$res.list[[1]]$data$DaysLeadTime, digits = 3)
  } else{
    stop("ncol not equal to total DaysLeadTime")
  }
  rownames(matrix)[i:nrow] = paste0(fcst.adj$res.list[[1]]$InitYear, "_", fcst.adj$res.list[[1]]$InitDayOfYear)
  cols = round(fcst.adj$res.list[[1]]$data$DaysLeadTime, digits = 3)
  matrix[i,as.character(cols)] = input.data 
  return(matrix)
}