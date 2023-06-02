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