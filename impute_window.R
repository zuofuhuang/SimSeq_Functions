get_next_candidates <- function(sequences, time, activity){
  this <- sequences[,c(time, time+1)]
  this <- this[this[,1] == activity,]
  return(as.vector(this[,2]))
}


# One problem with this is that there are more NAs toward the end of day, 
# but it doesn't matter if we treat NA as a valid state itself
impute_one_sequence <- function(start, sequences, window){
  imputed <- c(start,rep(NA, (window-1)))
  for (i in 1:(window-1)){
    imputed[i+1] <- sample(get_next_candidates(sequences, i, imputed[i]), 1)
  }
  return(imputed)
}


impute <- function(sequences, window){
  # Need to make an assumption that the activity at 23:59 continues into 12 am the next day
  # I think the only method not requiring this assumption is LSTM
  starts <- as.vector(sequences[,ncol(sequences)])
  imputed <- t(sapply(starts, impute_one_sequence, sequences = sequences, window = window))
}


imputed <- impute(sequences, 60)

imputed_seqs <- cbind(sequences, imputed)
