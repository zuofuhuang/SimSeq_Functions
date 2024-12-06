get_next_candidates <- function(sequences, time, activity){
  this <- sequences[,c(time, time+1)]
  this <- this[this[,1] == activity,]
  # return(as.vector(this[,2]))
  result <- as.vector(this[,2])
  if (length(result) == 0){
    return("N/A")
  } else {
    return(result)
  }
}


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


# set.seed(475)
# imputed <- impute(sequences, 60)
# imputed_seqs <- cbind(sequences, imputed)
# write.csv(imputed_seqs, "imputed_seqs.csv")
# 
# 
# set.seed(477)
# imputed_120 <- impute(sequences, 120)
# imputed_seqs_120 <- cbind(sequences, imputed_120)
# write.csv(imputed_seqs_120, "imputed_seqs_120.csv")


# set.seed(465)
# imputed <- impute(sequences_work, 60)
# imputed_seqs_work <- cbind(sequences_work, imputed)
# write.csv(imputed_seqs_work, "imputed_seqs_work.csv")

# set.seed(480)
# imputed <- impute(sequences_work_v2, 60) # v2 means that workdays have at least 30 min of work
# imputed_seqs_work_v2 <- cbind(sequences_work_v2, imputed)
# write.csv(imputed_seqs_work_v2, "imputed_seqs_work_v2.csv")


# set.seed(460)
# imputed2 <- impute(sequences, 60)
# imputed_seqs_s4k5 <- cbind(sequences, imputed)
# write.csv(imputed_seqs_s4k5, "imputed_seqs_s4k5.csv")


set.seed(460)
imputed2 <- impute(accel_subset, 60)
imputed_seqs_s4k5_updated <- cbind(accel_subset, imputed2)
write.csv(imputed_seqs_s4k5_updated, "imputed_seqs_s4k5_updated.csv")
