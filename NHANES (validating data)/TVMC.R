# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}



get_next_candidates <- function(sequences, time, activity){
  this <- sequences[,c(time, time+1)]
  this <- this[this[,1] == activity,]
  return(as.vector(this[,2]))
}



# Input sequences have to be vectorized already.
# If not, use vectorize_seq function before using this function
simulate_one_sequence <- function(sequences, cluster, margin = 60){
  len <- ncol(sequences)
  result <- rep(NA, len)

  # which cluster does this sequence that we are simulating belong to?
  cluster_assignment <- sample(cluster, 1)

  row_num <- which(cluster == cluster_assignment)
  this_cluster <- sequences %>%
    slice(row_num)
  current_activity <- sample(this_cluster[,1], 1)
  current_time <- 1
  
  result[1] <- current_activity
  
  while (current_time < len){
    candidates <- get_next_candidates(this_cluster, current_time, current_activity)
    current_time <- current_time + 1
    current_activity <- sample(candidates, 1)
    result[current_time] <- current_activity
  }

  return(result)
}


# simulate_multiple_sequences <- function(data, cluster, n){
#   simulated <- replicate(n, simulate_one_sequence(data, cluster))
#   result <- data.frame(t(simulated))
#   return(result)
# }


parallel_simulate_multiple_sequences <- function(data, cluster, n, seeds, margin = 60){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'MASS', 'clValid')) %dopar% {
    source("TVMC.R")
    set.seed(seeds[i])
    simulated <- simulate_one_sequence(data, cluster, margin) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}






