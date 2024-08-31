# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}


# 5 hidden states "Home", "Work", "Activity", "Transit", "Others"



# Input sequences have to be vectorized already.
# If not, use vectorize_seq function before using this function
simulate_one_sequence <- function(sequences, cluster, margin = 60){
  
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






