# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}






get_next_activity <- function(sequences, time, activity, margin){
  this <- sequences[,max(1,time-margin):min(1440,time+margin)]
  new <- as.matrix(this[,1:2])
  for (i in 2:(ncol(this) - 1)){
    new <- rbind(new, as.matrix(this[,c(i,i+1)]))
  }
  transMat <- prop.table(table(new[,1], new[,2]), 1)
  
  index <- which(rownames(transMat) == activity)
  next_act <- sample(1:ncol(transMat), size = 1, prob = transMat[index,])
  return(colnames(transMat)[next_act])
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
    current_time <- current_time + 1
    current_activity <- get_next_activity(this_cluster, current_time, current_activity, margin = margin)
    result[current_time] <- current_activity
  }

  return(result)
}


# Currently, don't use the non-parallelized version anymore. 
# To support it later, would modify the code and add seed argument etc.
# simulate_multiple_sequences <- function(data, cluster, n){
#   simulated <- replicate(n, simulate_one_sequence(data, cluster))
#   result <- data.frame(t(simulated))
#   return(result)
# }


parallel_simulate_multiple_sequences <- function(data, cluster, n, seeds, margin = 60){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)

  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'MASS', 'clValid',"stringr")) %dopar% {
    source("TVMC_neighbor.R")
    set.seed(seeds[i])
    simulated <- simulate_one_sequence(data, cluster, margin) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}






