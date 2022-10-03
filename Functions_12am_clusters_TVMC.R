# This version of Functions_start_with_clusters starts at 12am!
# Need to create the distance matrix first at once, which can be computationally intensive.



# When the user does not input their own clustering assignment
# create_clusters <- function(distances, clust = "hclust", method = "average"){
  # if (clust != "hclust"){
  #   warning("Other clustering methods are not supported. Please input your own clustering assignment or default to hierarchical clustering.")
  # }

#   hclust_real <- stats::hclust(as.dist(distances), method = method)
#   dunn_index <- rep(0, 100) # The first one is Inf; we will denote as 0.
#   for (i in 2:100){
#     fit_real <- cutree(hclust_real, k = i)
#     dunn_index[i] <- dunn(distance = distances, clusters = fit_real)
#   }
# 
#   num <- which.max(dunn_index)
#   final_cluster <- cutree(hclust_real, k = num)
#   return(final_cluster) # need to suppress the warning
# }

  

# When the user does not input their own clustering assignment
create_clusters_Dunn <- function(distances, clust = "hclust", method = "average"){
  # if (clust != "hclust"){
  #   warning("Other clustering methods are not supported. Please input your own clustering assignment or default to hierarchical clustering.")
  # }
  
  hclust_real <- stats::hclust(as.dist(distances), method = method)
  dunn_index <- matrix(0, nrow = 100, ncol = 8)
  
  for (i in 2:100){
    fit_real <- cutree(hclust_real, k = i)
    scatt <- cls.scatt.diss.mx(distances, fit_real)
    dunn_index[i,] <- clv.Dunn(scatt, intracls = c("complete", "average"),
                               intercls = c("single", "complete", "average", "hausdorff")) 
  }

  num <- apply(dunn_index, 2, which.max)
  final_cluster <- cutree(hclust_real, k = num[4])
  return(final_cluster)
}




# create_clusters_DB <- function(distances, clust = "hclust", method = "average"){
#   # if (clust != "hclust"){
#   #   warning("Other clustering methods are not supported. Please input your own clustering assignment or default to hierarchical clustering.")
#   # }
#   
#   hclust_real <- stats::hclust(as.dist(distances), method = method)
#   DB_index <- matrix(10, nrow = 100, ncol = 8)
#   
#   for (i in 2:100){
#     fit_real <- cutree(hclust_real, k = i)
#     scatt <- cls.scatt.diss.mx(distances, fit_real)
#     DB_index[i,] <- clv.Davies.Bouldin(scatt, intracls = c("complete", "average"),
#                                        intercls = c("single", "complete", "average", "hausdorff")) 
#   }
#   
#   num <- apply(DB_index, 2, which.min)
#   print(num)
#   final_cluster <- cutree(hclust_real, k = num[1])
#   return(final_cluster)
# }



  
# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}



# sequence: sequence for one day.
# activity: character
# time: column index from 1 to 1440 (not in xx:xx format)
# return a vector of three numbers (start time, end time, length of the subsequence)
# find_subsequence <- function(sequence, activity, time){
#   rle_oneday <- rle(sequence == activity)$length
#   cumSum <- c(0, cumsum(rle_oneday))
#   cutted <- cut(time, cumSum)
#   levels(cutted) <- 1:length(rle_oneday)
#   start <- cumSum[as.numeric(cutted)] + 1
#   end <- cumSum[as.numeric(cutted) + 1]
#   return(c(start, end, end - start + 1))
# }



# end_of_first_activity <- function(sequence){
#   return((rle(sequence)$length)[[1]])
# }



get_next <- function(sequences, time, activity){
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
  
  while (current_time < 1440){
    candidates <- get_next(this_cluster, current_time, current_activity)
    current_time <- current_time + 1
    current_activity <- sample(candidates, 1)
    result[current_time] <- current_activity
  }

  return(result)
}




simulate_multiple_sequences <- function(data, cluster, n){
  simulated <- replicate(n, simulate_one_sequence(data, cluster))
  result <- data.frame(t(simulated))
  return(result)
}


parallel_simulate_multiple_sequences <- function(data, cluster, n){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'MASS', 'clValid')) %dopar% {
    source("Functions_12am_clusters_TVMC.R")
    set.seed(i)
    simulated <- simulate_one_sequence(data, cluster) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}






