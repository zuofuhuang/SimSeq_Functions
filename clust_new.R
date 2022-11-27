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




create_clusters_DB <- function(distances, clust = "hclust", method = "average"){
  # if (clust != "hclust"){
  #   warning("Other clustering methods are not supported. Please input your own clustering assignment or default to hierarchical clustering.")
  # }
  
  hclust_real <- stats::hclust(as.dist(distances), method = method)
  DB_index <- matrix(10, nrow = 100, ncol = 8)
  
  for (i in 2:100){
    fit_real <- cutree(hclust_real, k = i)
    scatt <- cls.scatt.diss.mx(distances, fit_real)
    DB_index[i,] <- clv.Davies.Bouldin(scatt, intracls = c("complete", "average"),
                                       intercls = c("single", "complete", "average", "hausdorff")) 
  }
  
  num <- apply(DB_index, 2, which.min)
  print(num)
  final_cluster <- cutree(hclust_real, k = num[1])
  return(final_cluster)
}



  
# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}



# sequence: sequence for one day.
# activity: character
# time: column index from 1 to 1440 (not in xx:xx format)
# return a vector of three numbers (start time, end time, length of the subsequence)
find_subsequence <- function(sequence, activity, time){
  rle_oneday <- rle(sequence == activity)$length
  cumSum <- c(0, cumsum(rle_oneday))
  cutted <- cut(time, cumSum)
  levels(cutted) <- 1:length(rle_oneday)
  start <- cumSum[as.numeric(cutted)] + 1
  end <- cumSum[as.numeric(cutted) + 1]
  return(c(start, end, end - start + 1))
}



end_of_first_activity <- function(sequence){
  return((rle(sequence)$length)[[1]])
}



# from: from this activity to other activities.
# On the left, we are interested in how other activities transform to this activity, so from = FALSE.
activity_transition <- function(sequence, time, activity, margin = 60, from = TRUE){
  len <- length(sequence)
  window <- c(max(1, time - margin), min(time + margin, len))
  subseq <- sequence[window[1]:window[2]] # The subsequence (the interval) to focus on
  rle_result <- rle(subseq)
  activity_positions <- which(rle_result$values == activity)
  cumSum <- cumsum(rle_result$lengths)

  if (length(activity_positions) == 0) return() # The interval does not contain the corresponding activity

  if (from){
    if(activity_positions[length(activity_positions)] == length(rle(subseq)$length)){ # get rid of the last index if it is true
      activity_positions <- head(activity_positions, -1)
    }
    position <- window[1] + cumSum[activity_positions] # Position in the entire day
    activity <- (rle_result$values)[activity_positions + 1]
  } else {
    if(activity_positions[1] == 1){  # get rid of the first index if it is true
      activity_positions <- tail(activity_positions, -1)
    }
    position <- window[1] + cumSum[activity_positions - 1] - 1 # Position in the entire day
    activity <- rle_result$values[activity_positions - 1]
  }

  if (length(activity_positions) == 0) return() # After we delete the first/last one, interval does not contain corresponding activity

  if (!setequal(unname(sequence[position]),unname(activity))){
    stop("Activity state does not match")
  }

  # Find out the length of the subsequence that contains position
  lens <- rle(sequence)$length
  seq_cumSum <- c(0,cumsum(lens))
  cutted <- cut(position, seq_cumSum)
  levels(cutted) <- 1:length(lens)
  position <- lens[cutted]
  names(position) <- paste(activity, as.numeric(cutted == length(lens)))

  return(position)
}


create_freq_table <- function(sequences){
  states <- seqstatl(sequences)
  freq_table <- matrix(NA, nrow = length(states), ncol = ncol(sequences))
  
  for (i in 1:length(states)){
    dat <- (sequences == states[i])
    sums <- colSums(dat, na.rm = T)/nrow(dat)
    freq_table[i,] <- sums
  }
  freq_table <- data.frame(freq_table)
  rownames(freq_table) <- states
  return(freq_table)
}



# Input sequences have to be vectorized already.
# If not, use vectorize_seq function before using this function
simulate_one_sequence <- function(sequences, cluster, margin = 60){
  len <- ncol(sequences)
  result <- rep(NA, len)

  # which cluster does this sequence that we are simulating belong to?
  cluster_assignment <- sample(cluster, 1)

  row_num <- which(cluster == cluster_assignment)
  sequences <- sequences %>%
    slice(row_num)
  starting_activity <- sample(sequences[,1], 1)

  subset <- sequences %>%
    filter(.[[1]] == starting_activity)
  # subsequences <- data.frame(t(apply(subset, 1, find_subsequence, activity = starting_activity, time = 1)))
  ends <- apply(subset, 1, end_of_first_activity)

  # Draw one randomly
  # right <- subsequences$X2[sample(1:nrow(subsequences),1)]
  right <- sample(ends, 1)
  result[1:right] <- starting_activity
  right_activity <- starting_activity

  while (right < len){
    transitions <- unlist(apply(sequences, 1, activity_transition,
                                time = right, activity = right_activity, from = TRUE, margin = margin))

    if (is.null(transitions)){ # This is almost impossible to happen
      result[(right + 1) : min(right + margin, len)] <- right_activity
      right <- min(right + margin, len)
    } else {
      all_names <- unlist(strsplit(names(transitions), " "))
      all_act <- all_names[c(TRUE, FALSE)]
      till_end <- as.numeric(all_names[c(FALSE, TRUE)])
      
      rand <- sample(1:length(transitions), size = 1)
      activity <- all_act[rand]

      if (right < 1320 & till_end[rand]){
        dur <- len - right
      } else {
        dur <- transitions[rand]
        if (dur + right > len)  dur <- len - right
      }

      result[(right + 1):(right + dur)] <- rep(activity, dur)
      right <- right + dur
      right_activity <- activity
    }
  }

  return(result)
}



# Currently, we don't use the non-parallelized version anymore. 
# To support it later, modify the code and add seed argument etc.
# simulate_multiple_sequences <- function(data, cluster, n, margin = 60){
#   simulated <- replicate(n, simulate_one_sequence(data, cluster, margin = margin))
#   result <- data.frame(t(simulated))
#   return(result)
# }


parallel_simulate_multiple_sequences <- function(data, cluster, n, seeds, margin = 60){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'MASS', 'clValid')) %dopar% {
    set.seed(seeds[i])
    source("Functions_12am_clusters.R")
    simulated <- simulate_one_sequence(data, cluster, margin) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}



