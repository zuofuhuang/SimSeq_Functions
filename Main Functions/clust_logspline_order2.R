sample_one <- function(samp){
  if (length(samp) == 1){
    return(samp)
  } else {
    return(sample(samp, 1))
  }
}

sample_logspline <- function(potential_durations){
  # logspline needs at least 10 data points and 3 unique ones.
  if (length(potential_durations) > 9 & length(unique(potential_durations)) > 2){ 
    dur <- tryCatch({
      max(1,round(rlogspline(1, logspline(potential_durations, lbound = 1, ubound = 1440))))
    },
    warning = function(w){
      sample(potential_durations, 1)
    }, 
    error = function(e){
      sample(potential_durations, 1)
    })
  } else {
    dur <- sample(potential_durations, 1)
  }
  return(dur)
}



# Takes in a row of TraMineR sequence
vectorize_seq <- function(sequence){
  return(as.vector(unlist(sequence)))
}


end_of_first_activity <- function(sequence){
  return((rle(sequence)$length)[[1]])
}


# from: from this activity to other activities.
# On the left, we are interested in how other activities transform to this activity, so from = FALSE.
activity_transition_order2 <- function(sequence, time, activity, margin = 60, from = TRUE){
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
    next_activity <- unlist((rle_result$values)[activity_positions + 1])
  } else {
    if(activity_positions[1] == 1){  # get rid of the first index if it is true
      activity_positions <- tail(activity_positions, -1)
    }
    position <- window[1] + cumSum[activity_positions - 1] - 1 # Position in the entire day
    next_activity <- unlist(rle_result$values[activity_positions - 1])
  }
  
  if (length(activity_positions) == 0) return() # After we delete the first/last one, interval does not contain corresponding activity

  if (!setequal(unlist(sequence[position]),unname(next_activity))){
    stop("Activity state does not match")
  }
  
  # Find out the length of the subsequence that contains position
  lens <- rle(sequence)$length
  all_activities <- unlist(rle(sequence)$values)
  seq_cumSum <- c(0,cumsum(lens))
  cutted <- cut(position, seq_cumSum)
  levels(cutted) <- 1:length(lens)
  
  index_prior <- as.numeric(cutted) - 2    # Index of the previous activities by two
  names <- rep("", length(index_prior))
  for(i in 1:length(index_prior)){
    if (index_prior[i] == 0){
      names[i] <- paste0("None-", activity, "-", next_activity[i])
    }
    if (index_prior[i] != 0){
      names[i] <- paste0(all_activities[index_prior[i]], "-", activity, "-", next_activity[i])
    }
  }
  
  durations <- lens[cutted]
  names(durations) <- names
  
  return(durations)
}



# Risks overfitting
simulate_one_sequence_order2 <- function(sequences, cluster, margin = 60){
  len <- ncol(sequences)
  result <- rep(NA, len)
  
  # which cluster does this sequence that we are simulating belong to
  # could use sample_one, if some clusters have size 1
  cluster_assignment <- sample(cluster, 1)
  
  row_num <- which(cluster == cluster_assignment)
  sequences <- sequences %>%
    slice(row_num)
  starting_activity <- sample(sequences[,1], 1)
  
  subset <- sequences %>%
    filter(.[[1]] == starting_activity)
  ends <- apply(subset, 1, end_of_first_activity)
  right <- sample_one(ends) # Draw one end time randomly
  
  result[1:right] <- starting_activity
  right_activity <- starting_activity
  
  
  while (right < len){
    transitions <- unlist(apply(sequences, 1, activity_transition_order2, 
                                time = right, activity = right_activity, from = TRUE, margin = margin))
    
    activities <- names(transitions)
    previous <- sub(paste0("-",right_activity,"-.*"), "", activities)   
    after <- sub(paste0(".*-",right_activity,"-"), "", activities)  
    
    if (length(after) == 0){
      result[(right + 1) : min(right + margin, len)] <- right_activity
      right <- min(right + margin, len)
    } else {
      rle_res <- rle(result[1:right])$values
      
      if (length(rle_res) == 1){ # This means that there is no previous 2 activity, the current one is starting activity
        prev_activity <- "None"   # Pair with NONE in the previous activity
      } else {
        prev_activity <- rle_res[length(rle_res) - 1]   # Pair with previous activity
      }
      indices <- which(previous == prev_activity)
      
      if(length(indices) == 0){ # Cannot borrow information from previous activity
        # Go back to using this activity to predict next (like order 1).
        activity <- after[sample(1:length(transitions), size = 1)]
        potential_durations <- transitions[which(after == activity)]
        dur <- sample_logspline(potential_durations)
      } else { # We have order 2 matches
        activity <- after[sample_one(indices)]
        potential_durations <- transitions[which((previous == prev_activity) & (after == activity))]
        dur <- sample_logspline(potential_durations)
      }
      
      if (dur + right > len)  dur <- len - right
      
      result[(right + 1):(right + dur)] <- rep(activity, dur)
      right <- right + dur
      right_activity <- activity
    }
  }
  
  return(result)
}


# simulate_multiple_sequences <- function(data, cluster, n){
#   simulated <- replicate(n, simulate_one_sequence_order2(data, cluster))
#   result <- data.frame(t(simulated))
#   return(result)
# }


parallel_simulate_multiple_sequences <- function(data, cluster, n, seeds, margin = 60){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'MASS', 'clValid', 'logspline')) %dopar% {
    source("clust_logspline_order2.R")
    set.seed(seeds[i])
    simulated <- simulate_one_sequence_order2(data, cluster, margin) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}



