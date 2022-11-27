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
  end <- cumSum[as.numeric(cutted)+1]
  return(c(start, end, end - start + 1))
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

  if (!setequal(sequence[position],activity)){
    stop("Activity state does not match")
  }
  
  # Find out the length of the subsequence that contains position
  lens <- rle(sequence)$length
  seq_cumSum <- c(0,cumsum(lens))
  cutted <- cut(position, seq_cumSum)
  levels(cutted) <- 1:length(lens)
  position <- lens[cutted]
  names(position) <- activity
  
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
simulate_one_sequence <- function(sequences, margin = 60){
  len <- ncol(sequences)
  result <- rep(NA, len)
  
  # We could pass states/freq_table as a separate input variable. Otherwise, this would be run each time we generate a sequence
  states <- seqstatl(sequences)
  freq_table <- create_freq_table(sequences) 
  
  
  # Step 1
  starting_time <- sample(1:ncol(sequences), 1)
  probs <- freq_table[,starting_time]
  starting_activity <- sample(states, 1, prob = probs)
  
  # Step 2
  subset <- filter(sequences, sequences[,starting_time] == starting_activity)
  
  # Method #1: Generate activity from the midpoint and length (which follows normal distribution)
  subsequences <- data.frame(t(apply(subset, 1, find_subsequence, activity = starting_activity, time = starting_time)))
  subsequences$midpoint <- (subsequences$X1 + subsequences$X2)/2 # Don't need to round just yet.
  
  # X3 is the length
  cov_mat <- cov(subsequences[,c(3,4)])
  binorm <- rtmvnorm(n = 1, mean = c(mean(subsequences$X3), mean(subsequences$midpoint)),
                     sigma = (cov_mat + diag(2)*0.01), lower = c(1,1), upper = c(len,len)) 
  left <- round(max(binorm[2] - binorm[1]/2, 1))
  right <- round(min(binorm[2] + binorm[1]/2, len))
  result[left:right] <- starting_activity
  
  left_activity <- starting_activity
  right_activity <- starting_activity
  
  while (left > 1){
    transitions <- unlist(apply(sequences, 1, activity_transition, 
                                time = left, activity = left_activity, from = FALSE, margin = margin))
    
    # It's possible that none of the sequences contains a transition at this time interval at all.
    # If there is no transition, then we are staying with the same activity.
    if (is.null(transitions)){
      result[max(left - margin, 1):(left - 1)] <- left_activity
      left <- max(left - margin, 1)
    } else {
      # How many times other activities transition into this activity
      # activity_table <- table(names(transitions))
      # Currently, we will randomly pick one from the sequences
      rand <- sample(1:length(transitions), size = 1)
      dur <- transitions[rand]
      activity <- names(transitions)[rand]
      
      if (left - dur < 1)  dur <- left - 1
      
      result[(left - dur):(left - 1)] <- rep(activity, dur)
      left <- left - dur
      left_activity <- activity
    }
  }
  
  while (right < len){
    transitions <- unlist(apply(sequences, 1, activity_transition, 
                                time = right, activity = right_activity, from = TRUE, margin = margin))
    
    if (is.null(transitions)){
      result[(right + 1) : min(right + margin, len)] <- right_activity
      right <- min(right + margin, len)
    } else {
      rand <- sample(1:length(transitions), size = 1)
      dur <- transitions[rand]
      activity <- names(transitions)[rand]
      
      if (dur + right > len)  dur <- len - right
      
      result[(right + 1):(right + dur)] <- rep(activity, dur)
      right <- right + dur
      right_activity <- activity
    }
  }
  
  return(result)
}


simulate_multiple_sequences <- function(data, n){
  simulated <- replicate(n, simulate_one_sequence(data))
  result <- data.frame(t(simulated))
  return(result)
}


parallel_simulate_multiple_sequences <- function(data, n, seeds){
  cores <- detectCores()
  cl <- makeCluster(cores - 1)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:n, .combine = cbind, .packages = c('TraMineR','dplyr', 'tmvtnorm')) %dopar% {
    source("Functions.R")
    set.seed(seeds[i])
    simulated <- simulate_one_sequence(data) # calling a function
    simulated
  }
  
  result <- data.frame(t(result))
  stopCluster(cl)
  return(result)
}

