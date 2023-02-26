reorder_cluster_assignment <- function(assignment){
  tbl <- table(assignment)
  groups <- length(tbl)
  
  list <- list()
  for (i in 1:groups){
    list[[i]] <- which(assignment == i)
  }
  
  sorted_tbl <- sort(tbl, decreasing = TRUE)
  cluster_num_by_size <- as.numeric(names(sorted_tbl))
  
  for (j in 1:groups){
    assignment[list[[(cluster_num_by_size[j])]]] <- j
  }
  
  return(assignment)
}



# Distribution of number of transitions each day
transitions_per_day <- function(sequence){
  return(length((rle(sequence)$length)))
}


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
  final_cluster <- cutree(hclust_real, k = num[1])
  return(final_cluster)
}

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
