## Distribution plots

### Example:
# sequences[sequences == "N/A"] <- NA
# p_sequences <- seqplot(seqdef(sequences), border = NA, type = "d", cex.legend = 0.5)

# this <- as.data.frame(sim_list[[1]])
# this[this == "N/A"] <- NA
# seqplot(seqdef(this), border = NA, type = "d", cex.legend = 0.5)


## Percentage of each activity

### Example:
# a_sequences <- seqmeant(sequences_T)


## Distribution of duration of each activity

### Example:
dist_of_activity_dur <- function(traminer_seq){
  dist <- seqistatd(traminer_seq)
  means <- apply(dist, 2, mean)
  sds <- apply(dist, 2, sd)
  return(list(means,sds))
}

# dist1 <- seqistatd(sequences_T)
# c1 <- apply(dist1, 2, mean)
# d1 <- apply(dist1, 2, sd)


## The distribution of number of distinct states in each sequence

num_of_distinct_states <- function(traminer_seq){
  return(seqlength(seqdss(traminer_seq)))
}



## How much work/home does each day have?

length_of_certain_activity <- function(traminer_seq, activity){
  dist <- seqistatd(traminer_seq)
  return(dist[,activity])
}

### Example:
# work1 <- dist1[,"WORK"]


# dat <- as.data.frame(cbind(c(work1, work6, work7, work8, work9, work10, work11), c(rep("Original data",1929),rep("TVMC",1929*2), rep("Our Method",1929*2), rep("Our Method + order 2", 1929*2)), c(rep("Clustered", 1929), rep("Clustered",1929),rep("Not Clustered",1929),rep("Clustered",1929),rep("Not Clustered",1929),rep("Clustered",1929),rep("Not Clustered",1929))))
# 
# colnames(dat) <- c("Length", "Method","Clustered")
# dat$Length <- as.numeric(dat$Length)

# ggplot(data = dat, aes(x = Length)) +
#   facet_grid(Clustered ~ Method) +
#   geom_histogram(binwidth = 50) +
#   theme_bw()

# we can think of adding the individual level through clustering as well.


## Number of distinct activity periods per day
num_distinct_periods_perday <- function(traminer_seq, activity){
  return(table(rowSums(seqdss(traminer_seq) == activity)))
}

### Example:
# work_num1 <- table(rowSums(seqdss(sequences_T) == "WORK"))



## Sum of percentage of each state squared
# ref and data are traminer objects
sum_of_state_pct <- function(ref, data){
  ref_prob <- seqstatd(ref)[["Frequencies"]]
  data_prob <- seqstatd(data)[["Frequencies"]]
  
  dist <- 0
  for (i in 1:ncol(ref_prob)){
    this_ref <- as.vector(ref_prob[,i])
    this_data <- as.vector(data_prob[,i])
    this_dist <- sum((this_ref - this_data)^2)
    dist <- dist + this_dist
  }
  return(dist)
}


### Example: 
# sum_of_state_pct(sequences_T, cluster_T)



create_ecdf_data <- function(distances, len){
  vecs <- c()
  position <- c()
  category <- c()
  ref <- ecdf(distances[[1]])(0:len)
  for (i in 2:length(distances)){
    vecs <- c(vecs, ecdf(distances[[i]])(0:len) - ref)
    position <- c(position, 0:len)
    category <- c(category, rep(i, len+1))
  }
  dat <- as.data.frame(cbind(vecs, position, category))
  dat$category <- as.factor(dat$category)
  return(dat)
}




## Entropy

### Example:
# entropy1 <- TraMineR::seqient(sequences_T)


## Transition Matrix?

# tm1 <- seqtrate(sequences_T)
# tm6 <- seqtrate(sim_seqs_cluster_TVMC_T)
# tm7 <- seqtrate(sim_seqs_noCluster_TVMC_T)
# tm8 <- seqtrate(sim_seqs_imputed_cluster_T)
# tm9 <- seqtrate(sim_seqs_imputed_noCluster_T)
# tm10 <- seqtrate(sim_seqs_imputed_cluster_order2_T)
# tm11 <- seqtrate(sim_seqs_imputed_noCluster_order2_T)

# heatmap(tm1 - tm6, Rowv = NA, Colv = NA)
# heatmap(tm1 - tm7, Rowv = NA, Colv = NA)
# heatmap(tm1 - tm8, Rowv = NA, Colv = NA)
# heatmap(tm1 - tm9, Rowv = NA, Colv = NA)
# heatmap(tm1 - tm10, Rowv = NA, Colv = NA)
# heatmap(tm1 - tm11, Rowv = NA, Colv = NA)


## Durations of each state (to show that they are not exponentially distributed)

# find_lengths <- function(seq, state){
#   rled <- rle(seq)
#   lengths <- rled$lengths[which(rled$values == state)]
#   return(lengths)
# }
# 
# works <- unlist(apply(sequences, 1, find_lengths, state = "WORK"))
# homes <- unlist(apply(sequences, 1, find_lengths, state = "HOME"))
# cars <- unlist(apply(sequences, 1, find_lengths, state = "CAR"))
# 
# g1 <- ggplot(data = data.frame(works), aes(x = works)) + 
#   geom_density() +
#   scale_x_continuous(limits = c(0, 1440)) +
#   theme_bw() + 
#   xlab("Duration of work episodes") + 
#   ggtitle("Distribution of work episodes")
# 
# g2 <- ggplot(data = data.frame(homes), aes(x = homes)) + 
#   geom_density() +
#   scale_x_continuous(limits = c(0, 1440)) +
#   theme_bw() + 
#   xlab("Duration of home episodes") + 
#   ggtitle("Distribution of home episodes")
# 
# g3 <- ggplot(data = data.frame(cars), aes(x = cars)) + 
#   geom_density() + 
#   scale_x_continuous(limits = c(0, 1440)) +
#   theme_bw() + 
#   xlab("Duration of car episodes") + 
#   ggtitle("Distribution of car episodes")
#   
# 
# grid.arrange(g1,g2,g3, nrow = 3)





# sim_list <- list(sequences, cluster_TVMC, noCluster_TVMC, cluster, noCluster, cluster_order2, noCluster_order2, imputed_cluster, imputed_noCluster, imputed_cluster_order2, imputed_noCluster_order2, imputed_KDE_cluster, imputed_KDE_noCluster, imputed_KDE_cluster_order2, imputed_KDE_noCluster_order2, imputed_log_cluster, imputed_log_noCluster, imputed_log_cluster_order2, imputed_log_noCluster_order2)

# sim_list_T <- list(sequences_T, cluster_TVMC_T, noCluster_TVMC_T, cluster_T, noCluster_T, cluster_order2_T, noCluster_order2_T, imputed_cluster_T, imputed_noCluster_T, imputed_cluster_order2_T, imputed_noCluster_order2_T, imputed_KDE_cluster_T, imputed_KDE_noCluster_T, imputed_KDE_cluster_order2_T, imputed_KDE_noCluster_order2_T, imputed_log_cluster_T, imputed_log_noCluster_T, imputed_log_cluster_order2_T, imputed_log_noCluster_order2_T)


