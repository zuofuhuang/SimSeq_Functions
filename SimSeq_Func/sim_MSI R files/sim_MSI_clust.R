library(dplyr)
library(lubridate)
library(TraMineR)
library(performance)
library(GLMMadaptive)
library(pscl)
library(rbenchmark)
library(foreach)
library(doParallel)
library(fastcluster)
library(MASS)
library(clValid)
library(stringr)
library(clv)
library(TraMineRextras)
library(logspline)

source("clust.R")

cluster_assignment <- read.csv("cluster_assignment_for_MSI.csv")
cluster_assignment <- as.vector(cluster_assignment[,2])
sequences <- read.csv("original_sequences.csv")[,-1]

sim_seqs_noCluster <- parallel_simulate_multiple_sequences(sequences, rep(1,1929), n = 1929, seeds = 1:1929)
write.csv(sim_seqs_noCluster, "sim_seqs_noCluster.csv")

sim_seqs_cluster <- parallel_simulate_multiple_sequences(sequences, cluster_assignment, n = 1929, seeds = 1:1929 + 1929*1)
write.csv(sim_seqs_cluster, "sim_seqs_cluster.csv")


