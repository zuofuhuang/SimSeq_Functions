library(mhsmm)
library(TraMineR)

sequences <- read.csv("../sequences.csv")[,-1]
N <- nrow(sequences)
sequences[is.na(sequences)] <- "N/A"
seq_distances <- read.csv("../1929sequences_dist_om.csv")[,-1]
colnames(sequences) <- c(format(seq(as.POSIXct("00:00:00", format = "%T"), 
                                    as.POSIXct("23:59:00", format = "%T"), by = "1 min"), "%H:%M"))
seq_tminer <- seqdef(sequences)
states <- alphabet(seq_tminer)
transition_mat <- seqtrate(seq_tminer, with.missing = FALSE)
initial <- (table(seq_tminer[,1])[1:13])/N

all_activities <- read.csv("../all_activities.csv")
activity_dist <- cbind(states, NA, NA)
for (j in 1:length(states)){
  this <- all_activities[all_activities$X1 == (states[j]),]
  activity_dist[j,2] <- mean(this$lengths)
  activity_dist[j,3] <- sd(this$lengths)
}


# For d, nrow is the number of states; ncol is 1929?
d <- 

startval <- hmmspec(init = initial, 
                    trans = transition_mat, 
                    parms.emis = b0, 
                    sojourn = list(d = d, type='nonparametric'),
                    dens.emis = dnorm.hsmm)