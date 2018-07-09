library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("RESQ_GOF_workspace.RData")

mod <- loadObject("mod_RESQ_allcovs_TR_LP")
Y <- Y.LP.trem.all

## posterior predictive GOF wrt time removal ##
T.obs <- apply(Y, 2, sum, na.rm = T)

nsim <- mod$mcmc.info$n.samples
T.rep <- matrix(NA, nrow = nsim, ncol = length(T.obs))
for(i in 1:nsim) {
  # Predict new Y #
  Y.new <- array(NA, dim = dim(Y))
  pi <- mod$sims.list$p[i, ] %>% as.matrix()
  for(t in 2:nInt) pi <-
    cbind(pi, mod$sims.list$p[i, ] *
            apply((1 - mod$sims.list$p[i, ]) %>% as.matrix,
                  1, function(x) x^(t-1)))
  pcap <- apply(pi, 1, sum)
  n.new <- rbinom(nPoint, size = mod$sims.list$N[i, ], prob = pcap)
  pic <- array(NA, dim = dim(pi))
  pic <- pi / pcap
  for(n in 1:nrow(Y.new)) Y.new[n, ] <- rmultinom(1, n.new[n], pic[n, ]) # Might be able to do without the loop
  #Y.new <- rmultinom(nPoint, n.new, pic) # Try this next time??
  
  T.rep[i, ] <- apply(Y.new, 2, sum, na.rm = T)
}

Trep.md <- apply(T.rep, 2, median)
Trep.95lo <- apply(T.rep, 2, function(x) quantile(x, prob = 0.025, type = 8))
Trep.95hi <- apply(T.rep, 2, function(x) quantile(x, prob = 0.975, type = 8))
dat.plt <- data.frame(Time = 1:6, T.obs, Trep.md, Trep.95lo, Trep.95hi)

ggplot(dat.plt, aes(x = Time, y = Trep.md)) # Not finished