library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

## Lodgepole pine stratum ##
Y <- Y.LP.dist
nsites <- length(Y)

mods <- c("mod_RESQ_allcovs_HNdist_LP", "mod_RESQ_allNredpcovs_HNdist_LP",
          "mod_RESQ_allcovs_HZdist_LP", "mod_RESQ_allNredpcovs_HZdist_LP",
          "mod_RESQ_HZp(red_LP)_N(Pdead_pDeadXYSO_NdropXYSO)")
cols <- c("WAIC", "DIC")
out <- matrix(NA, nrow = length(mods), ncol = length(cols),
              dimnames = list(mods, cols))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  out[mods[m], "DIC"] <- mod$DIC
  lambda <- mod$sims.list$lambda
  p <- mod$sims.list$pcap
  lpd <- pwaic <- rep(NA, nsites)
  for(i in 1:nsites) {
    lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
    pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  }
  out[mods[m], "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))
}

out <- as.data.frame(out) %>%
  mutate(mods = mods) %>%
  mutate(dWAIC = WAIC - min(WAIC)) %>%
  mutate(dDIC = DIC - min(DIC)) %>%
  arrange(dWAIC) %>%
  select(mods, WAIC, dWAIC, DIC, dDIC)

write.csv(out, "Model_selection_RESQ_LP.csv", row.names = T)

## Spruce fir stratum ##
Y <- Y.SF.dist
nsites <- length(Y)

mods <- c("mod_RESQ_allcovs_HNdist_SF", "mod_RESQ_allNredpcovs_HNdist_SF",
          "mod_RESQ_allcovs_HZdist_SF", "mod_RESQ_allNredpcovs_HZdist_SF",
          "mod_RESQ_HZp(red_SF)_N(Pdead_pDeadXYSO_NdropXYSO)")
cols <- c("WAIC", "DIC")
out <- matrix(NA, nrow = length(mods), ncol = length(cols),
              dimnames = list(mods, cols))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  out[mods[m], "DIC"] <- mod$DIC
  lambda <- mod$sims.list$lambda
  p <- mod$sims.list$pcap
  SFd <- pwaic <- rep(NA, nsites)
  for(i in 1:nsites) {
    SFd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
    pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  }
  out[mods[m], "WAIC"] <- -2 * (sum(SFd) - sum(pwaic))
}

out <- as.data.frame(out) %>%
  mutate(mods = mods) %>%
  mutate(dWAIC = WAIC - min(WAIC)) %>%
  mutate(dDIC = DIC - min(DIC)) %>%
  arrange(dWAIC) %>%
  select(mods, WAIC, dWAIC, DIC, dDIC)

write.csv(out, "Model_selection_RESQ_SF.csv", row.names = T)
