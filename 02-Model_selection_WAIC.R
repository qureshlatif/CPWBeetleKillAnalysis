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

mod_HN <- loadObject("mod_RESQ_allcovs_HNdist_LP")
mod_HZ <- loadObject("mod_RESQ_allcovs_HZdist_LP")

rows <- c("Hazard rate", "Half normal")
cols <- c("WAIC", "DIC")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

mod <- mod_HN
out["Half normal", "DIC"] <- mod$DIC
lambda <- mod$sims.list$lambda
p <- mod$sims.list$pcap
lpd <- pwaic <- rep(NA, nsites)
for(i in 1:nsites) {
  lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
}
out["Half normal", "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))

mod <- mod_HZ
out["Hazard rate", "DIC"] <- mod$DIC
lambda <- mod$sims.list$lambda
p <- mod$sims.list$pcap
lpd <- pwaic <- rep(NA, nsites)
for(i in 1:nsites) {
  lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
}
out["Hazard rate", "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))

out <- as.data.frame(out)
out$dWAIC <- out$WAIC - min(out$WAIC)
out$dDIC <- out$DIC - min(out$DIC)
out <- out %>% select(WAIC, dWAIC, DIC, dDIC)

write.csv(out, "Model_selection_RESQ_LP.csv", row.names = T)

## Spruce fir stratum ##
Y <- Y.SF.dist
nsites <- length(Y)

mod_HN <- loadObject("mod_RESQ_allcovs_HNdist_SF")
mod_HZ <- loadObject("mod_RESQ_allcovs_HZdist_SF")

rows <- c("Hazard rate", "Half normal")
cols <- c("WAIC", "DIC")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

mod <- mod_HN
out["Half normal", "DIC"] <- mod$DIC
lambda <- mod$sims.list$lambda
p <- mod$sims.list$pcap
lpd <- pwaic <- rep(NA, nsites)
for(i in 1:nsites) {
  lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
}
out["Half normal", "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))

mod <- mod_HZ
out["Hazard rate", "DIC"] <- mod$DIC
lambda <- mod$sims.list$lambda
p <- mod$sims.list$pcap
lpd <- pwaic <- rep(NA, nsites)
for(i in 1:nsites) {
  lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
  pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * p[, i])))
}
out["Hazard rate", "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))

out <- as.data.frame(out)
out$dWAIC <- out$WAIC - min(out$WAIC)
out$dDIC <- out$DIC - min(out$DIC)
out <- out %>% select(WAIC, dWAIC, DIC, dDIC)

write.csv(out, "Model_selection_RESQ_SF.csv", row.names = T)
