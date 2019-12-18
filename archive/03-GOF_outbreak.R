library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(abind)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

## Outbreak models ##
# Script inputs #
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced2")
maxYSOForPD <- 9 # Set to 12 for LP and 9 for SF
#_______________#

# Detection data #
Y <- str_c("Y.", stratum) %>% as.name %>% eval
TPeriod <- str_c("T.", stratum) %>% as.name %>% eval
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
gridID <- Cov[, "gridIndex"]
n.grid <- max(gridID)
n.point <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Covariates #
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

YSO.b <- Cov[, "YSO"]
YSO.b[which(!outbreak_grids[gridID])] <- -1
YSO.check <- YSO.b %>% # This is for screening out late-outbreak PctDead values
  replace(which(!outbreak_grids[gridID]), -1) %>%
  replace(., which(is.na(.)), tapply(., gridID, median, na.rm = T)[gridID][which(is.na(.))] %>% round)
YSO.b <- YSO.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.d[is.na(YSO.d)] <- 0
YSO.b[is.na(YSO.b)] <- YSO.d[gridID][is.na(YSO.b)]

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(YSO.check > maxYSOForPD)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & !outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & !outbreak_grids[gridID])])
PctDead.sd <- sd(PctDead.b[which(!outbreak_grids[gridID])]) %>% rep(length(PctDead.b)) # SDs for imputing missing values  for non-outbreak grids
PctDead.b[which(is.na(PctDead.b) & outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & outbreak_grids[gridID])])
PctDead.sd[which(outbreak_grids[gridID])] <- sd(PctDead.b[which(outbreak_grids[gridID])]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values

RCovAS.b <- Cov[, "RCOV_AS"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovAS.d <- tapply(RCovAS.b, gridID, mean, na.rm = T) # Grid-level values
RCovAS.b.missing <- is.na(RCovAS.b) %>% as.integer # Index missing values to be imputed
RCovAS.sd <- tapply(RCovAS.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovAS.sd==0))
  RCovAS.sd[which(RCovAS.sd==0)] <- min(RCovAS.sd[which(RCovAS.sd>0)])
RCovAS.lower <- min(RCovAS.b, na.rm = T) # Lower bound for imputing missing values
RCovAS.b[which(is.na(RCovAS.b))] <- RCovAS.d[gridID][which(is.na(RCovAS.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovES.b <- Cov[, "RCOV_ES"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovES.d <- tapply(RCovES.b, gridID, mean, na.rm = T) # Grid-level values
RCovES.b.missing <- is.na(RCovES.b) %>% as.integer # Index missing values to be imputed
RCovES.sd <- tapply(RCovES.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovES.sd==0))
  RCovES.sd[which(RCovES.sd==0)] <- min(RCovES.sd[which(RCovES.sd>0)])
RCovES.lower <- min(RCovES.b, na.rm = T) # Lower bound for imputing missing values
RCovES.b[which(is.na(RCovES.b))] <- RCovES.d[gridID][which(is.na(RCovES.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovPine.b <- Cov[, "RCOV_Pine"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovPine.d <- tapply(RCovPine.b, gridID, mean, na.rm = T) # Grid-level values
RCovPine.b.missing <- is.na(RCovPine.b) %>% as.integer # Index missing values to be imputed
RCovPine.sd <- tapply(RCovPine.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovPine.sd==0))
  RCovPine.sd[which(RCovPine.sd==0)] <- min(RCovPine.sd[which(RCovPine.sd>0)])
RCovPine.lower <- min(RCovPine.b, na.rm = T) # Lower bound for imputing missing values
RCovPine.b[which(is.na(RCovPine.b))] <- RCovPine.d[gridID][which(is.na(RCovPine.b))] # Insert means for imputing PctDead for points in non-outbreak grids

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
if(stratum == "LP")
  Time.b[point.list.LP == "LP-057-06"] <-
  mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
# Impute missing time value in LP data based on survey times for point likely completed just before and just after.
# Inferred which points after viewing the other points and their associated times in ArcMap.

# Get parameters #
d0 <- mod$sims.list$d0

b0 <- mod$sims.list$b0
bb.RCovAS <- mod$sims.list$bb.RCovAS
if(stratum == "SF") bb.RCovES <- mod$sims.list$bb.RCovES
bb.RCovPine <- mod$sims.list$bb.RCovPine
bb.pdead <- mod$sims.list$bb.pdead
bb.YSO <- mod$sims.list$bb.YSO
bb.YSO2 <- mod$sims.list$bb.YSO2
bb.pdXYSO <- mod$sims.list$bb.pdXYSO

a0 <- mod$sims.list$a0
ba.Time <- mod$sims.list$ba.Time
ba.Time2 <- mod$sims.list$ba.Time2
ba.DOY <- mod$sims.list$ba.DOY
ba.DOY2 <- mod$sims.list$ba.DOY2
ba.pdead <- mod$sims.list$ba.pdead
ba.YSO <- mod$sims.list$ba.YSO

# Calculate GOF p for deviance #
test.dev <- test.Pears <- c()
for(i in 1:mod$mcmc.info$n.samples) {
  prob.y <- Y*0
  for(sp in 1:length(spp.list)) {
    psi <- expit(d0[i, sp])
    theta <- expit(b0[i, sp] + bb.RCovAS[i, sp]*RCovAS.b + bb.RCovPine[i, sp]*RCovPine.b +
                     bb.pdead[i, sp]*PctDead.b + bb.YSO[i, sp]*YSO.b +
                     bb.YSO2[i, sp]*(YSO.b^2) + bb.pdXYSO[i, sp]*PctDead.b*YSO.b)
    if(stratum == "SF") theta <- expit(logit(theta) + bb.RCovES[i, sp]*RCovES.b)
    p <- expit(a0[i, sp] + ba.Time[i, sp]*Time.b + ba.Time[i, sp]*(Time.b^2) +
                 ba.DOY[i, sp]*DOY.b + ba.DOY[i, sp]*(DOY.b^2) +
                 ba.pdead[i, sp]*PctDead.b + ba.YSO[i, sp]*YSO.b)
    prob.y[, sp] <- psi * theta * p
  }
  DD.obs <- -2*sum(log(dbinom(Y, TPeriod, prob.y)))
  DP.obs <- sum(((Y - prob.y)^2) / prob.y)
  yrep.arr <- abind(matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp), # Need to repeat this for number of minute intervals (awkward!)
                    matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp),
                    matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp),
                    matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp),
                    matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp),
                    matrix(rbinom(n.point * n.spp, 1, prob.y), n.point, n.spp), along = 3)
  yrep <- apply(yrep.arr, c(1, 2), max)
  TPrep <- matrix(6, nrow = n.point, ncol = n.spp)
  TPrep[which(yrep == 1)] <- apply(matrix(yrep.arr, n.point * n.spp, 6)[which(yrep == 1),], 1,
                                   function(x) which(x == 1)[1])
  DD.rep <- -2*sum(log(dbinom(yrep, TPrep, prob.y)))
  DP.rep <- sum(((yrep - prob.y)^2) / prob.y)
  test.dev <- c(test.dev, (DD.obs > DD.rep)*1)
  test.Pears <- c(test.Pears, (DP.obs > DP.rep)*1)
}

p <- sum(test.dev) / length(test.dev)
p
p <- sum(test.Pears) / length(test.Pears)
p
saveObject(p, str_c("GOFp_", stratum, "_outbreak"))