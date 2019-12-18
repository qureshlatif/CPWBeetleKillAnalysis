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
mod <- loadObject("mod_SFcommunity_habitat_reduced")
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

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > maxYSOForPD)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & !outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & !outbreak_grids[gridID])])
PctDead.sd <- sd(PctDead.b[which(!outbreak_grids[gridID])]) %>% rep(length(PctDead.b)) # SDs for imputing missing values  for non-outbreak grids
PctDead.b[which(is.na(PctDead.b) & outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & outbreak_grids[gridID])])
PctDead.sd[which(outbreak_grids[gridID])] <- sd(PctDead.b[which(outbreak_grids[gridID])]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed

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

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.d <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- ccov.d[gridID[is.na(ccov.b)]]

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.d <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- shcov.d[gridID[is.na(shcov.b)]]

RSC_Con.b <- Cov[, "RCShrb_UC"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RSC_Con.d <- tapply(RSC_Con.b, gridID, mean, na.rm = T) # Grid-level values
RSC_Con.b.missing <- is.na(RSC_Con.b) %>% as.integer # Index missing values to be imputed
RSC_Con.sd <- tapply(RSC_Con.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RSC_Con.sd==0))
  RSC_Con.sd[which(RSC_Con.sd==0)] <- min(RSC_Con.sd[which(RSC_Con.sd>0)])
RSC_Con.lower <- min(RSC_Con.b, na.rm = T) # Lower bound for imputing missing values
RSC_Con.b[which(is.na(RSC_Con.b))] <- RSC_Con.d[gridID][which(is.na(RSC_Con.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GHerb.b <- Cov[, "HerbCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GHerb.d <- tapply(GHerb.b, gridID, mean, na.rm = T) # Grid-level values
GHerb.b.missing <- is.na(GHerb.b) %>% as.integer # Index missing values to be imputed
GHerb.sd <- tapply(GHerb.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GHerb.sd==0))
  GHerb.sd[which(GHerb.sd==0)] <- min(GHerb.sd[which(GHerb.sd>0)])
GHerb.lower <- min(GHerb.b, na.rm = T) # Lower bound for imputing missing values
GHerb.b[which(is.na(GHerb.b))] <- GHerb.d[gridID][which(is.na(GHerb.b))] # Insert means for imputing PctDead for points in non-outbreak grids

Gwoody.b <- Cov[, "WoodyCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
Gwoody.d <- tapply(Gwoody.b, gridID, mean, na.rm = T) # Grid-level values
Gwoody.b.missing <- is.na(Gwoody.b) %>% as.integer # Index missing values to be imputed
Gwoody.sd <- tapply(Gwoody.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(Gwoody.sd==0))
  Gwoody.sd[which(Gwoody.sd==0)] <- min(Gwoody.sd[which(Gwoody.sd>0)])
Gwoody.lower <- min(Gwoody.b, na.rm = T) # Lower bound for imputing missing values
Gwoody.b[which(is.na(Gwoody.b))] <- Gwoody.d[gridID][which(is.na(Gwoody.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GDD.b <- Cov[, "DDCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GDD.d <- tapply(GDD.b, gridID, mean, na.rm = T) # Grid-level values
GDD.b.missing <- is.na(GDD.b) %>% as.integer # Index missing values to be imputed
GDD.sd <- tapply(GDD.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GDD.sd==0))
  GDD.sd[which(GDD.sd==0)] <- min(GDD.sd[which(GDD.sd>0)])
GDD.lower <- min(GDD.b, na.rm = T) # Lower bound for imputing missing values
GDD.b[which(is.na(GDD.b))] <- GDD.d[gridID][which(is.na(GDD.b))] # Insert means for imputing PctDead for points in non-outbreak grids

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
bb.CanCov <- mod$sims.list$bb.CanCov
bb.ShCov <- mod$sims.list$bb.ShCov
bb.RSC_Con <- mod$sims.list$bb.RSC_Con
bb.GHerb <- mod$sims.list$bb.GHerb
bb.Gwoody <- mod$sims.list$bb.Gwoody
bb.GDD <- mod$sims.list$bb.GDD

a0 <- mod$sims.list$a0
ba.Time <- mod$sims.list$ba.Time
ba.Time2 <- mod$sims.list$ba.Time2
ba.DOY <- mod$sims.list$ba.DOY
ba.DOY2 <- mod$sims.list$ba.DOY2
ba.ccov <- mod$sims.list$ba.ccov
ba.shcov <- mod$sims.list$ba.shcov

# Calculate GOF p for deviance #
test.dev <- test.Pears <- c()
for(i in 1:mod$mcmc.info$n.samples) {
  prob.y <- Y*0
  for(sp in 1:length(spp.list)) {
    psi <- expit(d0[i, sp])
    theta <- expit(b0[i, sp] + bb.RCovAS[i, sp]*RCovAS.b + bb.RCovPine[i, sp]*RCovPine.b +
                     bb.CanCov[i, sp]*ccov.b + bb.ShCov[i, sp]*shcov.b +
                     bb.RSC_Con[i, sp]*RSC_Con.b + bb.GHerb[i, sp]*GHerb.b +
                     bb.Gwoody[i, sp]*Gwoody.b + bb.GDD[i, sp]*GDD.b)
    if(stratum == "SF") theta <- expit(logit(theta) + bb.RCovES[i, sp]*RCovES.b)
    p <- expit(a0[i, sp] + ba.Time[i, sp]*Time.b + ba.Time[i, sp]*(Time.b^2) +
                 ba.DOY[i, sp]*DOY.b + ba.DOY[i, sp]*(DOY.b^2) +
                 ba.ccov[i, sp]*ccov.b + ba.shcov[i, sp]*shcov.b)
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