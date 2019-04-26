library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

#____________________________________________________ Script inputs _____________________________________________________________#
stratum <- "LP" # Set stratum (LP or SF)
maxYSOForPD <- 12 # Set to 12 for LP and 9 for SF
model.file <- "CPWBeetleKillAnalysis/model_RESQ_outbreak_HZdist_LP.jags"

# MCMC values.  Adjust as needed.
nc <- 3
nb <- 5000
ni <- 40000
nt <- 10

save.out <- "mod_RESQ_outbreak_HZdist_LP"
#________________________________________________________________________________________________________________________________#

data <- list("Y",
             "dclass", "mean.cl", "sd.cl", # needed for distance sampling
             "gridID", "nGrid", "nPoint",
             #"nInt", # needed for time removal
             "nInd", "nG", "area.band", "area.prop", "breaks", # needed for distance sampling
             "Time.b", "DOY.b",
             "PctDead.b", "YSO.b", "YSO.mins", "YSO.maxs", "YSO.missing",
             "PctDead.sd", "PctDead.lower", "PctDead.b.missing",
             "heatload.d", "TWI.d", "TWIP.d", "RDens.d", "WILD.d",
             "RCovAS.d", "RCovAS.b", "RCovAS.sd", "RCovAS.b.missing", "RCovAS.lower",
             "RCovES.d", "RCovES.b", "RCovES.sd", "RCovES.b.missing", "RCovES.lower",
             "RCovPine.d", "RCovPine.b", "RCovPine.sd", "RCovPine.b.missing", "RCovPine.lower") # Inputs for JAGS model.

parameters <- c("beta0.mean", "beta0.sd", "N.mean", "p.mean", # Assemble the parameters vector for JAGS (What we want to track).
                "beta0", "N", "cl.size",
                "bd.heatload", "bd.TWI", "bd.TWIP", "bd.RDens", "bd.WILD",
                "bl.pdead", "bl.YSO", "bl.YSO2", "bl.pdXYSO",
                "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
                ##___ Hazard rate parameters ___##
                "a0", "a.Time", "a.Time2",
                "a.DOY", "a.DOY2",
                "a.pdead", "a.YSO", "a.YSO2", "a.pdXYSO",
                "b",
                ##______________________________##
                ##___ Half-normal parameters or time removal ___##
                "bt.Int", # Time removal only - linear effect of minutes elapsed
                "bt.0", "bt.Time", "bt.Time2",
                "bt.DOY", "bt.DOY2",
                "bt.ccov", "bt.shcov",
                ##______________________________##
                "lambda", "a", # Needed for WAIC
                "test") # GOF

## For distance sampling ##
Y <- eval(as.name(str_c("Y.", stratum, ".dist"))) # For distance sampling
dclass <- eval(as.name(str_c("dclass.", stratum)))
mean.cl <- max(1.001, mean(dclass[, "CL_Count"]))
sd.cl <- max(0.001, sd(dclass[, "CL_Count"]))
dimnames(dclass) <- NULL
nInd <- nrow(dclass)
nPoint <- length(Y)
inits <- function() # Setting these based on posterior distribution from an initial successful run.
  list(N = Y, beta0.mean = rnorm(1, 1, 0.1), beta0.sd = rnorm(1, 0.66, 0.07),
       a0 = rnorm(1, 3.5, 0.5), b = rnorm(1, 3.16, 0.5)) # for hazard rate model
       #bt.0 = rnorm(1, 3.4, 0.03)) # for half-normal model

## For time removal ##
#Y <- eval(as.name(str_c("Y.", stratum, ".trem.all"))) # For time removal
#nInt <- dim(Y)[2]
#nPoint <- dim(Y)[1]
#inits <- function() # Setting these based on posterior distribution from an initial successful run.
#  list(N = apply(Y, 1, sum), n = apply(Y, 1, sum), beta0.mean = rnorm(1, 1, 0.1), beta0.sd = rnorm(1, 0.66, 0.07),
#       bt.0 = rnorm(1, 0, 1)) # for time removal model
#____________________________________________________________________________________________________________________________________#

# Detection data #
gridID <- eval(as.name(str_c("Cov.", stratum)))[, "gridIndex"]
nGrid <- max(gridID)

# Covariates #
Cov <- eval(as.name(str_c("Cov.", stratum)))

outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

YSO.b <- Cov[, "YSO"]
YSO.b[which(!outbreak_grids[gridID])] <- -1
YSO.check <- YSO.b %>% # This is for screening out late-outbreak PctDead values
  replace(which(!outbreak_grids[gridID]), -1) %>%
  replace(., which(is.na(.)), tapply(., gridID, median, na.rm = T)[gridID][which(is.na(.))] %>% round)
YSO.b <- YSO.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.mins <- tapply(YSO.b, gridID, min, na.rm = T) %>% as.numeric # Grid-level values
YSO.maxs <- tapply(YSO.b, gridID, max, na.rm = T) %>% as.numeric  # Grid-level values
ind <- which(YSO.mins == YSO.maxs)
YSO.mins[ind] <- YSO.mins[ind] - 0.01
YSO.maxs[ind] <- YSO.maxs[ind] + 0.01
rm(ind)
YSO.missing <- (is.na(YSO.b) & outbreak_grids[gridID]) %>% as.integer
YSO.b[is.na(YSO.b)] <- 0
YSO.d[is.na(YSO.d)] <- 0

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

heatload.d <- Cov[, "heatload"]  %>% # Grid-level values only
  tapply(gridID, mean, na.rm = T) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- Cov[, "TWI"]  %>% # Grid-level values only
  tapply(gridID, mean, na.rm = T) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWIP.d <- Cov[, "TWIP"]  %>% # Grid-level values only
  tapply(gridID, mean, na.rm = T) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

RDens.d <- Cov[, "Rd_dens1km"]  %>% # Grid-level values only
  tapply(gridID, mean, na.rm = T) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

WILD.d <- Cov[, "WILD"]  %>% # Grid-level values only
  tapply(gridID, function(x) (mean(x, na.rm = T) >= 0.5)*1)


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
if(stratum == "LP") Time.b[point.list.LP == "LP-057-06"] <- 
  mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
  # Imputed missing time value in LP data based on survey times for point likely completed just before and just after.
  # Inferred which points after viewing the other points and their associated times in ArcMap.

# Save workspace for model checking #
#save.image("RESQ_GOF_workspace.RData")

st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

max(out$summary[which(!is.na(out$summary[ ,"Rhat"])) ,"Rhat"])

min(out$summary[,"n.eff"])
#sort(out$summary[,"n.eff"])[1:200]

sum(out$sims.list$test) / out$mcmc.info$n.samples

library(R.utils)
saveObject(out, save.out)
