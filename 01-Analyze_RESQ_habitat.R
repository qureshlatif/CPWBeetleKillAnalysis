library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

#____________________________________________________ Script inputs _____________________________________________________________#
stratum <- "LP" # Set stratum (LP or SF)
model.file <- "CPWBeetleKillAnalysis/model_RESQ_habitat_HZdist_LP.jags"

# MCMC values.  Adjust as needed.
nc <- 3
nb <- 5000
ni <- 130000
nt <- 30

save.out <- "mod_RESQ_habitat_HZdist_LP"
#________________________________________________________________________________________________________________________________#

data <- list("Y",
             "dclass", "mean.cl", "sd.cl", # needed for distance sampling
             "gridID", "nGrid", "nPoint",
             #"nInt", # needed for time removal
             "nInd", "nG", "area.band", "area.prop", "breaks", # needed for distance sampling
             "Time.b", "DOY.b",
             "ccov.b", "shcov.b", "RSC_Con.b", "GHerb.b", "Gwoody.b",
             "GDD.b", "RCovAS.b", "RCovES.b", "RCovPine.b") # Inputs for JAGS model.

parameters <- c("beta0.mean", "beta0.sd", "N.mean", "p.mean", # Assemble the parameters vector for JAGS (What we want to track).
                "beta0", "N", "cl.size",
                "bl.ccov", "bl.shcov", "bl.RSC_Con", "bl.GHerb", "bl.Gwoody", "bl.GDD",
                "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
                ##___ Hazard rate parameters ___##
                "a0", "a.Time", "a.Time2",
                "a.DOY", "a.DOY2",
                "a.pdead", "a.YSO", "a.YSO2", "a.pdXYSO",
                "a.ccov", "a.shcov",
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

RCovAS.b <- Cov[, "RCOV_AS"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovAS.d <- tapply(RCovAS.b, gridID, mean, na.rm = T) # Grid-level values
RCovAS.b[which(is.na(RCovAS.b))] <- RCovAS.d[gridID][which(is.na(RCovAS.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovES.b <- Cov[, "RCOV_ES"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovES.d <- tapply(RCovES.b, gridID, mean, na.rm = T) # Grid-level values
RCovES.b[which(is.na(RCovES.b))] <- RCovES.d[gridID][which(is.na(RCovES.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovPine.b <- Cov[, "RCOV_Pine"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovPine.d <- tapply(RCovPine.b, gridID, mean, na.rm = T) # Grid-level values
RCovPine.b[which(is.na(RCovPine.b))] <- RCovPine.d[gridID][which(is.na(RCovPine.b))] # Insert means for imputing PctDead for points in non-outbreak grids

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.d <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.b[is.na(ccov.b)] <- ccov.d[gridID[is.na(ccov.b)]]

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.d <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.b[is.na(shcov.b)] <- shcov.d[gridID[is.na(shcov.b)]]

RSC_Con.b <- 100 - Cov[, "RCShrb_UD"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RSC_Con.d <- tapply(RSC_Con.b, gridID, mean, na.rm = T) # Grid-level values
RSC_Con.b[which(is.na(RSC_Con.b))] <- RSC_Con.d[gridID][which(is.na(RSC_Con.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GHerb.b <- Cov[, "HerbCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GHerb.d <- tapply(GHerb.b, gridID, mean, na.rm = T) # Grid-level values
GHerb.b[which(is.na(GHerb.b))] <- GHerb.d[gridID][which(is.na(GHerb.b))] # Insert means for imputing PctDead for points in non-outbreak grids

Gwoody.b <- Cov[, "WoodyCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
Gwoody.d <- tapply(Gwoody.b, gridID, mean, na.rm = T) # Grid-level values
Gwoody.b[which(is.na(Gwoody.b))] <- Gwoody.d[gridID][which(is.na(Gwoody.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GDD.b <- Cov[, "DDCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GDD.d <- tapply(GDD.b, gridID, mean, na.rm = T) # Grid-level values
GDD.b[which(is.na(GDD.b))] <- GDD.d[gridID][which(is.na(GDD.b))] # Insert means for imputing PctDead for points in non-outbreak grids

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
