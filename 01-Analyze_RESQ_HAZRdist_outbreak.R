library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

# Variables #
stratum <- "LP" # Set stratum

model.file <- "model_RESQ_allcovs_HZdist_outbreak.jags"

data <- list("Y", "dclass", "gridID", "nGrid", "nPoint", "nInd", "nG", "area.band", "area.prop", "breaks",
             "Time.b", "DOY.b", "ccov.b", "shcov.b", "PctDead.b", "YSO.b", "Outbrk.b", "TWIP.b",
             "PctDead.sd", "PctDead.lower", "TWIP.d", "TWIP.sd", "ccov.means", "ccov.sd", "shcov.means", "shcov.sd",
             "PctDead.b.missing", "TWIP.b.missing", "shcov.b.missing", "ccov.b.missing"
             ) # Inputs for JAGS model.

parameters <- c("beta0.mean", "beta0.sd", "N.mean", "p.mean", # Assemble the parameters vector for JAGS (What we want to track).
                "beta0", "N",
                "bl.pdead",
                "bl.pdead2", "bl.outbrk", "bl.YSO", "bl.YSO2", "bl.pdXYSO", "bl.TWIP",
                "a0", "a.Time", "a.Time2", "a.DOY", "a.DOY2", "a.ccov", "a.shcov", "b",
                "lambda", "pcap", # Needed for WAIC
                "test") # GOF

inits <- function() # Setting these based on posterior distribution from an initial successful run.
  list(N = Y, beta0.mean = rnorm(1, 1, 0.1), beta0.sd = rnorm(1, 0.66, 0.07), a0 = rnorm(1, 3.5, 0.5), b = rnorm(1, 3.16, 0.5))

# MCMC values.  Change to what works well for you.
nc <- 3
nb <- 5000
ni <- 10000
nt <- 10

save.out <- "mod_RESQ_allcovs_HZdist_LP"
#############

# Detection data #
Y <- eval(as.name(str_c("Y.", stratum, ".dist")))
dclass <- eval(as.name(str_c("dclass.", stratum)))
dimnames(dclass) <- NULL
gridID <- eval(as.name(str_c("Cov.", stratum)))[, "gridIndex"]
nGrid <- max(gridID)
nPoint <- length(Y)
nInd <- nrow(dclass)

# Covariates (save for later) #
Cov <- eval(as.name(str_c("Cov.", stratum)))
PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > 12)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead outside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & is.na(Cov[, "YSO"]))])
PctDead.b[which(is.na(PctDead.b) & !is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead inside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & !is.na(Cov[, "YSO"]))])
PctDead.sd <- sd(PctDead.b[which(is.na(Cov[, "YSO"]))]) %>% rep(length(PctDead.b)) # SDs for imputing missing values outside ADS polygons
PctDead.sd[which(!is.na(Cov[, "YSO"]))] <- sd(PctDead.b[which(!is.na(Cov[, "YSO"]))]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed

YSO.b <- Cov[, "YSO"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.b[is.na(YSO.b)] <- 0
YSO.d[is.na(YSO.d)] <- 0

Outbrk.b <- Cov[, "YSO"] %>% (function(x) as.integer(!is.na(x)))
Outbrk.d <- tapply(Outbrk.b, gridID, max)

TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- tapply(TWIP.b, gridID, mean, na.rm = T) # Grid-level values
TWIP.sd <- tapply(TWIP.b, gridID, sd, na.rm = T) # Grid-specific SDs for imputing missing values
if(stratum == "SF") { # For the two entire grids missing TWIP in SF data
  ind.fill <- which(is.na(TWIP.d))
  TWIP.d[ind.fill] <- mean(TWIP.d, na.rm = T)
  TWIP.sd[ind.fill] <- sd(TWIP.b, na.rm = T)
}
TWIP.b.missing <- is.na(TWIP.b) %>% as.integer # Index missing values to be imputed
TWIP.b[is.na(TWIP.b)] <- 0

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
if(stratum == "LP") Time.b[point.list.LP == "LP-057-06"] <- 
  mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
  # Imputed missing time value in LP data based on survey times for point likely completed just before and just after.
  # Inferred which points after viewing the other points and their associated times in ArcMap.

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.means <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- 0

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.means <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- 0

st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

max(out$summary[,"Rhat"])

min(out$summary[,"n.eff"])
#sort(out$summary[,"n.eff"])[1:200]

sum(out$sims.list$test) / out$mcmc.info$n.samples

library(R.utils)
saveObject(out, save.out)
